library(readxl)
library(msigdbr)
library(ggplot2)
library(edgeR)
library(dplyr)
library(sva)
library(viridis)
library(reshape2)
library(tidyr)

meta_filtered <- read_excel("filtered_metadata.xlsx") #metadata
raw_filtered <- read.table("raw_filtered.txt",sep="\t",header=TRUE) #raw count data

#data cleaning and normalization
cpm <- cpm(raw_filtered[,2:3061])
libsizes <- calcNormFactors(raw_filtered[,2:3061])
tpm <- cpm / libsizes
log_tpm <-log2(tpm+1)
rownames(log_tpm) <- raw_filtered$X
combat_data <- ComBat(log_tpm, batch = meta_filtered$Batch)

#subsetting
health_meta <- subset(meta_filtered)
age_sort <- health_meta[order(health_meta$Age),]
raw_sort <- combat_data[,age_sort$SRR.ID]

#--------------------------------------------------------------------------------------------
#make significant matrix

#retrieve gene sets 
all_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")

srr_info <- data.frame(
  SRR.ID = colnames(raw_sort),
  Condition = meta_filtered$Condition[match(colnames(raw_sort), meta_filtered$SRR.ID)]
)

healthy_indices <- which(meta_filtered$Condition == "Healthy")
disease_indices <- which(meta_filtered$Condition != "Healthy")

#differentially expressed genes
gene_set_expression <- function(gene_set, expression_data, indices) {
  genes_in_set <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
  valid_genes <- genes_in_set[genes_in_set %in% rownames(expression_data)]
  
  if (length(valid_genes) == 0) {
    return(rep(NA, length(indices)))  # Return NAs if no genes are found
  }
  
  expression_subset <- expression_data[valid_genes, indices, drop = FALSE]
  if (ncol(expression_subset) == 0) {
    return(rep(NA, length(indices)))  # Return NAs if no data columns
  }
  
  return(apply(expression_subset, 2, function(col) {
    if (all(is.na(col))) NA else mean(col, na.rm = TRUE)
  }))
}

deg_gene_sets <- list()
for (gene_set in unique(all_gene_sets$gs_name)) {
  expression_healthy <- gene_set_expression(gene_set, raw_sort, healthy_indices)
  expression_disease <- gene_set_expression(gene_set, raw_sort, disease_indices) 
  
  #check if both groups have non-null data and sufficient non-NA samples for a valid t-test
  if (!is.null(expression_healthy) && sum(!is.na(expression_healthy)) >= 2 &&
      !is.null(expression_disease) && sum(!is.na(expression_disease)) >= 2) {
    tryCatch({
      t_test_result <- t.test(expression_disease, expression_healthy, na.action = na.omit)
      
      #store results
      deg_gene_sets[[gene_set]] <- list(
        p_value = t_test_result$p.value,
        estimate = mean(expression_disease, na.rm = TRUE) - mean(expression_healthy, na.rm = TRUE)
      )
    }, error = function(e) {
      #log error when t-test fails
      deg_gene_sets[[gene_set]] <- list(
        p_value = NA,
        estimate = NA,
        note = paste("Error in t-test:", e$message)
      )
    })
  } else {
    #log or store information about gene sets with insufficient non-NA data
    deg_gene_sets[[gene_set]] <- list(
      p_value = NA,
      estimate = NA,
      note = "Insufficient non-NA data for t-test"
    )
  }
}
#print(deg_gene_sets)

#convert results list to dataframe
deg_gs_df <- do.call(rbind, lapply(names(deg_gene_sets), function(gene_set) {
  estimate = deg_gene_sets[[gene_set]]$estimate
  if (is.na(estimate) || estimate + 1 <= 0) {
    log2fold_change <- NA  #set as NA if the value is not suitable for log2(x+1)
  } else {
    log2fold_change <- log2(estimate + 1)
  }
  data.frame(
    GeneSet = gene_set,
    P_Value = deg_gene_sets[[gene_set]]$p_value,
    Log2FoldChange = log2fold_change
  )
}))
#print(deg_gs_df)
#output results
#print(significant_gene_sets)

#calculate the Bonferroni adjusted p-values
deg_gs_df$P_Bonferroni <- p.adjust(deg_gs_df$P_Value, method = "bonferroni")
#define new thresholds based on Bonferroni adjusted p-values
gs_bon_threshold = 0.05 / nrow(deg_gs_df)  
deg_gs_df$Significant <- ifelse(-log10(deg_gs_df$P_Bonferroni) > -log10(gs_bon_threshold), "Significant", "No")

#volcano plot disease independent
volcano_plot_bonferroni <- ggplot(deg_gs_df, aes(x = Log2FoldChange, y = -log10(P_Bonferroni), color = Significant)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("No" = "black", "Significant" = viridis(6)[3])) +
  geom_hline(yintercept = -log10(gs_bon_threshold), linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(legend.position = "none") +  # Correct placement of the '+'
  labs(x = "Log Fold Change Per Year", y = "-log10(p-value)", title = "Summary of Disease-Independent Differentially Expressed Gene Sets") +
  xlim(c(-2, 2))
volcano_plot_bonferroni
#ggsave(filename = "Disease_Gene_Sets.png",plot = last_plot(),dpi=1200)

#--------------------------------------------------------------------------------------------
#make significant matrix

#read files
gs_pcc_matrix <- read.csv("healthy_tissue_independent_pcc_matrix.csv", row.names = 1) 
gene_sets <- rownames(gs_pcc_matrix)
unique_conditions <- unique(meta_filtered$Condition[health_meta$Condition != "Healthy"])

#initialize matrices
gs_diff_exp_matrix <- matrix(NA, nrow = length(gene_sets), ncol = length(unique_conditions),
                          dimnames = list(gene_sets, unique_conditions))
gs_p_values_matrix <- matrix(NA, nrow = length(gene_sets), ncol = length(unique_conditions),
                          dimnames = list(gene_sets, unique_conditions))

#calculate differential expression for each disease vs healthy
for (gene_set in gene_sets) {
  healthy_expression <- gene_set_expression(gene_set, raw_sort, healthy_indices)
  
  for (disease in unique_conditions) {
    disease_indices <- which(meta_filtered$Condition == disease)
    disease_expression <- gene_set_expression(gene_set, raw_sort, disease_indices)
    
    #perform t-test if data is adequate
    if (length(na.omit(healthy_expression)) > 2 && length(na.omit(disease_expression)) > 2) {
      t_test_result <- t.test(disease_expression, healthy_expression, alternative = "two.sided")
      gs_diff_exp_matrix[gene_set, disease] <- mean(disease_expression, na.rm = TRUE) - mean(healthy_expression, na.rm = TRUE)
      gs_p_values_matrix[gene_set, disease] <- t_test_result$p.value
    } else {
      gs_diff_exp_matrix[gene_set, disease] <- NA
      gs_p_values_matrix[gene_set, disease] <- NA
    }
  }
}

#remove columns that are entirely NAs
gs_diff_exp_matrix <- gs_diff_exp_matrix[, !apply(is.na(gs_diff_exp_matrix), 2, all)]
gs_p_values_matrix <- gs_p_values_matrix[, !apply(is.na(gs_p_values_matrix), 2, all)]
#write.csv(diff_exp_matrix, file = "gs_diff_exp_matrix.csv", row.names = TRUE)
#write.csv(p_values_matrix, file = "gs_p_values_matrix.csv", row.names = TRUE)

#apply Bonferroni correction
bonferroni_threshold <- 0.05 / (length(gene_sets) * length(unique_conditions))
gs_significant_matrix <- ifelse(gs_p_values_matrix < bonferroni_threshold, 1, 0)
gs_significant_matrix[gs_diff_exp_matrix < 0 & gs_significant_matrix == 1] <- -1
#write.csv(significant_matrix, file = "gs_significant_matrix.csv", row.names = TRUE)

#--------------------------------------------------------------------------------------------
#bar plot

#long data frame
gs_long_df <- melt(gs_significant_matrix)
#naming columns 
colnames(gs_long_df) <- c("GeneSet", "Disease", "Regulation")
#number of upregulated and downregulated genes for each disease
gs_agg_data <- dcast(gs_long_df, Disease ~ Regulation, fun.aggregate = length)
colnames(gs_agg_data) <- c("Disease", "NonSignificant", "DownRegulated", "UpRegulated")
#gathering data to long format for ggplot2
gs_bar_data <- gather(gs_agg_data, RegulationType, Count, UpRegulated, DownRegulated)

#side-by-side bar plot
gs_bar <- ggplot(gs_bar_data, aes(x = Disease, y = Count, fill = RegulationType)) +
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("UpRegulated" = "#73d055ff", "DownRegulated" = "#481567ff")) +
  labs(title = "Differentially Expressed Gene Sets by Disease",
       x = "",
       y = "Number of Gene Sets",
       fill = "Regulation Type") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 14)) # Increase font size here
gs_bar
#ggsave(filename = "gs_disease_bar.png", plot = gs_bar, width = 12, height = 8, dpi = 1200)






