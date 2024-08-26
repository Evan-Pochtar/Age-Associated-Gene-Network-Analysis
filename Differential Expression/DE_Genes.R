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

#read file
g_pcc_matrix <- read.csv("gene_pcc_matrix.csv", row.names = 1)
gene_names <- rownames(g_pcc_matrix)
unique_conditions <- unique(meta_filtered$Condition[health_meta$Condition != "Healthy"])

#initialize matrices
g_diff_exp_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(unique_conditions),
                          dimnames = list(gene_names, unique_conditions))
g_p_values_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(unique_conditions),
                          dimnames = list(gene_names, unique_conditions))

healthy_indices <- which(meta_filtered$Condition == "Healthy")
disease_indices <- which(meta_filtered$Condition != "Healthy")

#differentially expressed genes
gene_expression <- function(gene, raw_sort, indices) {
  # Get the row index for the gene
  gene_index <- which(rownames(raw_sort) == gene)
  # Check if the gene is present in the data
  if (length(gene_index) == 0) {
    cat("Gene not found:", gene, "\n")
    return(NULL)
  }
  # Subset the expression values using the provided indices
  g_expression_values <- raw_sort[gene_index, indices]
  return(g_expression_values)
}

#calculate differential expression for each disease vs healthy
for (gene in gene_names) {
  g_healthy_expression <- gene_expression(gene, raw_sort, healthy_indices)
  
  for (disease in unique_conditions) {
    disease_indices <- which(meta_filtered$Condition == disease)
    g_disease_expression <- gene_expression(gene, raw_sort, disease_indices)
    
    #perform t-test if data is adequate
    if (length(na.omit(g_healthy_expression)) > 2 && length(na.omit(g_disease_expression)) > 2) {
      t_test_result <- t.test(g_disease_expression, g_healthy_expression, alternative = "two.sided")
      g_diff_exp_matrix[gene, disease] <- mean(g_disease_expression, na.rm = TRUE) - mean(g_healthy_expression, na.rm = TRUE)
      g_p_values_matrix[gene, disease] <- t_test_result$p.value
    } else {
      g_diff_exp_matrix[gene, disease] <- NA
      g_p_values_matrix[gene, disease] <- NA
    }
  }
}

#remove columns that are entirely NAs
g_diff_exp_matrix <- g_diff_exp_matrix[, !apply(is.na(g_diff_exp_matrix), 2, all)]
g_p_values_matrix <- g_p_values_matrix[, !apply(is.na(g_p_values_matrix), 2, all)]
#write.csv(diff_exp_matrix, file = "gene_diff_exp_matrix.csv", row.names = TRUE)
#write.csv(p_values_matrix, file = "gene_p_values_matrix.csv", row.names = TRUE)

#apply Bonferroni correction
bonferroni_threshold <- 0.05 / (length(gene_names) * length(unique_conditions))
g_significant_matrix <- ifelse(g_p_values_matrix < bonferroni_threshold, 1, 0)
g_significant_matrix[g_diff_exp_matrix < 0 & g_significant_matrix == 1] <- -1
#write.csv(significant_matrix, file = "gene_significant_matrix.csv", row.names = TRUE)

#--------------------------------------------------------------------------------------------
#bar plot

#long data frame
g_long_df <- melt(g_significant_matrix)
#naming columns 
colnames(g_long_df) <- c("Genet", "Disease", "Regulation")
#number of upregulated and downregulated genes for each disease
g_agg_data <- dcast(g_long_df, Disease ~ Regulation, fun.aggregate = length)
colnames(g_agg_data) <- c("Disease", "NonSignificant", "DownRegulated", "UpRegulated")
#gathering data to long format for ggplot2
g_bar_data <- gather(g_agg_data, RegulationType, Count, UpRegulated, DownRegulated)

#side-by-side bar plot
g_bar <- ggplot(g_bar_data, aes(x = Disease, y = Count, fill = RegulationType)) +
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.7) +  
  scale_fill_manual(values = c("UpRegulated" = "#73d055ff", "DownRegulated" = "#481567ff")) +
  labs(title = "Differentially Expressed Gene by Disease",
       x = "",
       y = "Number of Genes",
       fill = "Regulation Type") +
  theme_minimal() +
  coord_flip() +  
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 14)) 
g_bar
#ggsave(filename = "g_disease_bar.png", width = 12, height = 8, dpi=1200)




