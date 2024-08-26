library(readxl)
library(ggplot2)
library(ReactomePA)
library(reactome.db)
library(AnnotationDbi)
library(ReactomeContentService4R)
library(tidyr)
library(ggrepel)
library(edgeR)
library(dplyr)
library(sva)
library(DESeq2)
library(msigdbr)
library(tidyverse)
library(philentropy)
library(gridExtra)
library(car)
library(viridis)

# Import data (may need to change your working directory)
meta_filtered <- read_excel("filtered_metadata.xlsx") # metadata
raw_filtered <- read.table("raw_filtered.txt",sep="\t",header=TRUE) # raw count data


# Other potential pipeline set-ups, feel free to ignore
# merged_data <- inner_join(meta_filtered, raw_filtered, by = c("SRR.ID" = "X"))
# counts_matrix <- merged_data[, -c("SRR.ID", "Batch")]
# dge <- DGEList(raw_filtered)
# dge <- calcNormFactors(dge)
# pseudo_count <- 1
# log2_counts <- log2(dge$counts + pseudo_count)
# rownames(log2_counts) <- raw_filtered$X
# merged_data <- merge(log2_counts, age_sort_meta, by.x = "SRR.ID", by.y = "SRR.ID", all.x = TRUE)

# dds <- DESeqDataSetFromMatrix(countData = raw_filtered[,2:3061],
#                               colData = meta_filtered,
#                               design = ~ Condition + Tissue)

# Data cleaning and normalization
cpm <- cpm(raw_filtered[,2:3061])
libsizes <- calcNormFactors(raw_filtered[,2:3061])
tpm <- cpm / libsizes

log_tpm <-log2(tpm+1)
rownames(log_tpm) <- raw_filtered$X
combat_data <- ComBat(log_tpm, batch = meta_filtered$Batch)


# combat_first <- ComBat(raw_filtered[,2:length(colnames(raw_filtered))],batch = meta_filtered$Batch)
# cpm_second <- cpm(combat_first)
# libsizes_second <- calcNormFactors(cpm_second)
# tpm_second <- cpm_second/libsizes_second
# log_tpm_second <- log2(tpm_second+1)
# other_combat <- ComBat(tpm,batch=meta_filtered$Batch)
# min(combat_first)
# min(raw_filtered[,2:3061])

# Data subsetting
healthy_meta <- subset(meta_filtered, Healthy == TRUE)
non_healthy_meta <- subset(meta_filtered,Healthy == FALSE)
age_sort_non_healthy <- non_healthy_meta[order(non_healthy_meta$Age),]
age_sort_meta <- healthy_meta[order(healthy_meta$Age),]
raw_sort <- combat_data[,age_sort_meta$SRR.ID]
non_healthy_raw_sort <- combat_data[,age_sort_non_healthy$SRR.ID]

write.table(age_sort_meta,"age_sort_meta.tsv",sep="\t",row.names=FALSE)
age_sort_meta_t <- read_tsv("age_sort_meta.tsv")

median_vector <- c()
for (i in 1:length(age_sort_meta$Age)) {
  median_vector <- c(median_vector,median(age_sort_meta$Age[i:length(age_sort_meta$Age)]))
}
minimum_median <- which(median_vector == median(age_sort_non_healthy$Age))[1]
median_adjusted_healthy <- age_sort_meta[minimum_median:length(age_sort_meta$Age),]
median(non_healthy_meta$Age)
median(median_adjusted_healthy$Age)

# Proof of concept for gene set enrichment (by gene)
pathway_ids <- c("R-HSA-73884", "R-HSA-73893", "R-HSA-73942", "R-HSA-5693606",
                 "R-HSA-6783310", "R-HSA-5685942", "R-HSA-5693567", "R-HSA-5685939",
                 "R-HSA-5685938", "R-HSA-5358508", "R-HSA-5696398", "R-HSA-5693571",
                 "R-HSA-69618", "R-HSA-2484822")
pathway_test <- getParticipants("R-HSA-73884", retrieval = "EventsInPathways")
pathway_1 <- event2Ids(event.id = "R-HSA-73884")
pathway_1$geneSymbol

# Generate a linear regression model for change in log2(Expression) with age
if (exists("sig_table")) {
  rm(sig_table)
}
for (id in pathway_ids) {
  set_id <- event2Ids(event.id = id) # pull out pathway id
  relevant <- intersect(set_id$geneSymbol, raw_filtered$X) # subset count data by genes in that set 
  gene_data <- raw_sort[relevant, ]
  gene_data <- as.matrix(gene_data)
  gene_data_df <- as.data.frame(gene_data)
  gene_data_df$Gene <- rownames(gene_data_df)
  gene_data_long <- gather(gene_data_df, key = "SRR.ID", value = "Expression", -Gene) # arrange data into an optimal format for regression analysis
  gene_data_long_merged <- merge(gene_data_long, age_sort_meta, by.x = "SRR.ID", by.y = "SRR.ID", all.x = TRUE)
  regression_results <- by(gene_data_long_merged, list(gene_data_long_merged$Gene, gene_data_long_merged$Tissue), function(df) { # set up linear regression function
    model <- lm(Expression ~ Age, data = df) # set regression parameters
    coef_summary <- summary(model)$coefficients
    if (nrow(coef_summary) > 1) { # screen for proper outputs
      slope <- coef_summary[2, 1]
      p_value <- coef_summary[2, 4]
    } else {
      slope <- NA
      p_value <- NA
    }
    gene <- unique(df$Gene) # capture relevant metafata components
    tissue <- unique(df$Tissue)
    age <- unique(df$Age)
    return(list(Slope = slope, P_Value = p_value, Gene = gene, Tissue = tissue, Age = age))
  })
  
  # Capture by-gene regression results in a dataframe, then perform statistical analysis
  regression_summary <- do.call(rbind, regression_results)
  regression_summary <- as.data.frame(regression_summary)
  regression_summary$Log_P <- -log10(as.numeric(regression_summary$P_Value))
  regression_sorted <- regression_summary[order(as.numeric(regression_summary$P_Value)), ]
  regression_sorted$FDR <- p.adjust(regression_sorted$P_Value, method = "BH")
  regression_sorted <- regression_sorted[!is.na(regression_sorted$FDR), ]
  regression_sorted$Log_FDR <- -log10(as.numeric(regression_sorted$FDR))
  regression_sorted$ID <- id
  significant <- regression_sorted[regression_sorted$FDR < 0.05,]
  if (exists("sig_table")) { # identify significant genes
    sig_table <- rbind(sig_table, significant)
  } else {
    sig_table <- data.frame(significant)
  }
  
  # Plot regression data for each pathway
  print(ggplot(regression_sorted, aes(x = as.numeric(Slope), y = Log_FDR, label = ifelse(FDR < 0.05, Gene, ""))) +
          geom_point(size = 1, color = ifelse(regression_sorted$FDR < 0.05, "red", "black")) +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
          theme_minimal() +
          labs(x = "Log Fold Change per Year", y = "-log10(FDR)", title = paste("Volcano Plot of",id)) +
          xlim(-0.25, 0.25) +
          geom_text_repel(
            box.padding = 0.2,
            point.padding = 0.2,
            segment.color = NA,
            segment.size = 0,
            segment.alpha = 0,
            force = 0.2,
            nudge_x = 0,
            min.segment.length = 0.1
          )
  )
}
ggsave(filename = "Gene_set_example.png",plot = last_plot(),dpi=1200)
rownames(sig_table) <- NULL

# Proof of concept: it is possible to pull out significant pathways from by-gene data
gene_sorted_sig_table <- sig_table[order(as.character(sig_table$Gene)), ]
sig_pathways <- table(sig_table$ID)

# Import curated gene sets from human MSigDB
all_gene_sets <- msigdbr(species = "Homo sapiens",category = "C2") # might be worth splitting C2 into CGP and CP if there's time
write.csv(all_gene_sets,"all_gene_sets.csv",row.names=FALSE)
all_gene_sets_t <- read.csv("all_gene_sets.csv",sep=",")


# Set up function to perform linear regression analysis on gene sets (default option is tissue-independent, but can also set to tissue-specific)
perform_regression <- function(gene_set, gene_expression_data, age_tissue_data, tissue_specific = FALSE) {
  # Extract genes in the current gene set
  gene_symbols <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
  
  # Filter gene expression data for genes in the current gene set
  gene_set_expression <- gene_expression_data[rownames(gene_expression_data) %in% gene_symbols, ]
  
  result <- tryCatch({
    if (!is.null(gene_set_expression) && nrow(gene_set_expression) > 0 && ncol(gene_set_expression) > 0) { # if the gene set can be pulled out from the data
      # Aggregate gene expression data across genes within the gene set
      gene_set_expression_aggregated <- apply(gene_set_expression, 2, mean, na.rm = TRUE)
      
      if (tissue_specific) {
        # Merge with age and tissue metadata for tissue-specific regression
        merged_data <- merge(data.frame(Expression = gene_set_expression_aggregated, "SRR.ID" = names(gene_set_expression_aggregated)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
        
        # Perform linear regression within each tissue
        regression_results <- by(merged_data, merged_data$Tissue, function(df) {
          if (nrow(df) > 1) {
            model <- lm(Expression ~ Age, data = df)
            coefficients <- coef(summary(model))
            return(coefficients)
          } else {
            return(NULL)
          }
        })
      } else {
        # Merge with age metadata for tissue-independent regression
        merged_data <- merge(data.frame(Expression = gene_set_expression_aggregated, "SRR.ID" = names(gene_set_expression_aggregated)), age_tissue_data[, c("SRR.ID", "Age")], by = "SRR.ID", all.x = TRUE)
        
        # Perform tissue-independent linear regression
        model <- lm(Expression ~ Age, data = merged_data)
        regression_results <- list(all = coef(summary(model)))
      }
      
      return(list(gene_set = gene_set, regression_results = regression_results))
    } else {  # If gene_set_expression is empty or NULL, return NULL
      message("Gene set expression data is empty or NULL for gene set: ", gene_set)
      return(NULL)
    }
  }, error = function(e) {  # Catch any errors that occur
    message("Error occurred: ", e$message)
    return(NULL)
  })
  
  return(result)
}
# Note: some errors do occur, likely due to either empty sets or gene id discrepancies,
# but they account for less than 1% of gene sets, and the function is set up to handle them.

# Perform regression analysis for each gene set (tissue-independent)
regression_results_tissue_independent <- lapply(unique(all_gene_sets$gs_name), function(gene_set) {
  perform_regression(gene_set, raw_sort, age_sort_meta)
})

# Perform regression analysis for each gene set (tissue-specific)
regression_results_tissue_specific <- lapply(unique(all_gene_sets$gs_name), function(gene_set) {
  perform_regression(gene_set, raw_sort, age_sort_meta, tissue_specific = TRUE)
})

perform_regression_gene <- function(gene_expression_data, age_tissue_data, tissue_specific = FALSE) {
  # Initialize a list to store regression results for each gene
  regression_results <- list()
  
  # Loop over genes
  for (gene_symbol in rownames(gene_expression_data)) {
    # Filter gene expression data for the current gene
    gene_expression <- gene_expression_data[gene_symbol,]
    
    # If gene expression data is available
    if (!is.null(gene_expression)) {
      # If tissue-specific regression is requested
      if (tissue_specific) {
        # Merge with age and tissue metadata for tissue-specific regression
        merged_data <- merge(data.frame(Expression = gene_expression, SRR.ID = colnames(gene_expression_data)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
        
        # Perform linear regression within each tissue
        regression_results[[gene_symbol]] <- by(merged_data, merged_data$Tissue, function(df) {
          if (nrow(df) > 1) {
            model <- lm(Expression ~ Age, data = df)
            coefficients <- coef(summary(model))
            return(coefficients)
          } else {
            return(NULL)
          }
        })
      } else {  # If tissue-specific regression is not requested
        # Merge with age metadata for tissue-independent regression
        merged_data <- merge(data.frame(Expression = gene_expression, SRR.ID = colnames(gene_expression_data)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
        
        # Perform tissue-independent linear regression
        model <- lm(Expression ~ Age, data = merged_data)
        regression_results[[gene_symbol]] <- list(all = coef(summary(model)))
      }
    } else {  # If gene expression data is empty or NULL
      message("Gene expression data is empty or NULL for gene: ", gene_symbol)
      regression_results[[gene_symbol]] <- NULL
    }
  }
  
  return(regression_results)
}

# Perform regression analysis for each gene (tissue-specific)
gene_regression_results_tissue_specific <- perform_regression_gene(raw_sort, age_sort_meta, tissue_specific = TRUE)
simple_gene_regression_results_tissue_independent <- perform_regression_gene(raw_sort, age_sort_meta, tissue_specific = FALSE)


# Reformat the tissue_independent regression data for downstream analysis
regression_summary_independent <- do.call(rbind,regression_results_tissue_independent)
regression_summary_independent <- as.data.frame(regression_summary_independent)
for (i in 1:length(regression_summary_independent$gene_set)) {
  regression_summary_independent$intercept_estimate[i] <- regression_summary_independent$regression_results[[i]][[1]][1,1]
  regression_summary_independent$intercept_std_error[i] <- regression_summary_independent$regression_results[[i]][[1]][1,2]
  regression_summary_independent$intercept_t_value[i] <- regression_summary_independent$regression_results[[i]][[1]][1,3]
  regression_summary_independent$intercept_significance[i] <- regression_summary_independent$regression_results[[i]][[1]][1,4]
  regression_summary_independent$coefficient_estimate[i] <- regression_summary_independent$regression_results[[i]][[1]][2,1]
  regression_summary_independent$coefficient_std_error[i] <- regression_summary_independent$regression_results[[i]][[1]][2,2]
  regression_summary_independent$coefficient_t_value[i] <- regression_summary_independent$regression_results[[i]][[1]][2,3]
  regression_summary_independent$coefficient_significance[i] <- regression_summary_independent$regression_results[[i]][[1]][2,4]
}

# Statistical analysis
regression_summary_independent$log_coefficient_significance <- -log10(regression_summary_independent$coefficient_significance)
bon_threshold <- 0.05/length(regression_summary_independent$gene_set)

# Exporting significant gene sets for JSD analysis via MSI
regression_summary_export <- regression_summary_independent[, -2]
regression_summary_export <- as.data.frame(lapply(regression_summary_export, unlist))

write.csv(regression_summary_export, "regression_summary_independent.csv", row.names = FALSE)
regression_summary_independent_t <- read.csv("regression_summary_independent.csv",sep=",")

sig_gene_sets <- regression_summary_export$gene_set[abs(regression_summary_independent$coefficient_estimate) > 0.02]


# Graphing gene sets by log2(expression) change with age
ggplot(regression_summary_independent, aes(x = coefficient_estimate, y = log_coefficient_significance)) +
  geom_point(size = 1, color = ifelse(regression_summary_independent$log_coefficient_significance > log10(0.05/bon_threshold) & abs(regression_summary_independent$coefficient_estimate) > 0.02, viridis(3)[2], "black")) +
  geom_hline(yintercept = -log10(bon_threshold), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.02, 0.02), linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Log Fold Change per Year", y = "-log10(p-value)", title = paste("Summary of Tissue-Independent Gene Set Effects")) +
  xlim(-0.1, 0.1)
ggsave(filename = "Tissue_Independent_Gene_Set_Effects.png",plot = last_plot(),dpi=1200)

length(regression_summary_independent$gene_set[regression_summary_independent$coefficient_significance < bon_threshold])

# Reformat tissue-specific regression outputs into interpretable data
regression_summary_by_tissue <- do.call(rbind,regression_results_tissue_specific)
regression_summary_by_tissue <- as.data.frame(regression_summary_by_tissue)
for (i in 1:length(regression_summary_by_tissue$gene_set)) {
  for (j in 1:length(unique(age_sort_meta$Tissue))) {
    if (!is.null(regression_summary_by_tissue$regression_results[[i]][[j]])) {
      for (k in 1:2) {
        for (l in 1:4) {
          if (nrow(regression_summary_by_tissue$regression_results[[i]][[j]]) > 1) {
            column_name <- paste(rownames(regression_summary_by_tissue$regression_results[[i]][[j]])[k],
                              colnames(regression_summary_by_tissue$regression_results[[i]][[j]])[l],
                              names(regression_summary_by_tissue$regression_results[[i]][j]),
                              sep = "_")
            regression_summary_by_tissue[[column_name]][i] <- regression_summary_by_tissue$regression_results[[i]][[j]][k,l]
          } else {
            column_name <- paste(rownames(regression_summary_by_tissue$regression_results[[1]][[1]])[k],
                               colnames(regression_summary_by_tissue$regression_results[[1]][[1]])[l],
                               names(regression_summary_by_tissue$regression_results[[i]][j]),
                               sep = "_")
            if (k==1) {
              regression_summary_by_tissue[[column_name]][i] <- regression_summary_by_tissue$regression_results[[i]][[j]][k,l]
            } else {
              regression_summary_by_tissue[[column_name]][i] <- NA
            }
          }
        }
      }
    } else {
      for (k in 1:2) {
        for (l in 1:4) {
          column_name <- paste(rownames(regression_summary_by_tissue$regression_results[[1]][[1]])[k],
                             colnames(regression_summary_by_tissue$regression_results[[1]][[1]])[l],
                             names(regression_summary_by_tissue$regression_results[[i]][j]),
                             sep = "_")
          regression_summary_by_tissue[[column_name]][i] <- NA
        }
      }
    }
  }
}

# Reformat tissue-specific regression data for downstream analysis
tissue_specific_names <- grep("Age_Estimate", names(regression_summary_by_tissue), value = TRUE)
tissue_specific_changes <- regression_summary_by_tissue[, c("gene_set",tissue_specific_names)]
tissue_names <- sub("Age_Estimate_", "", names(tissue_specific_changes)[-1])
tissue_names_clean <- sub("Age_Estimate_", "", names(tissue_specific_changes)[-1])
tissue_specific_changes_long <- pivot_longer(tissue_specific_changes, 
                        cols = -gene_set, 
                        names_to = "tissue", 
                        values_to = "regression_value")
# Note there's no regression data for osteoblasts because there were only two patients
# and they were of the same age, making regression analysis impossible.
tissue_specific_p_value_names <- grep("Age_Pr(>|t|)", names(regression_summary_by_tissue), value = TRUE)
tissue_specific_p_values <- regression_summary_by_tissue[, c("gene_set",tissue_specific_p_value_names)]
changes_columns_to_remove <- c("Age_Estimate_Bone; osteoblasts", "Age_Estimate_Cerebrospinal fluid; exosomes")
p_value_columns_to_remove <- c("Age_Pr(>|t|)_Bone; osteoblasts", "Age_Pr(>|t|)_Cerebrospinal fluid; exosomes")
tissue_specific_changes_clean <- tissue_specific_changes[, !colnames(tissue_specific_changes) %in% changes_columns_to_remove]
tissue_specific_p_value_clean <- tissue_specific_p_values[, !colnames(tissue_specific_p_values) %in% p_value_columns_to_remove]
tissue_names_clean <- sub("Age_Estimate_", "", names(tissue_specific_changes_clean)[-1])

volcano_plots <- list()

# Loop over columns and create volcano plots
for (i in 2:ncol(tissue_specific_changes_clean)) {
  # Get column names
  col_name_changes <- colnames(tissue_specific_changes_clean)[i]
  col_name_p_value <- colnames(tissue_specific_p_value_clean)[i]
  
  # Create volcano plot
  volcano_plot <- ggplot(data = NULL, aes(x = tissue_specific_changes_clean[, i], y = -log10(tissue_specific_p_value_clean[, i]))) +
    geom_point() +
    labs(x = paste("Coefficient -", col_name_changes),
         y = paste("-log10(p-value) -", col_name_p_value),
         title = paste("Volcano Plot for", col_name_changes)) +
    theme_minimal()
  print(volcano_plot)
  # Add the plot to the list
  volcano_plots[[i]] <- volcano_plot
}
grid.arrange(grobs = volcano_plots, ncol = 3)

# Attempt to assess significant gene sets by tissue
tissue_bon_threshold <- 0.05 / (nrow(tissue_specific_changes_clean) * ncol(tissue_specific_changes_clean))
sig_gene_sets_by_tissue <- list()

for (i in 2:length(tissue_names_clean)) {
  tissue <- tissue_names_clean[i]
  tissue_range <- range(meta_filtered$Age[which(meta_filtered$Tissue == tissue)])
  coefficient_threshold <- 1 / tissue_range
  sig_sets_for_tissue <- tissue_specific_changes_clean$gene_set[abs(tissue_specific_changes_clean[, i]) > coefficient_threshold & tissue_specific_p_value_clean[, i] < tissue_bon_threshold]
  
  if (length(sig_sets_for_tissue) > 0) {
    sig_gene_sets_by_tissue[[tissue]] <- sig_sets_for_tissue
  } else {
    sig_gene_sets_by_tissue[[tissue]] <- NA
  }
}

# Identify tissue-independent age-associated gene sets
sig_gene_sets <- unlist(regression_summary_independent$gene_set[abs(regression_summary_independent$coefficient_estimate) > 0.02 & regression_summary_independent$coefficient_significance < bon_threshold])

sig_gene_set_expression <- do.call(rbind, lapply(sig_gene_sets, function(gene_set) {
  gene_symbols <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
  gene_set_expression <- raw_sort[rownames(raw_sort) %in% gene_symbols, ]
  gene_set_means <- data.frame(t(colMeans(gene_set_expression)))
  row.names(gene_set_means) <- gene_set
  return(gene_set_means)
}))

# Compute Pearson Correlation coefficients
pearson_correlation <- function(x, y) {
  cor(x, y)
}
pcc_matrix <- apply(sig_gene_set_expression, 1, function(x) {
  apply(sig_gene_set_expression, 1, function(y) {
    pearson_correlation(x, y)
  })
})
write.csv(pcc_matrix, "healthy_tissue_independent_pcc_matrix.csv", row.names = TRUE)

# More by-tissue analysis (not included in final paper)
gene_regression_summary_by_tissue <- do.call(rbind,gene_regression_results_tissue_specific)
gene_regression_summary_by_tissue <- as.data.frame(gene_regression_summary_by_tissue)

gene_slope_by_tissue <- matrix(NA,nrow=nrow(gene_regression_summary_by_tissue),ncol=length(tissue_names_clean),dimnames = list(rownames(gene_regression_summary_by_tissue),tissue_names_clean))
gene_p_value_by_tissue <- matrix(NA,nrow=nrow(gene_regression_summary_by_tissue),ncol=length(tissue_names_clean),dimnames = list(rownames(gene_regression_summary_by_tissue),tissue_names_clean))
for (gene in rownames(gene_regression_summary_by_tissue)) {
  for (tissue in tissue_names_clean) {
    entry_slope <- gene_regression_summary_by_tissue[[gene,tissue]][2,1]
    entry_p_value <- gene_regression_summary_by_tissue[[gene,tissue]][2,4]
    gene_slope_by_tissue[gene,tissue] <- entry_slope
    gene_p_value_by_tissue[gene,tissue] <- entry_p_value
  }
}


library(matrixStats)  # Required for rowMeans2 function

# Define a function to compute JSD between two probability distributions (not used)
# jsd <- function(p, q) {
#   m <- 0.5 * (p + q)
#   kl_pm <- sum(p * log2(p / m))
#   kl_qm <- sum(q * log2(q / m))
#   sqrt(0.5 * (kl_pm + kl_qm))
# }

mean_jsd <- matrix(NA, nrow = length(sig_gene_sets), ncol = length(sig_gene_sets),
                   dimnames = list(sig_gene_sets, sig_gene_sets))
sd_jsd <- matrix(NA, nrow = length(sig_gene_sets), ncol = length(sig_gene_sets),
                 dimnames = list(sig_gene_sets, sig_gene_sets))

raw_sort_positive <- raw_sort
raw_sort_positive[raw_sort_positive<0] <- 0
raw_sort_positive <- cbind(data.frame(rownames(raw_sort_positive)),raw_sort_positive)
write.csv(raw_sort_positive, "raw_sort_positive.csv", row.names = TRUE)
raw_sort_positive_t <- read.csv("raw_sort_positive.csv",sep=",")


# Compute JSD between pairs of gene sets across age (test script; real one was submitted to MSI)
start_time <- Sys.time()
jsd_values <- lapply(sig_gene_sets[1:3], function(gene_set1) {
  gene_symbols1 <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set1]
  gene_set_expression1 <- raw_sort_positive[rownames(raw_sort_positive) %in% gene_symbols1, ]
  
  lapply(sig_gene_sets[1:3], function(gene_set2) {
    gene_symbols2 <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set2]
    gene_set_expression2 <- raw_sort_positive[rownames(raw_sort_positive) %in% gene_symbols2, ]
    
    jsd_values_across_age <- sapply(colnames(raw_sort), function(age_column_name) {
      # Estimate PDFs for the two continuous distributions
      expr1 <- gene_set_expression1[, age_column_name]
      expr2 <- gene_set_expression2[, age_column_name]
      expr_table <- t(cbind(expr1,expr2))
      suppressMessages({gJSD(expr_table,est.prob = "empirical")})
    })
    mean_jsd[gene_set1, gene_set2] <<- mean(jsd_values_across_age)
    sd_jsd[gene_set1, gene_set2] <<- sd(jsd_values_across_age)
  })
})
end_time <- Sys.time()
elapsed_time <- end_time - start_time
elapsed_time

# Linking linear regression and network outputs to hallmarlks of aging
hoa_genes <- read.delim("/Users/administrator/Downloads/hallmarks_of_aging.gmt", header = FALSE, sep = "\t")
row.names(hoa_genes) <- hoa_genes[,1]

hoa_component <- matrix(0, nrow = length(sig_gene_sets), ncol = length(rownames(hoa_genes)),
                     dimnames = list(sig_gene_sets, rownames(hoa_genes)))

for (gene_set in sig_gene_sets) {
  gene_symbols <- intersect(all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set],rownames(raw_sort))
  for (gene in gene_symbols) {
    row_indices <- which(hoa_genes == gene, arr.ind = TRUE)[, 1]
    row_names <- names(row_indices)
    hoa_component[gene_set,row_names] <- hoa_component[gene_set,row_names] + 1
  }
}

expected_hoa <- rowSums(!is.na(hoa_genes))/(ncol(hoa_genes)-2)
gene_normalized_hoa_component <- hoa_component
for (gene_set in rownames(gene_normalized_hoa_component)) {
  genes <- all_gene_sets$gene_symbol[all_gene_sets$gs_name==gene_set]
  true_genes <- intersect(genes,rownames(raw_sort))
  gene_normalized_hoa_component[gene_set,] <- gene_normalized_hoa_component[gene_set,]/length(true_genes)
}
hoa_enrichment <- colMeans(gene_normalized_hoa_component)
hoa_baseline <- rowSums(!is.na(hoa_genes))/nrow(raw_sort)
summary(aov(hoa_enrichment~hoa_baseline))

set_hoa_df <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"),
                                    baseline = hoa_baseline,
                                    enrichment = hoa_enrichment)

y_limit <- max(hoa_enrichment)
ggplot(set_hoa_df, aes(x = hoa, y = baseline, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "HOA Enrichment",fill = "Hallmarks of Aging", title = "Baseline HOA Prediction") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, y_limit)) +
  theme(axis.text.x = element_blank())
ggsave(filename = "baseline_hoa_prediction.png",plot = last_plot(),dpi=1200)

ggplot(set_hoa_df, aes(x = hoa, y = enrichment, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "HOA Enrichment",fill = "Hallmarks of Aging", title = "Network HOA Prediction") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, y_limit)) +
  theme(axis.text.x = element_blank())
ggsave(filename = "network_hoa_prediction.png",plot = last_plot(),dpi=1200)

gene_hoa_component <- matrix(0, nrow = length(simple_sig_genes), ncol = length(rownames(hoa_genes)),
                        dimnames = list(simple_sig_genes, rownames(hoa_genes)))

for (gene in simple_sig_genes) {
  row_indices <- which(hoa_genes == gene, arr.ind = TRUE)[, 1]
  row_names <- names(row_indices)
  gene_hoa_component[gene,row_names] <- gene_hoa_component[gene,row_names] + 1
}

gene_hoa_enrichment <- colMeans(gene_hoa_component)
gene_hoa_baseline <- rowSums(!is.na(hoa_genes))/nrow(raw_sort)
summary(aov(gene_hoa_enrichment~gene_hoa_baseline))

rowSums(!is.na(hoa_genes))/nrow(raw_sort)
adjusted_hoa_component <- gene_normalized_hoa_component/expected_hoa

write.csv(hoa_component,"hoa_component.csv")
write.csv(gene_normalized_hoa_component, "gene_normalized_hoa_component.csv")
write.csv(adjusted_hoa_component,"adjusted_hoa_component.csv")

adjusted_hoa_100_sums <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"),
                                    values = c(332.434710,9.871978,125.380199,9.747676,23.285416,123.635788,177.125010,39.187062,19.220703))
adjusted_hoa_100_sums$norm_values <- adjusted_hoa_100_sums$values/sum(adjusted_hoa_100_sums$values)
ggplot(adjusted_hoa_100_sums, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "Bar Plot with Viridis Color Palette") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

adjused_hoa_100_edge_comps <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"), 
                                         values = c(50410.876763, 1471.642774, 19819.278091, 1623.170900, 3712.419403, 20091.403441, 27452.751902, 6276.385142, 3047.227362))
adjused_hoa_100_edge_comps$norm_values <- adjused_hoa_100_edge_comps$values/sum(adjused_hoa_100_edge_comps$values)
ggplot(adjused_hoa_100_edge_comps, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "adjused_hoa_100_edge_comps") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave(filename = "adjused_hoa_100_edge_comps.png",plot = last_plot(),dpi=1200)

gn_hoa_100_sums <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"), 
                                         values = c(12.115787, 0.344165, 4.450203, 0.403539, 0.772078, 4.831132, 6.209809, 1.233794, 0.664789))
gn_hoa_100_sums$norm_values <- gn_hoa_100_sums$values/sum(gn_hoa_100_sums$values)
ggplot(gn_hoa_100_sums, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "Bar Plot with Viridis Color Palette") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

gn_hoa_100_edge_comps <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"), 
                              values = c(1909.048655, 54.124039, 710.747350, 63.972691, 125.600177, 758.188026, 968.537631, 198.687025, 106.064457))
gn_hoa_100_edge_comps$norm_values <- gn_hoa_100_edge_comps$values/sum(gn_hoa_100_edge_comps$values)
ggplot(gn_hoa_100_edge_comps, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "gn_hoa_100_edge_comps") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave(filename = "gn_hoa_100_edge_comps.png",plot = last_plot(),dpi=1200)

hoa_100_sums <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"), 
                              values = c(1561, 42, 556, 61, 103, 659, 593, 296, 61))
hoa_100_sums$norm_values <- hoa_100_sums$values/sum(hoa_100_sums$values)
ggplot(hoa_100_sums, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "Bar Plot with Viridis Color Palette") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

hoa_100_edge_comps <- data.frame(hoa = c("Altered Intercellular Communication","Cellular Senescence","Deregulated Nutrient Sensing","Epigenetic Alterations","Genomic Instability","Loss of Proteostasis","Mitochondrial Dysfunction","Stem Cell Exhaustion","Telomere Attrition"), 
                           values = c(261008, 7189, 92828, 9993, 17593, 109704, 98599, 50332, 10501))
hoa_100_edge_comps$norm_values <- hoa_100_edge_comps$values/sum(hoa_100_edge_comps$values)
ggplot(hoa_100_edge_comps, aes(x = hoa, y = norm_values, fill = hoa)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "", y = "Fraction of Critical Edges",fill = "Hallmarks of Aging", title = "hoa_100_edge_comps") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave(filename = "hoa_100_edge_comps.png",plot = last_plot(),dpi=1200)

# Analyzing JSD Data
mean_jsd_csv_directory <- "/Users/administrator/Downloads/RawTables/JSD_Mean_Output"
sd_jsd_csv_directory <- "/Users/administrator/Downloads/RawTables/JSD_SD_Output"

mean_jsd_csv_files <- list.files(mean_jsd_csv_directory, pattern = "\\.csv$", full.names = TRUE)
sd_jsd_csv_files <- list.files(sd_jsd_csv_directory, pattern = "\\.csv$", full.names = TRUE)

mean_jsd_combined_data <- lapply(mean_jsd_csv_files, read.csv)
mean_second_columns <- lapply(mean_jsd_combined_data, function(df) df[, 2])
mean_jsd_combined_matrix <- do.call(cbind, mean_second_columns)
mean_jsd_column_order <- c()
for (name in 1:length(sig_gene_sets)) {
  col_gene_set <- names(mean_jsd_combined_data[[name]][2])
  mean_jsd_column_order <- c(mean_jsd_column_order,col_gene_set)
}
colnames(mean_jsd_combined_matrix) <- mean_jsd_column_order
rownames(mean_jsd_combined_matrix) <- sig_gene_sets
mean_jsd_combined_matrix <- mean_jsd_combined_matrix[, rownames(mean_jsd_combined_matrix)]

sd_jsd_combined_data <- lapply(sd_jsd_csv_files, read.csv)
sd_second_columns <- lapply(sd_jsd_combined_data, function(df) df[, 2])
sd_jsd_combined_matrix <- do.call(cbind, sd_second_columns)
sd_jsd_column_order <- c()
for (name in 1:length(sig_gene_sets)) {
  col_gene_set <- names(sd_jsd_combined_data[[name]][2])
  sd_jsd_column_order <- c(sd_jsd_column_order,col_gene_set)
}
colnames(sd_jsd_combined_matrix) <- sd_jsd_column_order
rownames(sd_jsd_combined_matrix) <- sig_gene_sets
sd_jsd_combined_matrix <- sd_jsd_combined_matrix[, rownames(sd_jsd_combined_matrix)]

mean_jsd_vector <- c(mean_jsd_combined_matrix)
sd_jsd_vector <- c(sd_jsd_combined_matrix)

diag_indices <- seq(1, length(mean_jsd_vector), by = nrow(mean_jsd_combined_matrix) + 1)
mean_jsd_vector_no_diag <- mean_jsd_vector[-diag_indices]
sd_jsd_vector_no_diag <- sd_jsd_vector[-diag_indices]

jsd_info <- data.frame(mean = mean_jsd_vector_no_diag, sd = sd_jsd_vector_no_diag, coeff_var = sd_jsd_vector_no_diag/mean_jsd_vector_no_diag)

mean_jsd_quantile_20 <- quantile(mean_jsd_vector_no_diag, probs = 0.2)
sd_jsd_quantile_20 <- quantile(sd_jsd_vector_no_diag, probs = 0.2)
coeff_var_quantile_20 <- quantile(jsd_info$coeff_var, na.rm = TRUE, probs = 0.2)
edge_indices <- which(jsd_info$mean < mean_jsd_quantile_20 & jsd_info$sd < sd_jsd_quantile_20)

ggplot(data = jsd_info, aes(x = mean, y = sd)) +
  geom_point(data = jsd_info[setdiff(1:nrow(jsd_info), edge_indices), ], color = "black", na.rm = TRUE) +
  geom_point(data = jsd_info[edge_indices, ], color = "blue", na.rm = TRUE)
  

write.csv(mean_jsd_combined_matrix, "/Users/administrator/Downloads/RawTables/mean_jsd_combined_matrix.csv")
write.csv(sd_jsd_combined_matrix, "/Users/administrator/Downloads/RawTables/sd_jsd_combined_matrix.csv")

coeff_jsd_csv_directory <- "/Users/administrator/Downloads/RawTables/JSD_Coeff_Output"
p_value_jsd_csv_directory <- "/Users/administrator/Downloads/RawTables/JSD_P_Value_Output"

coeff_jsd_csv_files <- list.files(coeff_jsd_csv_directory, pattern = "\\.csv$", full.names = TRUE)
p_value_jsd_csv_files <- list.files(p_value_jsd_csv_directory, pattern = "\\.csv$", full.names = TRUE)

coeff_jsd_combined_data <- lapply(coeff_jsd_csv_files, read.csv)
coeff_second_columns <- lapply(coeff_jsd_combined_data, function(df) df[, 2])
coeff_jsd_combined_matrix <- do.call(cbind, coeff_second_columns)
coeff_jsd_column_order <- c()
for (name in 1:length(sig_gene_sets)) {
  col_gene_set <- names(coeff_jsd_combined_data[[name]][2])
  coeff_jsd_column_order <- c(coeff_jsd_column_order,col_gene_set)
}
colnames(coeff_jsd_combined_matrix) <- coeff_jsd_column_order
rownames(coeff_jsd_combined_matrix) <- sig_gene_sets
coeff_jsd_combined_matrix <- coeff_jsd_combined_matrix[, rownames(coeff_jsd_combined_matrix)]

p_value_jsd_combined_data <- lapply(p_value_jsd_csv_files, read.csv)
p_value_second_columns <- lapply(p_value_jsd_combined_data, function(df) df[, 2])
p_value_jsd_combined_matrix <- do.call(cbind, p_value_second_columns)
p_value_jsd_column_order <- c()
for (name in 1:length(sig_gene_sets)) {
  col_gene_set <- names(p_value_jsd_combined_data[[name]][2])
  p_value_jsd_column_order <- c(p_value_jsd_column_order,col_gene_set)
}
colnames(p_value_jsd_combined_matrix) <- p_value_jsd_column_order
rownames(p_value_jsd_combined_matrix) <- sig_gene_sets
p_value_jsd_combined_matrix <- p_value_jsd_combined_matrix[, rownames(p_value_jsd_combined_matrix)]

coeff_jsd_vector <- c(coeff_jsd_combined_matrix)
p_value_jsd_vector <- c(p_value_jsd_combined_matrix)

diag_indices <- seq(1, length(coeff_jsd_vector), by = nrow(coeff_jsd_combined_matrix) + 1)
coeff_jsd_vector_no_diag <- coeff_jsd_vector[-diag_indices]
p_value_jsd_vector_no_diag <- p_value_jsd_vector[-diag_indices]

jsd_lm_info <- data.frame(coeff = coeff_jsd_vector_no_diag, p_value = p_value_jsd_vector_no_diag)
jsd_all_info <- cbind(jsd_info,jsd_lm_info)

jsd_bon_threshold <- 0.05/(length(sig_gene_sets)*(length(sig_gene_sets)-1)/2)
jsd_high_p <- jsd_all_info[jsd_all_info$p_val > jsd_bon_threshold,]

ggplot(data = jsd_high_p, aes(x = mean, y = p_value)) +
  geom_point()

write.csv(coeff_jsd_combined_matrix, "/Users/administrator/Downloads/RawTables/coeff_jsd_combined_matrix.csv")
write.csv(p_value_jsd_combined_matrix, "/Users/administrator/Downloads/RawTables/p_value_jsd_combined_matrix.csv")

# Compute variance across the different correlation matrices
var(c(pcc_matrix/2))
var(c(mean_jsd_combined_matrix))

## Other linear regression models (not incorporated into the final paper)
# perform_regression_gene_tissue_factor <- function(gene_expression_data, age_tissue_data, tissue_specific = FALSE) {
#   # Initialize a list to store regression results for each gene
#   regression_results <- list()
#   
#   # Loop over genes
#   for (gene_symbol in rownames(gene_expression_data)) {
#     # Filter gene expression data for the current gene
#     gene_expression <- gene_expression_data[gene_symbol,]
#     
#     # If gene expression data is available
#     if (!is.null(gene_expression)) {
#       # If tissue-specific regression is requested
#       if (tissue_specific) {
#         # Merge with age and tissue metadata for tissue-specific regression
#         merged_data <- merge(data.frame(Expression = gene_expression, SRR.ID = colnames(gene_expression_data)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
#         
#         # Perform linear regression within each tissue
#         regression_results[[gene_symbol]] <- by(merged_data, merged_data$Tissue, function(df) {
#           if (nrow(df) > 1) {
#             model <- lm(Expression ~ Age + Tissue, data = df)
#             coefficients <- coef(summary(model))
#             return(coefficients)
#           } else {
#             return(NULL)
#           }
#         })
#       } else {  # If tissue-specific regression is not requested
#         # Merge with age metadata for tissue-independent regression
#         merged_data <- merge(data.frame(Expression = gene_expression, SRR.ID = colnames(gene_expression_data)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
#         
#         # Perform tissue-independent linear regression
#         model <- lm(Expression ~ Age + Tissue, data = merged_data)
#         regression_results[[gene_symbol]] <- list(all = coef(summary(model)))
#       }
#     } else {  # If gene expression data is empty or NULL
#       message("Gene expression data is empty or NULL for gene: ", gene_symbol)
#       regression_results[[gene_symbol]] <- NULL
#     }
#   }
#   
#   return(regression_results)
# }
# 
# factor_gene_regression_results_tissue_independent <- perform_regression_gene_tissue_factor(raw_sort, age_sort_meta, tissue_specific = FALSE)
# factor_gene_regression_summary_tissue_independent <- do.call(rbind,factor_gene_regression_results_tissue_independent)
# factor_gene_regression_summary_tissue_independent <- as.data.frame(factor_gene_regression_summary_tissue_independent)
# for (i in 1:length(factor_gene_regression_summary_tissue_independent$all)) {
#   factor_gene_regression_summary_tissue_independent$intercept_estimate[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][1,1]
#   factor_gene_regression_summary_tissue_independent$intercept_std_error[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][1,2]
#   factor_gene_regression_summary_tissue_independent$intercept_t_value[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][1,3]
#   factor_gene_regression_summary_tissue_independent$intercept_significance[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][1,4]
#   factor_gene_regression_summary_tissue_independent$coefficient_estimate[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][2,1]
#   factor_gene_regression_summary_tissue_independent$coefficient_std_error[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][2,2]
#   factor_gene_regression_summary_tissue_independent$coefficient_t_value[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][2,3]
#   factor_gene_regression_summary_tissue_independent$coefficient_significance[i] <- factor_gene_regression_summary_tissue_independent$all[[i]][2,4]
# }
# 
# # Set up function to perform linear regression analysis on gene sets (default option is tissue-independent, but can also set to tissue-specific)
# perform_regression_factor <- function(gene_set, gene_expression_data, age_tissue_data, tissue_specific = FALSE) {
#   # Extract genes in the current gene set
#   gene_symbols <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
#   
#   # Filter gene expression data for genes in the current gene set
#   gene_set_expression <- gene_expression_data[rownames(gene_expression_data) %in% gene_symbols, ]
#   
#   result <- tryCatch({
#     if (!is.null(gene_set_expression) && nrow(gene_set_expression) > 0 && ncol(gene_set_expression) > 0) { # if the gene set can be pulled out from the data
#       # Aggregate gene expression data across genes within the gene set
#       gene_set_expression_aggregated <- apply(gene_set_expression, 2, mean, na.rm = TRUE)
#       
#       if (tissue_specific) {
#         # Merge with age and tissue metadata for tissue-specific regression
#         merged_data <- merge(data.frame(Expression = gene_set_expression_aggregated, "SRR.ID" = names(gene_set_expression_aggregated)), age_tissue_data, by = "SRR.ID", all.x = TRUE)
#         
#         # Perform linear regression within each tissue
#         regression_results <- by(merged_data, merged_data$Tissue, function(df) {
#           if (nrow(df) > 1) {
#             model <- lm(Expression ~ Age + Tissue, data = df)
#             coefficients <- coef(summary(model))
#             return(coefficients)
#           } else {
#             return(NULL)
#           }
#         })
#       } else {
#         # Merge with age metadata for tissue-independent regression
#         merged_data <- merge(data.frame(Expression = gene_set_expression_aggregated, "SRR.ID" = names(gene_set_expression_aggregated)), age_tissue_data[, c("SRR.ID", "Age","Tissue")], by = "SRR.ID", all.x = TRUE)
#         
#         # Perform tissue-independent linear regression
#         model <- lm(Expression ~ Age + Tissue, data = merged_data)
#         regression_results <- list(all = coef(summary(model)))
#       }
#       
#       return(list(gene_set = gene_set, regression_results = regression_results))
#     } else {  # If gene_set_expression is empty or NULL, return NULL
#       message("Gene set expression data is empty or NULL for gene set: ", gene_set)
#       return(NULL)
#     }
#   }, error = function(e) {  # Catch any errors that occur
#     message("Error occurred: ", e$message)
#     return(NULL)
#   })
#   
#   return(result)
# }
# # Note: some errors do occur, likely due to either empty sets or gene id discrepancies,
# # but they account for less than 1% of gene sets, and the function is set up to handle them.
# 
# # Perform regression analysis for each gene set (tissue-independent)
# factor_regression_results_tissue_independent <- lapply(unique(all_gene_sets$gs_name), function(gene_set) {
#   perform_regression_factor(gene_set, raw_sort, age_sort_meta)
# })
# 
# factor_regression_summary_tissue_independent <- do.call(rbind,factor_regression_results_tissue_independent)
# factor_regression_summary_tissue_independent <- as.data.frame(factor_regression_summary_tissue_independent)
# for (i in 1:length(factor_regression_summary_tissue_independent$regression_results)) {
#   factor_regression_summary_tissue_independent$intercept_estimate[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,1]
#   factor_regression_summary_tissue_independent$intercept_std_error[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,2]
#   factor_regression_summary_tissue_independent$intercept_t_value[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,3]
#   factor_regression_summary_tissue_independent$intercept_significance[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,4]
#   factor_regression_summary_tissue_independent$coefficient_estimate[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,1]
#   factor_regression_summary_tissue_independent$coefficient_std_error[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,2]
#   factor_regression_summary_tissue_independent$coefficient_t_value[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,3]
#   factor_regression_summary_tissue_independent$coefficient_significance[i] <- factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,4]
# }
# true_sig_gene_sets <- unlist(factor_regression_summary_tissue_independent$gene_set[which(factor_regression_summary_tissue_independent$coefficient_significance < bon_threshold)])
# 
# non_healthy_samples_factor_regression_results_tissue_independent <- lapply(unique(all_gene_sets$gs_name), function(gene_set) {
#   perform_regression_factor(gene_set, non_healthy_raw_sort, age_sort_non_healthy)
# })
# 
# non_healthy_samples_factor_regression_summary_tissue_independent <- do.call(rbind,non_healthy_samples_factor_regression_results_tissue_independent)
# non_healthy_samples_factor_regression_summary_tissue_independent <- as.data.frame(non_healthy_samples_factor_regression_summary_tissue_independent)
# for (i in 1:length(non_healthy_samples_factor_regression_summary_tissue_independent$regression_results)) {
#   non_healthy_samples_factor_regression_summary_tissue_independent$intercept_estimate[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,1]
#   non_healthy_samples_factor_regression_summary_tissue_independent$intercept_std_error[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,2]
#   non_healthy_samples_factor_regression_summary_tissue_independent$intercept_t_value[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,3]
#   non_healthy_samples_factor_regression_summary_tissue_independent$intercept_significance[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][1,4]
#   non_healthy_samples_factor_regression_summary_tissue_independent$coefficient_estimate[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,1]
#   non_healthy_samples_factor_regression_summary_tissue_independent$coefficient_std_error[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,2]
#   non_healthy_samples_factor_regression_summary_tissue_independent$coefficient_t_value[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,3]
#   non_healthy_samples_factor_regression_summary_tissue_independent$coefficient_significance[i] <- non_healthy_samples_factor_regression_summary_tissue_independent$regression_results[[i]][[1]][2,4]
# }
# non_healthy_samples_sig_gene_sets <- unlist(non_healthy_samples_factor_regression_summary_tissue_independent$gene_set[which(non_healthy_samples_factor_regression_summary_tissue_independent$coefficient_significance < bon_threshold)])
# 
# # Set up function to perform linear regression analysis on gene sets (default option is tissue-independent, but can also set to tissue-specific)
# perform_regression_factor <- function(gene_set, gene_expression_data, disease_data, age_tissue_disease_metadata) {
#   # Extract genes in the current gene set
#   gene_symbols <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
#   
#   # Filter gene expression data for genes in the current gene set
#   gene_set_expression <- gene_expression_data[rownames(gene_expression_data) %in% gene_symbols, ]
#   
#   result <- tryCatch({
#     if (!is.null(gene_set_expression) && nrow(gene_set_expression) > 0 && ncol(gene_set_expression) > 0) { # if the gene set can be pulled out from the data
#       # Aggregate gene expression data across genes within the gene set
#       gene_set_expression_aggregated <- apply(gene_set_expression, 2, mean, na.rm = TRUE)
#       
#       # Merge with age metadata for tissue-independent regression
#       merged_data <- merge(data.frame(Expression = gene_set_expression_aggregated, "SRR.ID" = names(gene_set_expression_aggregated)), age_tissue_disease_metadata[, c("SRR.ID", "Age","Tissue","Condition")], by = "SRR.ID", all.x = TRUE)
#       # Perform tissue-independent linear regression
#       model <- lm(Expression ~ Age*Condition + Tissue, data = merged_data)
#       regression_results <- list(all = coef(summary(healthy_model)))
#       
#       
#       return(list(gene_set = gene_set, regression_results = regression_results))
#     } else {  # If gene_set_expression is empty or NULL, return NULL
#       message("Gene set expression data is empty or NULL for gene set: ", gene_set)
#       return(NULL)
#     }
#   }, error = function(e) {  # Catch any errors that occur
#     message("Error occurred: ", e$message)
#     return(NULL)
#   })
#   
#   return(result)
# }
# # Note: some errors do occur, likely due to either empty sets or gene id discrepancies,
# # but they account for less than 1% of gene sets, and the function is set up to handle them.
# 
# # Perform regression analysis for each gene set (tissue-independent)
# factor_regression_results_tissue_independent <- lapply(unique(all_gene_sets$gs_name), function(gene_set) {
#   perform_regression_factor(gene_set, combat_data, meta_filtered)
# })
# 
# 
# condition_healthy <- which(names(coefficients(model))=="ConditionHealthy")
# age_condition_healthy <- which(names(coefficients(model))=="Age:ConditionHealthy")
# ctr <- rep(0, length(names(coefficients(model))))
# ctr[condition_healthy] <- -1
# ctr[condition_healthy+1] <- 1
# v_names <- paste(names(coefficients(model))[condition_healthy],names(coefficients(model))[condition_healthy+1],sep="-")
# ctr <- rbind(v_names = ctr)      
# summary(glht(model,ctr))
# woob <- linearHypothesis(model, "TissuePancreatic islet = TissueBlood;PBMC",singular.ok=TRUE)
