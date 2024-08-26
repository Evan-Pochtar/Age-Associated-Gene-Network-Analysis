task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(philentropy)
raw_sort_positive <- read.csv("raw_sort_positive.csv",sep=",")
regression_summary_independent <- read.csv("regression_summary_independent.csv",sep=",")
all_gene_sets <- read.csv("all_gene_sets.csv",sep=",")
sig_gene_sets <- unlist(regression_summary_independent$gene_set[abs(regression_summary_independent$coefficient_estimate) > 0.02])
sig_gene_set_expression <- do.call(rbind, lapply(sig_gene_sets, function(gene_set) {
  gene_symbols <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set]
  gene_set_expression <- raw_sort_positive[rownames(raw_sort_positive) %in% gene_symbols, ]
  gene_set_means <- data.frame(t(colMeans(gene_set_expression)))
  row.names(gene_set_means) <- gene_set
  return(gene_set_means)
}))
mean_jsd <- matrix(NA, nrow = length(sig_gene_sets), ncol = 1,
                   dimnames = list(sig_gene_sets, sig_gene_sets[task_id]))
sd_jsd <- matrix(NA, nrow = length(sig_gene_sets), ncol = 1,
                 dimnames = list(sig_gene_sets, sig_gene_sets[task_id]))
jsd_values <- lapply(sig_gene_sets, function(gene_set1) {
  gene_symbols1 <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set1]
  gene_set_expression1 <- raw_sort_positive[rownames(raw_sort_positive) %in% gene_symbols1, ]
  gene_set2 <- sig_gene_sets[task_id]
  gene_symbols2 <- all_gene_sets$gene_symbol[all_gene_sets$gs_name == gene_set2]
  gene_set_expression2 <- raw_sort_positive[rownames(raw_sort_positive) %in% gene_symbols2, ]
  jsd_values_across_age <- sapply(colnames(raw_sort_positive), function(age_column_name) {
    expr1 <- gene_set_expression1[, age_column_name]
    expr2 <- gene_set_expression2[, age_column_name]
    if (all(expr1 == 0) || all(expr2 == 0)) {
      return(NA)
    } else {
      expr_table <- t(cbind(expr1, expr2))
      suppressMessages({gJSD(expr_table, est.prob = "empirical")})
    }
  })
  mean_jsd[gene_set1, gene_set2] <<- mean(jsd_values_across_age, na.rm = TRUE)
  sd_jsd[gene_set1, gene_set2] <<- sd(jsd_values_across_age, na.rm = TRUE)
})
mean_variable_name <- paste("mean_jsd_", task_id, ".csv", sep = "")
sd_variable_name <- paste("sd_jsd_", task_id, ".csv", sep = "")
write.csv(mean_jsd, mean_variable_name, row.names = TRUE)
write.csv(sd_jsd, sd_variable_name, row.names = TRUE)