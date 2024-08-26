task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
library(philentropy)
raw_sort_positive <- read.csv("raw_sort_positive.csv",sep=",")
regression_summary_independent <- read.csv("regression_summary_independent.csv",sep=",")
all_gene_sets <- read.csv("all_gene_sets.csv",sep=",")
age_sort_meta <- read.csv("age_sort_meta.tsv",sep="\t")
sig_gene_sets <- unlist(regression_summary_independent$gene_set[abs(regression_summary_independent$coefficient_estimate) > 0.02])
jsd_coeff <- matrix(NA, nrow = length(sig_gene_sets), ncol = 1,
                   dimnames = list(sig_gene_sets, sig_gene_sets[task_id]))
jsd_p_value <- matrix(NA, nrow = length(sig_gene_sets), ncol = 1,
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
  age_sort_meta$jsd <- jsd_values_across_age
  jsd_lm <- lm(jsd~Age,data=age_sort_meta)
  regression_results <- list(all = coef(summary(jsd_lm)))
  jsd_coeff[gene_set1, gene_set2] <<- regression_results[[1]][2,1]
  jsd_p_value[gene_set1, gene_set2] <<- regression_results[[1]][2,4]
})
coeff_variable_name <- paste("coeff_jsd_", task_id, ".csv", sep = "")
p_value_variable_name <- paste("p_value_jsd_", task_id, ".csv", sep = "")
write.csv(jsd_coeff, coeff_variable_name, row.names = TRUE)
write.csv(jsd_p_value, p_value_variable_name, row.names = TRUE)