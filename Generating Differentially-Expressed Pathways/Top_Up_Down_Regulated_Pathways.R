## Examine top 10 up- and down-regulated pathways in the tissue-independent dataset -------

# Run pipeline to acquire data tables
source('Expression_by_Age_Pipeline.R')

# Sort by coefficient_estimate, look at top up/down gene sets
top_up = regression_summary_independent[order(regression_summary_independent$coefficient_estimate, decreasing = TRUE), ]

top_up = top_up[, c('gene_set','coefficient_estimate')]
str(top_up[1:10, 1])
print(top_up[1:10, ])


top_down = regression_summary_independent[order(regression_summary_independent$coefficient_estimate, decreasing = FALSE), ]

top_down = top_down[,c('gene_set','coefficient_estimate')]
str(top_down[1:10, 1])
print(top_down[1:10, ])
