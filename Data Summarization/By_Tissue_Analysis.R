## By-Tissue Analysis --------

library(tidyverse)

### Run pipeline to acquire data tables
#source('Expression_by_Age_Pipeline.R')

meta_filtered <- read_excel("filtered_metadata.xlsx") # metadata

# Number of each tissue type
n_tissues = meta_filtered %>% group_by(Tissue) %>% summarise(Count = n())
n_tissues = n_tissues %>% arrange(desc(Count))

print(n_tissues, n=32)

tissue_names = unique(n_tissues$Tissue)

## Summary Tables --------------------------------------------------------------
# number of old/adult/young samples for a given tissue
table(meta_filtered$label[meta_filtered$Tissue =='Liver;liver hepatocytes'])


# list of ages of a given tissue
meta_filtered$Age[meta_filtered$Tissue=='Liver;liver hepatocytes']


summary = meta_filtered %>% dplyr::select(Age, Tissue, label)
summary %>% group_by(Tissue)

# age distribution of each tissue
summary(meta_filtered$Age[meta_filtered$Tissue=='Liver;liver hepatocytes'])

# age values included within each age label (young/adult/old)
range(summary$Age[summary$label=='young'])
range(summary$Age[summary$label=='adult'])
range(summary$Age[summary$label=='old'])

summary$Tissue <- factor(summary$Tissue , levels=tissue_names)


# Summary table with all tissues, sorted by decreasing # of samples
summary %>% filter(Tissue == tissue_names[1] |
                     Tissue == tissue_names[2] |
                     Tissue == tissue_names[3] |
                     Tissue == tissue_names[4] |
                     Tissue == tissue_names[5] |
                     Tissue == tissue_names[6] |
                     Tissue == tissue_names[7] |
                     Tissue == tissue_names[8] |
                     Tissue == tissue_names[9] |
                     Tissue == tissue_names[10] |
                     Tissue == tissue_names[11] |
                     Tissue == tissue_names[12] |
                     Tissue == tissue_names[13] |
                     Tissue == tissue_names[14] |
                     Tissue == tissue_names[15] |
                     Tissue == tissue_names[16] |
                     Tissue == tissue_names[17] |
                     Tissue == tissue_names[18] |
                     Tissue == tissue_names[19] |
                     Tissue == tissue_names[20] |
                     Tissue == tissue_names[21] |
                     Tissue == tissue_names[22] |
                     Tissue == tissue_names[23] |
                     Tissue == tissue_names[24] |
                     Tissue == tissue_names[25] |
                     Tissue == tissue_names[26] |
                     Tissue == tissue_names[27] |
                     Tissue == tissue_names[28] |
                     Tissue == tissue_names[29] |
                     Tissue == tissue_names[30] |
                     Tissue == tissue_names[31] |
                     Tissue == tissue_names[32]) %>% tbl_summary(by = Tissue)


## Box Plots -------------------------------------------------------------------
png(filename = 'Age~Tissue Box Plots.png', 
    units = 'px',
    width = 1370*(1200/72),
    height = 600*(1200/72),
    res = 1200)
par(mar = c(12,4,1,1))
b = boxplot(Age~Tissue, data=summary, las=2, cex=0.8, cex.axis=0.5, 
            xlab='', main = 'Age Distributions of Each Tissue', col = '#22A884FF')
text( 
  x=c(1:length(tissue_names)), 
  y=b$stats[nrow(b$stats),] + 2, 
  paste("n = ",table(summary$Tissue),sep=""),
  cex = 0.6
)
dev.off()


# horizontal
png(filename = 'Age~Tissue Box Plots, Horizontal.png', 
    units = 'px',
    width = 1370*(1200/72),
    height = 600*(1200/72),
    res = 1200)
par(mar = c(4,12,1,1))
b = boxplot(Age~Tissue, data=summary, las=2, cex=0.5, cex.axis=0.5, 
            ylab='', main = 'Age Distributions of Each Tissue', 
            horizontal = TRUE, col = '#22A884FF')
text( 
  y=c(1:length(tissue_names)), 
  x=-2, 
  paste("n = ",table(summary$Tissue),sep=""),
  cex = 0.6
)
dev.off()

