## By-Disease Analysis --------

library(tidyverse)

### Run pipeline to acquire data tables
#source('Expression_by_Age_Pipeline.R')

meta_filtered <- read_excel("filtered_metadata.xlsx") # metadata


# Number of each disease type
n_diseases = meta_filtered %>% group_by(Condition) %>% summarise(Count = n())
n_diseases = n_diseases %>% arrange(desc(Count))

print(n_diseases, n=28)

disease_names = unique(n_diseases$Condition)


## Summary Tables --------------------------------------------------------------
# number of old/adult/young samples for a given condition
table(meta_filtered$label[meta_filtered$Condition =="Alzheimer's"])


# list of ages of a given condition
meta_filtered$Age[meta_filtered$Condition=="Alzheimer's"]


summary = meta_filtered %>% dplyr::select(Age, Condition, label)
summary %>% group_by(Condition)

# age distribution of each tissue
summary(meta_filtered$Age[meta_filtered$Condition=="Alzheimer's"])

# age values included within each age label (young/adult/old)
range(summary$Age[summary$label=='young'])
range(summary$Age[summary$label=='adult'])
range(summary$Age[summary$label=='old'])

summary$Condition <- factor(summary$Condition , levels=disease_names)

# Summary table with first 9 (based on # samples) conditions

# Summary table with all conditions, sorted by decreasing # of samples
summary %>% filter(Condition == disease_names[1] |
                     Condition == disease_names[2] |
                     Condition == disease_names[3] |
                     Condition == disease_names[4] |
                     Condition == disease_names[5] |
                     Condition == disease_names[6] |
                     Condition == disease_names[7] |
                     Condition == disease_names[8] |
                     Condition == disease_names[9] |
                     Condition == disease_names[10] |
                     Condition == disease_names[11] |
                     Condition == disease_names[12] |
                     Condition == disease_names[13] |
                     Condition == disease_names[14] |
                     Condition == disease_names[15] |
                     Condition == disease_names[16] |
                     Condition == disease_names[17] |
                     Condition == disease_names[18] |
                     Condition == disease_names[19] |
                     Condition == disease_names[20] |
                     Condition == disease_names[21] |
                     Condition == disease_names[22] |
                     Condition == disease_names[23] |
                     Condition == disease_names[24] |
                     Condition == disease_names[25] |
                     Condition == disease_names[26] |
                     Condition == disease_names[27] |
                     Condition == disease_names[28]) %>% tbl_summary(by = Condition)


## Box Plots -------------------------------------------------------------------
png(filename = 'Age~Condition Box Plots.png', 
    units = 'px',
    width = 1370*(1200/72),
    height = 600*(1200/72),
    res = 1200)
par(mar = c(12,4,1,1))
b = boxplot(Age~Condition, data=summary, las=2, cex=0.8, cex.axis=0.6, 
            xlab = '', main = 'Age Distributions of Each Condition', col = '#414487FF')

y = b$stats[nrow(b$stats),] + 2.5
y[1] = 103 # adjust where the sample # falls for Healthy so not obstructed by outliers

text( 
  x=c(1:length(disease_names)), 
  y=y, 
  paste("n = ",table(summary$Condition),sep=""),
  cex = 0.6
)
dev.off()


# horizontal
png(filename = 'Age~Condition Box Plots, Horizontal.png', 
    units = 'px',
    width = 1370*(1200/72),
    height = 600*(1200/72),
    res = 1200)
par(mar = c(4,12,1,1))
b = boxplot(Age~Condition, data=summary, las=2, cex=0.5, cex.axis=0.6, 
            ylab = '', main = 'Age Distributions of Each Condition', 
            horizontal = TRUE, col = '#414487FF')
text( 
  y=c(1:length(disease_names)), 
  x=-2, 
  paste("n = ",table(summary$Condition),sep=""),
  cex = 0.6
)
dev.off()

