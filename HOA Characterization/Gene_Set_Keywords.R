## List of descriptions of the significant gene sets ---------------------------

# Run pipeline to acquire data tables
source('Expression_by_Age_Pipeline (v2).R')

# Collect the names of the significant gene sets
sig_gene_set_names = rownames(sig_gene_set_expression)

index = c()
for (name in 1:length(sig_gene_set_names)){
  index = append(index, which(all_gene_sets$gs_name == sig_gene_set_names[name]))
}

# get descriptions for significant gene sets from the all_gene_sets table
sig_gene_set_description = all_gene_sets[index,c('gs_name','gs_description')]
sig_gene_set_description = unique(sig_gene_set_description)

# modify descriptions to be all lowercase, remove special characters
sig_gene_set_description$gs_description = tolower(sig_gene_set_description$gs_description)
sig_gene_set_description$gs_description = str_replace_all(sig_gene_set_description$gs_description, "[[:punct:]]", "")

# function to find pathways that include certain keywords in the brief descriptions
Pways = function (keyword){
  noquote(print(paste('Keyword =',keyword)))
  sig_gene_set_description[grep(keyword,sig_gene_set_description$gs_description), ]
}

# Ex. Pways('cancer') 



## Search full descriptions of each pathway for keywords of interest ----------------

library(rvest) 

# fixing some names to be the same as the msigdb website 
sig_gene_set_names[510] = 'WP_BONE_MORPHOGENIC_PROTEIN_SIGNALING_AND_REGULATION'
sig_gene_set_names[513] = 'WP_CELL_TYPE_DEPENDENT_SELECTIVITY_OF_CCK2R_SIGNALING'
sig_gene_set_names[522] = 'WP_FOXP3_IN_COVID_19'
sig_gene_set_names[531] = 'WP_LINOLEIC_ACID_METABOLISM_AFFECTED_BY_SARS_COV_2'
sig_gene_set_names[532] = 'WP_LNCRNA_MEDIATED_MECHANISMS_OF_THERAPEUTIC_RESISTANCE'
sig_gene_set_names[540] = 'WP_NO_CGMP_PKG_MEDIATED_NEUROPROTECTION'
sig_gene_set_names[542] = 'WP_PATHOGENESIS_OF_SARS_COV_2_MEDIATED_BY_NSP9_NSP10_COMPLEX'
sig_gene_set_names[544] = 'WP_PKC_GAMMA_CALCIUM_SIGNALING_PATHWAY_IN_ATAXIA'
sig_gene_set_names[545] = 'WP_PRADER_WILLI_AND_ANGELMAN_SYNDROME'
sig_gene_set_names[548] = 'WP_SARS_COV_2_ALTERING_ANGIOGENESIS_VIA_NRP1'
sig_gene_set_names[549] = 'WP_SELECTIVE_EXPRESSION_OF_CHEMOKINE_RECEPTORS_DURING_T_CELL_POLARIZATION'
sig_gene_set_names[556] = 'WP_VITAMIN_D_SENSITIVE_CALCIUM_SIGNALING_IN_DEPRESSION'

# Retrieve full descriptions of each significant gene set from msigdb website
long_descriptions = c()
for (set in sig_gene_set_names){
  url = paste('https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/',set,'.html', sep='')
  # Read the HTML content of the website 
  webpage <- read_html(url) 
  # Select the table using CSS selector 
  table_node <- html_nodes(webpage, "table") 
  # Extract the table content 
  table_content <- html_table(table_node)[[2]] 
  extracted_text = table_content[[4,2]]
  long_descriptions = append(long_descriptions, extracted_text)
}

# add long descriptions to table
sig_gene_set_description$long_description = long_descriptions

# if long description is blank, repeat brief description
for (i in (1:length(sig_gene_set_names))){
  if (sig_gene_set_description$long_description[i] == ''){
    sig_gene_set_description$long_description[i] = sig_gene_set_description$gs_description[i]
  }
}

# modify long descriptions to be all lowercase, remove special characters
sig_gene_set_description$long_description = tolower(sig_gene_set_description$long_description)
sig_gene_set_description$long_description = str_replace_all(sig_gene_set_description$long_description, "[[:punct:]]", "")


# function to find pathways that include certain keyword in the long descriptions
# can print a list of the pathways+descriptions, and export to .csv
LongPways = function (keyword){
  noquote(print(paste('Keyword =',keyword)))
  pways = sig_gene_set_description[grep(keyword,sig_gene_set_description$long_description), ]
  name = paste('Pathways_relating_to_',keyword, sep='')
  return(pways)
  #write.csv(pways, paste(name,'.csv', sep=''))
}

# ex. LongPways('cell cycle')

# keywords related to known age-associated pathways
Keywords = c(
  'cell cycle',
  'mitosis|mitotic',
  'mitochondria',
  'lifespan|healthspan',
  'nucleus|nuclear|lamin',
  'telomere|shelterin',
  'methylation|acetylation', # big hitter
  'epigenetic|chromatin', 
  'insulin|igf|akt|mtor|pi3k', 
  'sirt',
  'nfkb|wasp|wave|cgas|sting',
  'transcription|gene expression',
  'mirna',
  'progeroid|progeria',
  'heat shock|folding|folded|chaperone',
  'degeneration',
  'phagy', # as in autophagy
  'growth factor|growth hormone',
  'ros|pgc',
  'senescence|ink4|arf|rb|p16|p53',
  'stem cell',
  'inflam', # as in inflammation
  'immun|immune',
  'neur',
  'signaling|cytokine|chemokine',
  'microbe|microbiome|bacteria|virus',
  'dna damage|lesion',
  'transposon',
  'ecm|junctions|adhesion',
  'translation',
  'lysosome|degradation|proteolysis',
  'foxo|e26',
  'cardiac|diabetes|alzheimers', # plus any other non-cancer chronic diseases
  'cancer|oncogen',
  'adipose|adipocytes|fat|triglyceride|lipid',
  'reprogramming|differentiation|pluripotency',
  'wound', # repair/healing
  'ras|raf|mek|erk|notch|tgf|smad|hedgehog|gli',
  'ifn|tnf',
  'dna repair',
  'age|aging|ageing', 
  'proteostasis')


# calculate the number of outputs in each category
# run loop for each input, output number of pathways in each category
hits = c()
for (entry in Keywords){
  hits = append(hits, length(LongPways(entry)$gs_name))
}

keyword_hits = data.frame(Keywords = Keywords, Hits = hits)

# Make keywords more legible for reading on graph axis
keyword_hits$Keywords = gsub('\\|', '/', keyword_hits$Keywords)
keyword_hits$Keywords = toupper(keyword_hits$Keywords)

# graph of categories and # of outputs
png(filename = 'GeneSetKeywords.png', 
    units = 'px',
    width = 800*(1200/72),
    height = 600*(1200/72),
    res = 1200)
par(mar = c(4,12,1,1))
barplot(keyword_hits$Hits, 
        xlab = 'Number of Gene Sets',
        xlim = c(0, 200),
        main = "Number of Significant Gene Sets Per Annotation Category", 
        horiz = TRUE,
        names.arg = keyword_hits$Keywords, 
        las = 1,
        cex.names = 0.5,
        col = '#7AD151FF')
dev.off()


# Function for pulling out pathways and descriptions of a given keyword annotation 
DescribePways = function(keyword){
  LongPways(keyword)[,c('gs_name','long_description')]
}






