library(ggplot2)
library(tidyr)
library(dplyr)

DIR= "/users/genomics/marta/BLCA"
biomart = read.csv("/genomics/users/marta/genomes/transcript_gene_v41.txt")
names(biomart) = c("gene_id","gene_name", "gene_type")

## tumor-expressed
data = read.csv(file.path(DIR,"analysis/07_quantification/straightforward/tumor_1TPM.csv"))
data %>% head
### tumorspecific
tumorsp = read.csv(file.path(DIR,"analysis/08_tumor_specific/straightforward/without_abundant_case3_3TPM.csv"))

nrow(tumorsp)
table(tumorsp$gene_type)
tumorsp_genename = merge(tumorsp, data, by=c("gene_name","gene_type"))
tumorsp_genename = tumorsp_genename %>% unique()
tumorsp_genename = tumorsp_genename %>% select(-c(gene_id))
tumorsp_genename = tumorsp_genename %>% unique()
nrow(tumorsp_genename)
tumorsp_genename %>% head

############################# TABLE OF COUNTS ###############################
############################# projects and files ###############################
projects_names = c("UC-GENOME","IMvigor210","SNY-2017","UNC-108","HdM-BLCA-1")
################################################################################

for(i in 1:length(projects_names)) {
  print(projects_names[i])
  
  # read file
  tableofcounts = read.csv(file.path(DIR,projects_names[i],"analysis/01_counts/TPMs_genenames.csv"))
  names(tableofcounts) = sub("X","",names(tableofcounts))
  
  # select from table of counts the tumor-specific ones
  tableofcounts_tumorspecific = tableofcounts %>% subset(gene_name %in% tumorsp_genename$gene_name)
  write.csv(tableofcounts_tumorspecific, file=file.path(DIR,projects_names[i],"tableofcounts_3TPM_tumorspecific_nocase3.csv"), row.names = FALSE, quote = FALSE) #specify the path and filename for the output to be saved
  
  # count how many patients have the gene expressed > 1 TPM
  tableofcounts_tumorspecific$n = apply(tableofcounts_tumorspecific[, -1] > 3, 1, sum)
  tableofcounts_tumorspecific = tableofcounts_tumorspecific %>% select(gene_name, n)
  tableofcounts_tumorspecific %>% head
  
  # rename with the name of the dataset
  names(tableofcounts_tumorspecific)[2] = projects_names[i]
  
  # add new columns to count how many tumor-specific genes are tumor-specific in how many cohorts
  tumorsp_genename = merge(tumorsp_genename, tableofcounts_tumorspecific, by="gene_name", all.x = T)
  tumorsp_genename[is.na(tumorsp_genename)] <- 0
  
}  

tsp_in_cohorts = tumorsp_genename %>% select(gene_name, gene_type, projects_names)
tsp_in_cohorts$total_n = rowSums(tsp_in_cohorts[,3:7])
tsp_in_cohorts$total_percentage = round((tsp_in_cohorts$total_n / (16+20+235+84+64) )*100,2)
write.csv(tsp_in_cohorts, file=file.path(DIR,"RNASeq_expression/tumorspecific_3TPM_cohorts.csv"), row.names = F, quote = F)

tsp_in_cohorts %>% head
table(tsp_in_cohorts$gene_type)
