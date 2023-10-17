library(tidyverse)
library(dplyr)
library(ggplot2)
library(rcartocolor)
require(QualityGraphs)
require(RColorBrewer)
require(stringr)
library("pheatmap")

projects=c("BLCA")
cohorts=c("BLCA")

wd = "/users/genomics/marta"
# wd = "~/Downloads/IMIM"

# biomart = read.csv("/genomics/users/marta/genomes/transcript_gene.txt")
biomart = read.csv("/genomics/users/marta/genomes/transcript_gene_v41.txt")
names(biomart) = c("gene_id","gene_name", "gene_type")

# get all lncRNA and protein coding transcripts tumor-specific
tumorspecific = data.frame("gene_name" = as.character(),
                           stringsAsFactors = F)

for (p in 1:length(projects)) {
  proteincoding = read.csv(file.path(wd,projects[p],"analysis/08_tumor_specific/straightforward/common_tumorspecific_coding_genes_GTExvalidated_testis+ovary05.csv"))
  proteincoding$project = cohorts[p]
  proteincoding$biotype = "Protein Coding"
  tumorspecific = rbind(tumorspecific, proteincoding)
  
  lncRNA = read.csv(file.path(wd,projects[p],"analysis/08_tumor_specific/straightforward/common_tumorspecific_noncoding_genes_GTExvalidated_testis+ovary05.csv"))
  lncRNA$project = cohorts[p] 
  lncRNA$biotype = "lncRNA"
  tumorspecific = rbind(tumorspecific, lncRNA)
}

table(tumorspecific$biotype)

### patients information

patients = read.csv(file.path(wd,projects[p],"results/paired_patients.csv"))

patients$project = "BLCA"
full_patients = patients


### HOW MANY TUMOR-SPECIFIC TRANSCRIPTS HAVE EXPRESSION IN PATIENTS WHERE IT IS NOT TUMOR-SPECIFIC
tc = read.csv(file.path(wd,"BLCA/analysis/07_quantification/straightforward/TPMs_genenames.csv"))

# pivot longer whole table of counts (both normal and tumor)
tc_long = tc %>% pivot_longer(cols=!c(gene_name), names_to = "sample", values_to = "TPM") 
# tc_long = tc_long %>% select(!gene_name)

# patients info
long_patients = full_patients %>% pivot_longer(cols=!c(patient, project), names_to = "normal_tumor", values_to = "sample")
long_patients$sample = gsub("-",".",long_patients$sample)
long_patients$patient = gsub("-",".",long_patients$patient)

# merge whole table of counts with patient information
tc_toanalyse_patients = merge(tc_long, long_patients, on="sample")

# select transcripts that are tumor-specific in some patient
tumorsp_in_some = merge(tc_toanalyse_patients, tumorspecific, on="gene_name")

tumorsp_in_some = tumorsp_in_some %>% select(-sample)
tumorsp_in_some = unique(tumorsp_in_some)
tumorsp_in_some = na.omit(tumorsp_in_some)
tumorsp_in_some_wide = tumorsp_in_some %>% pivot_wider(names_from = "normal_tumor", values_from = "TPM")

tumorsp_in_some_wide %>% head

# # which are tumor-specific? classify again
tumorsp_in_some_wide = tumorsp_in_some_wide %>% mutate(tumor_specific_1TPM = case_when((normal <= "0.1") & (tumor >= "1") ~ "YES",
                                                                                  (normal < "1") & (tumor < "1") ~ "case1",
                                                                                  (normal > "0.1") & (normal < "1") & (tumor >= "1") ~ "case2",
                                                                                  (normal > "1")  ~ "case3"))
tumorsp_in_some_wide = tumorsp_in_some_wide %>% mutate(tumor_specific_3TPM = case_when((normal <= "0.1") & (tumor >= "3") ~ "YES",
                                                                                        (normal < "1") & (tumor < "3") ~ "case1",
                                                                                        (normal > "0.1") & (normal < "1") & (tumor >= "3") ~ "case2",
                                                                                        (normal > "1")  ~ "case3"))
table(tumorsp_in_some_wide$tumor_specific_1TPM)
table(tumorsp_in_some_wide$tumor_specific_3TPM)

tumorsp_in_some_wide_genename = merge(tumorsp_in_some_wide, biomart, on="gene_name")
tumorsp_in_some_wide_genename %>% head

# tumorsp_in_some_wide_genename = tumorsp_in_some_wide_genename %>% select(patient, gene_name, forced_tumor_expression)
tumorsp_in_some_wide_genename = tumorsp_in_some_wide_genename %>% select(patient, gene_name, tumor_specific_3TPM, gene_type)

new_table_cases = tumorsp_in_some_wide_genename %>% pivot_wider(names_from = patient, values_from = tumor_specific_3TPM)
new_table_cases = as.data.frame(new_table_cases)


new_table_cases$case1 <- apply(new_table_cases, 1, function(x) length(which(x=="case1")))
new_table_cases$case2 <- apply(new_table_cases, 1, function(x) length(which(x=="case2")))
new_table_cases$case3 <- apply(new_table_cases, 1, function(x) length(which(x=="case3")))
new_table_cases$YES <- apply(new_table_cases, 1, function(x) length(which(x=="YES")))

new_table_cases = new_table_cases %>% select(gene_name, gene_type, case1, case2, case3, YES)
nrow(new_table_cases)

# > 10 % patients
new_table_cases = new_table_cases %>% subset(YES > 1)
nrow(new_table_cases)
write.csv(new_table_cases[order(new_table_cases$YES, decreasing = TRUE), ], "/users/genomics/marta/BLCA/analysis/08_tumor_specific/straightforward/tumorspecific_cases_proportions_3TPM_10percent.csv")


new_table_cases$n = rowSums(new_table_cases %>% select(3,4,5,6))
new_table_cases = new_table_cases %>% subset(n == 18)
## calculate proportions
new_table_cases$proportion_case1 <- new_table_cases$case1/18
new_table_cases$proportion_case2 <- new_table_cases$case2/18
new_table_cases$proportion_case3 <- new_table_cases$case3/18
new_table_cases$proportion_YES <- new_table_cases$YES/18
new_table_cases %>% head

# get those whose proportion of case 3 is > 5%
abundant_case3 = new_table_cases %>% subset(proportion_case3 > 0.05) 
write.csv(abundant_case3[order(abundant_case3$YES, decreasing = TRUE), ], "/users/genomics/marta/BLCA/analysis/08_tumor_specific/straightforward/removed_abundant_case3_3TPM.csv", row.names = F, quote = F)

no_abundant_case3 = new_table_cases %>% subset(!proportion_case3 > 0.05)
write.csv(no_abundant_case3[order(no_abundant_case3$YES, decreasing = TRUE), ], "/users/genomics/marta/BLCA/analysis/08_tumor_specific/straightforward/without_abundant_case3_3TPM.csv", row.names = F, quote = F)
table(no_abundant_case3$gene_type)

