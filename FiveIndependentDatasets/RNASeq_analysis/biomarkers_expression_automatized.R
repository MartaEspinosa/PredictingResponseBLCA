# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 6th Oct 2023
# This scripts automatically checks the expression of 3 biomarker genes for which we know their expression in Responders and NonResponders
# Ideally, TFGB1 should be overexpressed in NR while the other two genes should be overexpressed in R

######################## libraries to be loaded #################################
library(dplyr)
library(edgeR)
library(tidyr)
library(ggplot2)
library(ggthemes)
################################################################################

################### info that always needs to be loaded ########################
DIR="/users/genomics/marta/BLCA"
GENOMEDIR = "/genomics/users/marta/genomes"
################################################################################

############################# projects and files ###############################
projects_names = c("UC-GENOME","IMvigor210","SNY-2017","UNC-108","HdM-BLCA-1")
wilcoxon_test_biomarkers_file = file.path(DIR,"RNASeq_expression/biomarkers/biomarkers_wilcoxontest.txt")
if (file.exists(wilcoxon_test_biomarkers_file)) {
  file.remove(wilcoxon_test_biomarkers_file)
}
################################################################################

############################# selected biomarkers ##############################
### CD8A - ENSG00000153563
### TGFB1 - ENSG00000105329
### CXCL9 - ENSG00000138755
biomarkers = data.frame("gene_id" = c('ENSG00000153563', 'ENSG00000105329', 'ENSG00000138755'),
                        "gene_name" = c("CD8A","TGFB1","CXCL9"))
all_biomarkers_expr = data.frame("gene_name" = character(),
                                 "patient" = character(),
                                 "TPM" = numeric(),
                                 "logTPM" = numeric(),
                                 "Dataset" = factor())
biomarkers_pvalues = data.frame("gene_name" = character(),
                                "Dataset" = factor(),
                                "Wilcox-pvalues" = numeric())
################################################################################

for(i in 1:length(projects_names)) {
  print(projects_names[i])
  patient_file = read.csv(file.path(DIR,projects_names[i],"patients_response.csv"))

  tableofcounts = read.csv(file.path(DIR,projects_names[i],"analysis/01_counts/TPMs_genenames.csv"))
  
  # subset biomarkers expression
  biomarkers_expression = tableofcounts %>% subset(gene_name %in% biomarkers$gene_name)
  biomarkers_expression = biomarkers_expression %>% pivot_longer(cols=-c("gene_name"), names_to = "patient", values_to = "TPM")
  # correct patients IDs if needed
  biomarkers_expression$patient = gsub("X","", biomarkers_expression$patient)
  if(projects_names[i] == "HdM-BLCA-1") {
    biomarkers_expression$patient = gsub("_.*","", biomarkers_expression$patient)
  }
  # merge data with patients information to know whether a patient is R or NR
  biomarkers_expression = merge(biomarkers_expression, patient_file, by="patient")
  biomarkers_expression$logTPM <- log10(biomarkers_expression$TPM)
  biomarkers_expression$Dataset = projects_names[i]
  
  all_biomarkers_expr = rbind(all_biomarkers_expr, biomarkers_expression)
  ggplot(biomarkers_expression, aes(x=Response, y=logTPM)) +
    geom_boxplot(aes(fill=Response),width=0.5) +
    ggtitle(projects_names[i]) +
    scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size=16)) +
    facet_wrap(~gene_name, scales="free",nrow=3, ncol=1)
  ggsave(file.path(DIR,projects_names[i],"/results/biomarkers.pdf"), width=4.06, height=6.15)
  ggsave(file.path(DIR,projects_names[i],"/results/biomarkers.png"), width=4.06, height=6.15)
  
  # compute a Wilcoxon test to see if there is a statistical association between the expression of each biomarker with the response to immunotherapy
  CD8A = biomarkers_expression[biomarkers_expression$gene_name == "CD8A",]
  CD8A_test = wilcox.test(TPM ~ Response, data=CD8A)

  TGFB1 = biomarkers_expression[biomarkers_expression$gene_name == "TGFB1",]
  TGFB1_test = wilcox.test(TPM ~ Response, data=TGFB1)
  
  CXCL9 = biomarkers_expression[biomarkers_expression$gene_name == "CXCL9",]
  CXCL9_test = wilcox.test(TPM ~ Response, data=CXCL9)
  
  temp_tests = data.frame("gene_name" = c("CD8A","TGFB1","CXCL9"),
                          "Dataset" = c(projects_names[i],projects_names[i],projects_names[i]),
                          "Wilcox-pvalues" = c(CD8A_test$p.value,TGFB1_test$p.value,CXCL9_test$p.value))
  biomarkers_pvalues = rbind(biomarkers_pvalues, temp_tests)
}
# save wilcoxon p values in a text file
write.csv(biomarkers_pvalues, file=wilcoxon_test_biomarkers_file, row.names = F, quote = F) 


ggplot(all_biomarkers_expr, aes(x=Dataset, y=logTPM)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("Biomarkers expression") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16)) +
  facet_wrap(~gene_name, scales="free",nrow=3, ncol=1)
ggsave(file.path(DIR,"RNASeq_expression/biomarkers/biomarkers.png"), width=6.99, height=5.85)
ggsave(file.path(DIR,"RNASeq_expression/biomarkers/biomarkers.pdf"), width=6.99, height=5.85)

temp = all_biomarkers_expr %>% select(patient, Response, Dataset) %>% unique()
table(temp$Response, temp$Dataset)
