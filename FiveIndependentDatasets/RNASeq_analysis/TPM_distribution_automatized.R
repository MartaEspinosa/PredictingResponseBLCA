# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 28th Sept 2023
# This scripts automatically checks TPM distribution of lncRNA + pseudogenes versus protein-coding genes
# Ideally, the non-coding genes selected should be less expressed than the protein-coding ones
# They usually cross arround logTPM = 0

######################## libraries to be loaded #################################
library(dplyr)
library(edgeR)
library(janitor)
library(tidyr)
library(ggplot2)
library(ggthemes)
################################################################################

################### info that always needs to be loaded ########################
DIR="/users/genomics/marta/BLCA"
GENOMEDIR = "/genomics/users/marta/genomes"
biomart = read.csv(file.path(GENOMEDIR,"transcript_gene_v41.txt"))
names(biomart) = c("gene_id","gene_name", "gene_type")

lncRNA = biomart %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene")
nrow(lncRNA)

proteincoding = biomart %>% subset(gene_type == "protein_coding")
nrow(proteincoding)

################################################################################

############################# projects and files ###############################
projects_sergiofolders = c("UC-Genome","Mariathasan","Snyder","UNC-108","HdM")
projects_names = c("UC-GENOME","IMvigor210","SNY-2017","UNC-108","HdM-BLCA-1")
################################################################################


for(i in 1:length(projects_names)) {
  print(projects_names[i])

  # read file
  tableofcounts = read.csv(file.path(DIR,projects_names[i],"analysis/01_counts/TPMs_genenames.csv"))
  
  # print(nrow(tableofcounts))
  # merge with biomart in order to have gene names and gene types as well
  tableofcounts_biomart = merge(tableofcounts, biomart, on="gene_name")
  tableofcounts_biomart = subset(tableofcounts_biomart, !duplicated(gene_name))
  print(nrow(tableofcounts_biomart))
  
  # make it longer
  tableofcounts_biomart = tableofcounts_biomart %>% pivot_longer(cols=-c(gene_id, gene_name, gene_type), names_to = "samples", values_to = "TPM")
  # logtransform
  tableofcounts_biomart$logTPM = log(tableofcounts_biomart$TPM)
  # we only check lncRNA + pseudogenes & protein-coding genes expressions' distribution
  tableofcounts_biomart = tableofcounts_biomart %>% subset(gene_type == "processed_pseudogene" | gene_type == "lncRNA" | gene_type == "protein_coding") %>% mutate(biotype = case_when(gene_type == "processed_pseudogene" ~ "lncRNA",
                                                                                                                                                                            TRUE ~ as.character(gene_type)))
  # avoid redundancies
  tableofcounts_biomart = tableofcounts_biomart %>% unique()
  lncRNA_tableofcounts_biomart = tableofcounts_biomart %>% subset(biotype == "lncRNA") %>% select(gene_name) %>% unique()
  proteincoding_tableofcounts_biomart = tableofcounts_biomart %>% subset(biotype == "protein_coding") %>% select(gene_name) %>% unique()
  
  # plot
  ggplot(tableofcounts_biomart, aes(x=logTPM, fill=biotype)) +
    geom_density(alpha=.5) +
    scale_fill_manual("Gene type", values=c("#006400","#756bb1")) +
    ggtitle(paste0(projects_names[i]," TPM distribution")) +
    geom_vline(xintercept=0) +
    theme_minimal() +
    theme(plot.title = element_text(size=16))
  # ggsave(file.path(DIR,projects_names[i],"/results/TPM_distribution.pdf"), width=5.35, height=2.86)
  # ggsave(file.path(DIR,projects_names[i],"/results/TPM_distribution.png"), width=5.35, height=2.86)
  
  temp = tableofcounts_biomart %>% select(gene_id, biotype) %>% unique() 
  # print(table(temp$biotype))
  
  ## testing cut-off lncRNA/pseudogenes
  TPM1 = tableofcounts_biomart %>% subset(TPM > 1) %>% select(gene_name,biotype) %>% unique()
  print(table(TPM1$biotype))
  TPM1_lncRNA = TPM1 %>% subset(biotype == "lncRNA")
  print(nrow(TPM1_lncRNA) / nrow(lncRNA_tableofcounts_biomart) )*100
  
  TPM1_proteincoding = TPM1 %>% subset(biotype != "lncRNA")
  print(nrow(TPM1_proteincoding) / nrow(proteincoding_tableofcounts_biomart) )*100
  # print(nrow(TPM1_proteincoding) / 20023 )*100
  
}

