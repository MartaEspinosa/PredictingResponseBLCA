# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 9th Oct 2023
# This scripts automatically checks the expression of PD1 and PDL1 genes in Responders and NonResponders

######################## libraries to be loaded #################################
library(dplyr)
library(edgeR)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(data.table)

################################################################################

################### info that always needs to be loaded ########################
DIR="/users/genomics/marta/BLCA"
GENOMEDIR = "/genomics/users/marta/genomes"
# biomart = read.csv(file.path(GENOMEDIR,"transcript_gene_v41.txt"))
# names(biomart) = c("gene_id","gene_name", "gene_type")
################################################################################

############################# projects and files ###############################
projects_names = c("UC-GENOME","IMvigor210","SNY-2017","UNC-108","HdM-BLCA-1")
wilcoxon_test_biomarkers_file = file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1PDL1_wilcoxontest.txt")
metaanalysis_biomarkers_file = file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1PDL1_metaanalysis.txt")
################################################################################

############################# selected biomarkers ##############################
### PD1 - PDCD1 - ENSG00000188389
### PDL1 - CD274 - ENSG00000120217

PDs = data.frame("gene_id" = c('ENSG00000188389', 'ENSG00000120217'),
                 "gene_name" = c("PDCD1","CD274"))
PDs_expr = data.frame("gene_name" = character(),
                                 "patient" = character(),
                                 "TPM" = numeric(),
                                 "logTPM" = numeric(),
                                 "Dataset" = factor())
PDs_zscores = data.frame("gene_name" = character(),
                      "patient" = character(),
                      "TPM" = numeric(),
                      "logTPM" = numeric(),
                      "Dataset" = factor())
PDs_pvalues = data.frame("gene_name" = character(),
                                "Dataset" = factor(),
                                "Wilcox-pvalues" = numeric())
PD_forest = data.frame("gene_name" = character(),
                        "Heterogeneity" = numeric(),
                        "Randomeffects-pvalues" = numeric(),
                        "Commoneffects-pvalues" = numeric())
################################################################################

for(i in 1:length(projects_names)) {
  print(projects_names[i])
  patient_file = read.csv(file.path(DIR,projects_names[i],"patients_response.csv"))

  # read file
  tableofcounts = read.csv(file.path(DIR,projects_names[i],"analysis/01_counts/TPMs_genenames.csv"))
  zscores = fread(file.path(DIR,projects_names[i],"analysis/01_counts/zscores_geneasrows_patientsascols.csv"), sep=",", header=T)
  
  # subset biomarkers expression
  biomarkers_expression = tableofcounts %>% subset(gene_name %in% PDs$gene_name)
  biomarkers_expression %>% head
  biomarkers_expression = biomarkers_expression %>% pivot_longer(cols=-c("gene_name"), names_to = "patient", values_to = "TPM")
  
  zscores_selected = zscores %>% subset(gene_name %in% PDs$gene_name)
  zscores_selected = zscores_selected %>% pivot_longer(cols=-c("gene_name"), names_to = "patient", values_to = "zscore")
  
  # correct patients IDs if needed
  biomarkers_expression$patient = gsub("X","", biomarkers_expression$patient)
  zscores_selected$patient = gsub("X","", zscores_selected$patient)
  
  if(projects_names[i] == "HdM-BLCA-1") {
    biomarkers_expression$patient = gsub("_.*","", biomarkers_expression$patient)
    zscores_selected$patient = gsub("_.*","", zscores_selected$patient)
    
  }
  # merge data with patients information to know whether a patient is R or NR
  biomarkers_expression = merge(biomarkers_expression, patient_file, by="patient")
  zscores_selected = merge(zscores_selected, patient_file, by="patient")
  biomarkers_expression$logTPM <- log10(biomarkers_expression$TPM)
  biomarkers_expression$Dataset = projects_names[i]
  zscores_selected$Dataset = projects_names[i]
  
  PDs_expr = rbind(PDs_expr, biomarkers_expression)
  PDs_zscores = rbind(PDs_zscores,zscores_selected)
  ggplot(biomarkers_expression, aes(x=Response, y=logTPM)) +
    geom_boxplot(aes(fill=Response),width=0.5) +
    ggtitle(projects_names[i]) +
    scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size=16),
          strip.text = element_text(size=12)) +
    facet_wrap(~gene_name, scales="free",nrow=3, ncol=1)
  ggsave(file.path(DIR,projects_names[i],"/results/PD1-PDL1.pdf"), width=4.06, height=6.15)
  ggsave(file.path(DIR,projects_names[i],"/results/PD1-PDL1.png"), width=4.06, height=6.15)
  
  # compute a Wilcoxon test to see if there is a statistical association between the expression of each biomarker with the response to immunotherapy
  PDL1 = biomarkers_expression[biomarkers_expression$gene_name == "CD274",]
  PDL1_test = wilcox.test(TPM ~ Response, data=PDL1)
  
  PD1 = biomarkers_expression[biomarkers_expression$gene_name == "PDCD1",]
  PD1_test = wilcox.test(TPM ~ Response, data=PD1)
  
  temp_tests = data.frame("gene_name" = c("PD1","PDL1"),
                          "Dataset" = c(projects_names[i],projects_names[i]),
                          "Wilcox-pvalues" = c(PDL1_test$p.value,PD1_test$p.value))
  PDs_pvalues = rbind(PDs_pvalues, temp_tests)
}
# save wilcoxon p values in a text file
write.csv(PDs_pvalues, file=wilcoxon_test_biomarkers_file, row.names = F, quote = F) 

PDs_expr = PDs_expr %>% mutate(gene_name_PD = case_when(gene_name == "CD274" ~ "PDL1",
                                                        gene_name == "PDCD1" ~ "PD1"))
ggplot(PDs_expr, aes(x=Dataset, y=logTPM)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("PD1 and PDL1 expression") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap(~gene_name_PD, scales="free",nrow=3, ncol=1)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1.png"), width=6.99, height=5.85)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1.pdf"), width=6.99, height=5.85)

temp = PDs_expr %>% select(patient, Response, Dataset) %>% unique()
table(temp$Response, temp$Dataset)


#### zscores
PDs_zscores %>% head
PDs_zscores = PDs_zscores %>% mutate(gene_name_PD = case_when(gene_name == "CD274" ~ "PDL1",
                                                        gene_name == "PDCD1" ~ "PD1"))
ggplot(PDs_zscores, aes(x=Dataset, y=zscore)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("Expression of PD-1 and PD-L1") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap( ~ gene_name_PD, scales="free", ncol=1)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1_zscore.png"), width=6.99, height=5.85)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1_zscore.pdf"), width=6.99, height=5.85)

ggplot(PDs_zscores, aes(x=Response, y=zscore)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("Expression of PD-1 and PD-L1") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap( ~ gene_name_PD, scales="free", ncol=1)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1_zscore_420.png"), width=6.99, height=5.85)
ggsave(file.path(DIR,"RNASeq_expression/PD1-PDL1/PD1_PDL1_zscore_420.pdf"), width=6.99, height=5.85)

zscores_PD1_values = PDs_zscores%>% subset(gene_name == "PDCD1") %>% select(zscore)
shapiro.test(zscores_PD1_values$zscore)

zscores_PDL1_values = PDs_zscores%>% subset(gene_name == "CD274") %>% select(zscore)
shapiro.test(zscores_PDL1_values$zscore)

PD1_zscore = PDs_zscores[PDs_zscores$gene_name == "PDCD1",]
PD1_zscore_test = wilcox.test(zscore ~ Response, data=PD1_zscore)
PD1_zscore_test

PDL1_zscore = PDs_zscores[PDs_zscores$gene_name == "CD274",]
PDL1_zscore_test = wilcox.test(zscore ~ Response, data=PDL1_zscore)
PDL1_zscore_test

####### METAANALYSIS
print("Starting Metaanalysis")
library(meta)
PDs_expr %>% head

for(u in 1:length(unique(PDs_expr$gene_name))) {
  R_to_metaanalysis = PDs_expr %>% 
    subset(gene_name == unique(PDs_expr$gene_name)[u]) %>% 
    subset(Response == "R") %>%
    group_by(Dataset) %>% 
    mutate(mean.R=mean(TPM), sd.R = sd(TPM)) %>%
    mutate(n.R = length(Response == "R")) %>%
    select(Dataset, n.R, mean.R, sd.R) %>% 
    unique()
  
  NR_to_metaanalysis = PDs_expr %>% 
    subset(gene_name == unique(PDs_expr$gene_name)[u]) %>% 
    subset(Response == "NR") %>%
    group_by(Dataset) %>% 
    mutate(mean.NR=mean(TPM), sd.NR = sd(TPM)) %>%
    mutate(n.NR = length(Response == "NR")) %>%
    select(Dataset, n.NR, mean.NR, sd.NR) %>% 
    unique()
  
  subset_to_metaanalysis = merge(NR_to_metaanalysis, R_to_metaanalysis, by="Dataset")
  
  res =  metacont(n.R, mean.R, sd.R,
                  n.NR, mean.NR, sd.NR,
                  comb.fixed = T, comb.random = T, studlab = Dataset,
                  data = subset_to_metaanalysis, sm = "SMD")
  print(res)
  
  temp_forest = data.frame("gene_name" = c(unique(PDs_expr$gene_name[u])),
                           "Heterogeneity" = res$pval.Q,
                           "Randomeffects-pvalues" = res$pval.random,
                           "Commoneffects-pvalues" = res$pval.common)
  print(temp_forest)
  PD_forest = rbind(PD_forest, temp_forest)
  
  svg(file=paste0(DIR,'/RNASeq_expression/PD1-PDL1/forestplot_',unique(PDs_expr$gene_name)[u],'.svg'), width=12) # Open PDF device with specific file name
  forest(res, leftcols = c('Dataset'))
  dev.off()
}

# save metaanalysis p values in a text file
write.csv(PD_forest, file=metaanalysis_biomarkers_file, row.names = F, quote = F) 




PD1_R = PDs_expr %>% 
  subset(gene_name_PD == "PD1") %>% subset(Response == "R") %>%
  group_by(Dataset) %>% 
  mutate(mean.R=mean(TPM), sd.R = sd(TPM)) %>%
  mutate(n.R = length(Response == "R")) %>%
  select(Dataset, n.R, mean.R, sd.R) %>% 
  unique()

PD1_NR = PDs_expr %>% 
  subset(gene_name_PD == "PD1") %>% subset(Response == "NR") %>%
  group_by(Dataset) %>% 
  mutate(mean.NR=mean(TPM), sd.NR = sd(TPM)) %>%
  mutate(n.NR = length(Response == "NR")) %>%
  select(Dataset, n.NR, mean.NR, sd.NR) %>% 
  unique()
PD1 = merge(PD1_NR, PD1_R, by="Dataset")

# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.PD1 =  metacont(n.R, mean.R, sd.R,
                    n.NR, mean.NR, sd.NR,
                    comb.fixed = T, comb.random = T, studlab = Dataset,
                    data = PD1, sm = "SMD")
res.PD1

svg(file=file.path(DIR,'RNASeq_expression/PD1-PDL1/forestplot_PD1.svg'), width=12) # Open PDF device with specific file name
forest(res.PD1, leftcols = c('Dataset'))
dev.off()


PDL1_R = PDs_expr %>% 
  subset(gene_name_PD == "PDL1") %>% subset(Response == "R") %>%
  group_by(Dataset) %>% 
  mutate(mean.R=mean(TPM), sd.R = sd(TPM)) %>%
  mutate(n.R = length(Response == "R")) %>%
  select(Dataset, n.R, mean.R, sd.R) %>% 
  unique()

PDL1_NR = PDs_expr %>% 
  subset(gene_name_PD == "PDL1") %>% subset(Response == "NR") %>%
  group_by(Dataset) %>% 
  mutate(mean.NR=mean(TPM), sd.NR = sd(TPM)) %>%
  mutate(n.NR = length(Response == "NR")) %>%
  select(Dataset, n.NR, mean.NR, sd.NR) %>% 
  unique()
PDL1 = merge(PDL1_NR, PDL1_R, by="Dataset")

# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.PDL1 =  metacont(n.R, mean.R, sd.R,
                     n.NR, mean.NR, sd.NR,
                     comb.fixed = T, comb.random = T, studlab = Dataset,
                     data = PDL1, sm = "SMD")
res.PDL1

svg(file=file.path(DIR,'RNASeq_expression/PD1-PDL1/forestplot_PDL1.svg'), width=12) # Open PDF device with specific file name
forest(res.PDL1, leftcols = c('Dataset'))
dev.off()



