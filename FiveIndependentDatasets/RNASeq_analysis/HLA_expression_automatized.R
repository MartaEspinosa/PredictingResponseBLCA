# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 10th Oct 2023
# This scripts automatically checks the expression of HLA genes in Responders and NonResponders

######################## libraries to be loaded #################################
library(dplyr)
library(edgeR)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(data.table)
library(ggforce)
library(ggsignif)
################################################################################

################### info that always needs to be loaded ########################
DIR="/users/genomics/marta/BLCA"
GENOMEDIR = "/genomics/users/marta/genomes"
biomart = read.csv(file.path(GENOMEDIR,"transcript_gene_v41.txt"))
names(biomart) = c("gene_id","gene_name", "gene_type")
################################################################################

############################# projects and files ###############################
projects_names = c("UC-GENOME","IMvigor210","SNY-2017","UNC-108","HdM-BLCA-1")
wilcoxon_test_biomarkers_file = file.path(DIR,"RNASeq_expression/HLAgenes/HLA_wilcoxontest.txt")
metaanalysis_biomarkers_file = file.path(DIR,"RNASeq_expression/HLAgenes/HLA_metaanalysis.txt")
################################################################################

############################# selected biomarkers ##############################
## only protein_coding ones
HLA_genes = biomart %>% subset(grepl("^HLA", gene_name) & gene_type == "protein_coding") %>% select(gene_id, gene_name)


HLA_expr = data.frame("gene_name" = character(),
                      "patient" = character(),
                      "TPM" = numeric(),
                      "logTPM" = numeric(),
                      "Dataset" = factor())
HLA_zscores = data.frame("gene_name" = character(),
                         "patient" = character(),
                         "TPM" = numeric(),
                         "logTPM" = numeric(),
                         "Dataset" = factor())
HLA_pvalues = data.frame("gene_name" = character(),
                         "Dataset" = factor(),
                         "Wilcox-pvalues" = numeric())
HLA_forest = data.frame("gene_name" = character(),
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
  biomarkers_expression = tableofcounts %>% subset(gene_name %in% HLA_genes$gene_name)
  biomarkers_expression %>% head
  biomarkers_expression = biomarkers_expression %>% pivot_longer(cols=-c("gene_name"), names_to = "patient", values_to = "TPM")
  
  zscores_selected = zscores %>% subset(gene_name %in% HLA_genes$gene_name)
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
  # compute logTPM 
  biomarkers_expression$logTPM <- log10(biomarkers_expression$TPM)
  biomarkers_expression$Dataset = projects_names[i]
  zscores_selected$Dataset = projects_names[i]
  
  HLA_expr = rbind(HLA_expr, biomarkers_expression)
  HLA_zscores = rbind(HLA_zscores,zscores_selected)
  
  ######################## plot per dataset #############################
  ggplot(biomarkers_expression, aes(x=Response, y=logTPM)) +
    geom_boxplot(aes(fill=Response),width=0.5) +
    ggtitle(projects_names[i]) +
    scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size=16)) +
    facet_wrap(~gene_name, scales="free")
  ggsave(file.path(DIR,projects_names[i],"/results/HLAgenes_expression.pdf"), width=6.85, height=6.15)
  ggsave(file.path(DIR,projects_names[i],"/results/HLAgenes_expression.png"), width=6.85, height=6.15)
  
  # compute a Wilcoxon test to see if there is a statistical association between the expression of each biomarker with the response to immunotherapy
  for(r in 1:length(unique(HLA_expr$gene_name))) {
    print(projects_names[i])
    gene_subset = HLA_expr[HLA_expr$gene_name == unique(HLA_expr$gene_name)[r],]
    gene_subset_test = wilcox.test(TPM ~ Response, data=gene_subset)
    
    temp = data.frame("gene_name" = c(unique(HLA_expr$gene_name[r])),
                      "Dataset" = c(projects_names[i]),
                      "Wilcox-pvalues" = c(gene_subset_test$p.value))
    HLA_pvalues = rbind(HLA_pvalues, temp)
  }
}
# save wilcoxon p values in a text file
write.csv(HLA_pvalues, file=wilcoxon_test_biomarkers_file, row.names = F, quote = F) 


######################## plot per HLA gene #############################
p <- ggplot(HLA_expr, aes(x=Dataset, y=logTPM)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("HLA genes expression") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap_paginate(~ gene_name, scales="free", ncol = 1, nrow = 1)
for(a in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ gene_name, ncol = 1, nrow = 1, page = a)
  ggsave(plot = p_save, filename = paste0(DIR,'/RNASeq_expression/HLAgenes/HLAgenes_page', a, '.png'), width=6.99, height=5.85)
  ggsave(plot = p_save, filename = paste0(DIR,'/RNASeq_expression/HLAgenes/HLAgenes_page', a, '.pdf'), width=6.99, height=5.85)

}

temp = HLA_expr %>% select(patient, Response, Dataset) %>% unique()
table(temp$Response, temp$Dataset)


######################## Z-SCORES #############################
HLA_zscores %>% head
table(HLA_zscores$Dataset)

############### plot per HLA gene w/ zscores ##################
ggplot(HLA_zscores, aes(x=Dataset, y=zscore)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("HLA genes expression") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap_paginate(~ gene_name, scales="free", ncol = 1, nrow = 1)
for(a in 1:n_pages(p)){
  p_save <-  p + 
    facet_wrap_paginate(~ gene_name, ncol = 1, nrow = 1, page = a)
  ggsave(plot = p_save, filename = paste0(DIR,'/RNASeq_expression/HLAgenes/HLAgenes_zscore_page', a, '.png'), width=6.99, height=5.85)
  ggsave(plot = p_save, filename = paste0(DIR,'/RNASeq_expression/HLAgenes/HLAgenes_zscore_page', a, '.pdf'), width=6.99, height=5.85)
  
}

# compute a Wilcoxon test to see if there is a statistical association between the expression of each biomarker with the response to immunotherapy
HLA_zscore_pvalues = data.frame("gene_name" = character(),
                                "Wilcox-pvalues" = numeric())

for(r in 1:length(unique(HLA_zscores$gene_name))) {
  print(projects_names[i])
  gene_subset = HLA_zscores[HLA_zscores$gene_name == unique(HLA_zscores$gene_name)[r],]
  gene_subset_test = wilcox.test(zscore ~ Response, data=gene_subset)
  
  temp = data.frame("gene_name" = c(unique(HLA_zscores$gene_name[r])),
                    "Wilcox-pvalues" = c(gene_subset_test$p.value))
  print(temp)
  HLA_zscore_pvalues = rbind(HLA_zscore_pvalues, temp)
}

print(HLA_zscore_pvalues)

############### general plot with z-scores ##################
ggplot(HLA_zscores, aes(x=Response, y=zscore)) +
  geom_boxplot(aes(fill=Response),width=0.5) +
  ggtitle("HLA genes expression") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12, margin = margin(b = 15))) +
  facet_wrap(~ gene_name, scales="free") +
  geom_label(data = HLA_zscore_pvalues, aes(label = paste("p =", round(Wilcox.pvalues, 3))),
            x = Inf, y = Inf, hjust = 2, vjust =1, size = 3,
            fill = "white", color = "black")
ggsave(file.path(DIR,"RNASeq_expression/HLAgenes/HLAgenes_expression_zscores_420_pvalues.pdf"), width=8.23, height=7.35)
ggsave(file.path(DIR,"RNASeq_expression/HLAgenes/HLAgenes_expression_zscores_420_pvalues.png"), width=8.23, height=7.35)




############################# META ANALYSIS ##############################
print("Starting Metaanalysis")
library(meta)

for(u in 1:length(unique(HLA_expr$gene_name))) {
  # select and compute needed info for Responders
  R_to_metaanalysis = HLA_expr %>% 
    subset(gene_name == unique(HLA_expr$gene_name)[u]) %>% 
    subset(Response == "R") %>%
    group_by(Dataset) %>% 
    mutate(mean.R=mean(TPM), sd.R = sd(TPM)) %>%
    mutate(n.R = length(Response == "R")) %>%
    select(Dataset, n.R, mean.R, sd.R) %>% 
    unique()
  
  # select and compute needed info for Non-Responders
  NR_to_metaanalysis = HLA_expr %>% 
    subset(gene_name == unique(HLA_expr$gene_name)[u]) %>% 
    subset(Response == "NR") %>%
    group_by(Dataset) %>% 
    mutate(mean.NR=mean(TPM), sd.NR = sd(TPM)) %>%
    mutate(n.NR = length(Response == "NR")) %>%
    select(Dataset, n.NR, mean.NR, sd.NR) %>% 
    unique()
  
  # merge both
  subset_to_metaanalysis = merge(NR_to_metaanalysis, R_to_metaanalysis, by="Dataset")
  
  # calculate the meta analysis
  res =  metacont(n.R, mean.R, sd.R,
                      n.NR, mean.NR, sd.NR,
                      comb.fixed = T, comb.random = T, studlab = Dataset,
                      data = subset_to_metaanalysis, sm = "SMD")
  print(res)
  
  # save the values in a dataframe to later save them in a text file
  temp_forest = data.frame("gene_name" = c(unique(HLA_expr$gene_name[u])),
                    "Heterogeneity" = res$pval.Q,
                    "Randomeffects-pvalues" = res$pval.random,
                    "Commoneffects-pvalues" = res$pval.common)
  print(temp_forest)
  HLA_forest = rbind(HLA_forest, temp_forest)
  
  ############### forest plot per dataset ##################
  svg(file=paste0(DIR,'/RNASeq_expression/HLAgenes/forestplot_',unique(HLA_expr$gene_name)[u],'.svg'), width=12) # Open PDF device with specific file name
  forest(res, leftcols = c('Dataset'))
  dev.off()
}

# save metaanalysis p values in a text file
write.csv(HLA_forest, file=metaanalysis_biomarkers_file, row.names = F, quote = F) 


############################# HEATMAP ##############################
HLA_zscores %>% head

# compute mean z-score per HLA gene and dataset
df_zscores = HLA_zscores %>% select(gene_name, zscore, Dataset, Response) %>% group_by(gene_name, Response, Dataset) %>% summarise(mean_zscore=mean(zscore))
df_zscores %>% head
df_zscores = df_zscores %>% pivot_wider(names_from = "Response", values_from = "mean_zscore")

# compute the fold change
df_zscores$FC = df_zscores$R - df_zscores$NR
df_zscores_dataset = df_zscores %>% select(-c(R, NR)) %>% pivot_wider(names_from = Dataset, values_from = FC)
df_zscores_dataset = as.data.frame(df_zscores_dataset)
rownames(df_zscores_dataset) = df_zscores_dataset$gene_name
# make two groups: class I and class II
HLAI = df_zscores_dataset %>% subset(nchar(gene_name) == 5)
HLAII = df_zscores_dataset %>% subset(nchar(gene_name) > 5)

HLAI$gene_name = NULL
HLAII$gene_name = NULL

# function to save the heatmap
save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


breaksList = seq(-1.5, 1.5, by=.05)
library(pheatmap)

##### HLA-II
matrix_HLAII = data.matrix(HLAII)
plot = pheatmap(matrix_HLAII, cluster_cols = F,cluster_rows = T,
         treeheight_col = 0,
         display_numbers = F,
         breaks = breaksList,
         color = colorRampPalette(c("#1e375f","#2a6aa9","#95c3d7","#deebef","#f9fbfa","#f9dacb","#eda389","#d55e53","#a71530","#670220"))(length(breaksList)),
         main="HLA expression Responders/NonResponders")
save_pheatmap_pdf(plot, file.path(DIR,"RNASeq_expression/HLAgenes/heatmap_HLAII.pdf"), 6.10, 7.35)

##### HLA-I
matrix_HLAI = data.matrix(HLAI)
plot = pheatmap(matrix_HLAI, cluster_cols = F,cluster_rows = T,
         treeheight_col = 0,
         display_numbers = F,
         breaks = breaksList,
         color = colorRampPalette(c("#1e375f","#2a6aa9","#95c3d7","#deebef","#f9fbfa","#f9dacb","#eda389","#d55e53","#a71530","#670220"))(length(breaksList)),
         main="HLA expression Responders/NonResponders")
save_pheatmap_pdf(plot, file.path(DIR,"RNASeq_expression/HLAgenes/heatmap_HLAI.pdf"), 6.10, 4.4)



#### statistical test for heatmap - t.test
HLA_ttest = data.frame("gene_name" = character(),
                       "Dataset" = factor(),
                                "pvalues" = numeric())
total_HLA_expression = data.frame("Dataset" = factor(),
                                  "pvalues" = numeric())

for(i in 1:length(projects_names)) {
  print(projects_names[i])
  
  for(r in 1:length(unique(HLA_zscores$gene_name))) {
    print(unique(HLA_zscores$gene_name)[r])
    # subset per dataset and gene
    gene_subset = HLA_zscores %>% subset(Dataset == projects_names[i]) %>% subset(gene_name == unique(HLA_zscores$gene_name)[r])
    # run the chosen test (t.test)
    gene_subset_test = t.test(zscore ~ Response, data=gene_subset)
    
    temp = data.frame("gene_name" = c(unique(HLA_zscores$gene_name[r])),
                      "Dataset" = projects_names[i],
                      "pvalues" = c(gene_subset_test$p.value))
    HLA_ttest = rbind(HLA_ttest, temp)
    
  }
  # subset only per dataset to have the whole effect of the HLA genes per dataset
  total_subset = HLA_zscores %>% subset(Dataset == projects_names[i])
  # run the chosen test (t.test)
  total_subset_test = t.test(zscore ~ Response, data=total_subset)
  temp_total = data.frame("Dataset" = projects_names[i],
                    "pvalues" = c(total_subset_test$p.value))
  total_HLA_expression = rbind(total_HLA_expression, temp_total)
}
print(HLA_ttest)

# save heatmap p values from t.test in a text file
HLA_ttest = HLA_ttest %>% pivot_wider(names_from = "Dataset", values_from = "pvalues")                                     
write.csv(HLA_ttest, file.path(DIR,"RNASeq_expression/HLAgenes/ttest_zscores.csv"), quote=F, row.names = F)

total_HLA_expression


######## strategies to reduce the number of variables regarding HLA genes ######### 
###################### mean ###################### 
HLA_zscores = HLA_zscores %>% mutate(class = case_when(nchar(gene_name) == 5 ~ "HLA-I",
                                                       nchar(gene_name) > 5 ~ "HLA-II"))
mean_HLA_zscores = HLA_zscores %>% group_by(patient, Response, Dataset, class) %>% summarise(mean_HLA_zscore = mean(zscore))
mean_HLA_zscores %>% head
mean_HLA_zscores_wide = mean_HLA_zscores %>% ungroup() %>% select(-c(Response, Dataset))
mean_HLA_zscores_wide = mean_HLA_zscores_wide %>% pivot_wider(names_from = "class", values_from = "mean_HLA_zscore")
mean_HLA_zscores_wide %>% head
mean_HLA_zscores_wide = as.data.frame(mean_HLA_zscores_wide)
row.names(mean_HLA_zscores_wide) = mean_HLA_zscores_wide$patient
mean_HLA_zscores_wide$patient = NULL
write.csv(mean_HLA_zscores_wide, file.path(DIR,"RNASeq_expression/HLAgenes/mean_HLAgenes.csv"), quote=F)

ggplot(mean_HLA_zscores %>% subset(class == "classI"), aes(x=patient, y=mean_HLA_zscore, fill=Response)) +
  geom_histogram(stat="identity") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  ggtitle("Average HLA-I genes expression") +
  labs(x="Patients",
       y="Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap(Dataset ~ Response, scales = "free", ncol = 2)

ggplot(mean_HLA_zscores %>% subset(class == "classII"), aes(x=patient, y=mean_HLA_zscore, fill=Response)) +
  geom_histogram(stat="identity") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  ggtitle("Average HLA-II genes expression") +
  labs(x="Patients",
       y="Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap(Dataset ~ Response, scales = "free", ncol = 2)


mean_of_means_HLA_zscores = mean_HLA_zscores %>% group_by(Dataset, Response, class) %>% summarise(mean_of_means = mean(mean_HLA_zscore))
mean_of_means_HLA_zscores_wide = mean_of_means_HLA_zscores %>% pivot_wider(names_from = "Response", values_from = "mean_of_means")
mean_of_means_HLA_zscores_wide %>% head

ggplot(mean_of_means_HLA_zscores, aes(x=class, y=mean_of_means, fill=Response)) +
  geom_col(position="dodge") +
  scale_fill_manual(values=c("#5ab4ac","#d8b365")) +
  ggtitle("HLA genes expression") +
  labs(x="Class HLA genes",
       y="Z-score") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(size=16),
        strip.text = element_text(size=12)) +
  facet_wrap( ~ Dataset, scales = "free", ncol = 2)
ggsave(file.path(DIR,"RNASeq_expression/HLAgenes/meanHLAgenes_zscores.png"), width=8.27, height=8.03)
ggsave(file.path(DIR,"RNASeq_expression/HLAgenes/meanHLAgenes_zscores.pdf"), width=8.27, height=8.03)


###################### GSVA ###################### 
library(GSVA)

## all of them
HLA_zscores_gsva = HLA_zscores %>% select(gene_name, patient, zscore) %>% pivot_wider(names_from = "patient", values_from = "zscore")
HLA_zscores_gsva = as.data.frame(HLA_zscores_gsva)
# make rownames as gene_names
row.names(HLA_zscores_gsva) = c(HLA_zscores_gsva$gene_name)
HLA_genes = as.vector(HLA_zscores_gsva$gene_name)
HLA_zscores_gsva$gene_name = NULL

## HLA-I
HLAI = HLA_zscores %>% subset(class == "HLA-I")
HLAI_zscores_gsva = HLAI %>% select(gene_name, patient, zscore) %>% pivot_wider(names_from = "patient", values_from = "zscore")
HLAI_zscores_gsva = as.data.frame(HLAI_zscores_gsva)
# make rownames as gene_names
row.names(HLAI_zscores_gsva) = c(HLAI_zscores_gsva$gene_name)
HLAI_genes = as.vector(HLAI_zscores_gsva$gene_name)
HLAI_zscores_gsva$gene_name = NULL

## HLA-II
HLAII = HLA_zscores %>% subset(class == "HLA-II")
HLAII_zscores_gsva = HLAII %>% select(gene_name, patient, zscore) %>% pivot_wider(names_from = "patient", values_from = "zscore")
HLAII_zscores_gsva = as.data.frame(HLAII_zscores_gsva)
# make rownames as gene_names
row.names(HLAII_zscores_gsva) = c(HLAII_zscores_gsva$gene_name)
HLAII_genes = as.vector(HLAII_zscores_gsva$gene_name)
HLAII_zscores_gsva$gene_name = NULL

gs = list("HLA-I" = HLAI_genes,
          "HLA-II" = HLAII_genes)
gs
gsva_HLA = gsva(as.matrix(HLA_zscores_gsva),gs,method="gsva", verbose=T)

write.csv(t(gsva_HLA), file.path(DIR,"RNASeq_expression/HLAgenes/GSVA_scores.csv"), quote=F)
