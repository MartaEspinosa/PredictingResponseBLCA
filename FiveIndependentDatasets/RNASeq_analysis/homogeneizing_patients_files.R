# Author: Marta Espinosa Camarena
# email: mespinosa@imim.es
# Date: 28th Sept 2023
# This scripts homogeneizes patients file data for RNA-Seq datasets of BLCA study

######################## libraries to be loaded #################################
library(dplyr)
library(tidyr)
library(data.table)
################################################################################

############################### UC-GENOME ######################################

UCG_patients = read.csv("/datasets/sergio/UC-GENOME/4_Results/response_index.csv") # only RNA and no outliers
lUCG_patients = read.csv("/datasets/sergio/UC-GENOME/mapper.csv") # only RNA and no outliers
# UCG_patients = read.csv("/genomics/users/lilian/Comparing_Cohorts/clinicalData_compareGroups/UCgenome/UC-GENOME_Immunotherapy_clinicaldata.csv")
# UCG_patients = UCG_patients %>% select(submitted_subject_id, Response, RECIST) %>% subset(RECIST == "PR" | RECIST == "CR" | RECIST == "PD") %>% drop_na()
UCG_patients$X = NULL
names(UCG_patients) = c("patient","Response")
table(UCG_patients$Response)
write.csv(UCG_patients, "/users/genomics/marta/BLCA/UC-GENOME/patients_response.csv", quote=F, row.names = F)
UCG_patients$Dataset = "UC-GENOME"

merged = merge(UCG_patients, lUCG_patients, by.x = "patient", by.y="Sample.Name")
merged %>% head
merged = merged %>% select(patient, Run, Response)
write.csv(merged, "/users/genomics/marta/BLCA/UC-GENOME/patients_run_response.csv", quote=F, row.names = F)

############################### IMvigor210 #####################################
# MHS_patients = read.csv("/datasets/sergio/Mariathasan/4_Results/tmb_response_tcga.csv")
MHS_patients = read.csv("/genomics/users/lilian/download_data/EGA_Mariathasan/IMvigor210CoreBiologies_2.0.0/Mariathasan_FMOne_Response.csv")
MHS_patients = MHS_patients%>% select(Patient.ID, Response) %>% drop_na()
names(MHS_patients) = c("patient","Response")
table(MHS_patients$Response)
write.csv(MHS_patients, "/users/genomics/marta/BLCA/IMvigor210/patients_response.csv", quote=F, row.names = F)
MHS_patients$Dataset = "IMvigor210"

############################## SNY-2017 #####################################
identifiers = read.csv("/genomics/users/lilian/download_data/dbGaP.proj32726/rna_paths_sergiov.csv")
SNY_patients = read.csv(file.path("/projects_eg/projects/lilian/Snyder/Snyder-BLCA_clinicalData.csv"))
SNY_patients = merge(SNY_patients, identifiers, by.x ="Tumor_Sample_Barcode", by.y = c("Patient.ID"))
SNY_patients = SNY_patients %>% select(Run, index, Tumor_Sample_Barcode, Response)
# SNY_patients = SNY_patients %>% mutate(Response = case_when(Best.Response.RECIST.1.1 == "PD" ~ "NR",
#                                                             Best.Response.RECIST.1.1 == "CR" ~ "R",
#                                                             Best.Response.RECIST.1.1 == "PR" ~ "R"))
SNY_patients = SNY_patients %>% subset((Response == "NR") | (Response == "R")) %>% drop_na()
names(SNY_patients)[3] = "patient"
SNY_patients = SNY_patients %>% select(index, patient, Response)
SNY_patients$patient = gsub("_T","", SNY_patients$patient)
write.csv(SNY_patients, "/users/genomics/marta/BLCA/SNY-2017/patients_index_response.csv", quote=F, row.names = F)

SNY_patients = SNY_patients %>% select(patient, Response)
table(SNY_patients$Response)
write.csv(SNY_patients, "/users/genomics/marta/BLCA/SNY-2017/patients_response.csv", quote=F, row.names = F)
SNY_patients$Dataset = "SNY-2017"

############################## HdM-BLCA-1 #####################################
# should be 20, check Lili's paper, Julia used 20
HM_patients = read.csv(file.path("/users/genomics/marta/BLCA/HdM-BLCA-1/metadata_RNAseq_sergio.csv")) #R16 fuera
# HM_patients = fread("/projects_eg/projects/lilian/HdM/HdM-BLCA-1_clinicalData.csv") # Lili's
HM_patients = HM_patients %>% subset(grepl("RNA",NGS))
HM_patients = HM_patients%>% select(2,3) 
names(HM_patients) = c("patient","Response")
table(HM_patients$Response)
write.csv(HM_patients, "/users/genomics/marta/BLCA/HdM-BLCA-1/patients_response.csv", quote=F, row.names = F)
HM_patients$Dataset = "HdM-BLCA-1" 

################################# UNC-108 ######################################
# UNC_patients = read.csv(file.path("/datasets/sergio/UNC-108/4_Results/unc-108_metadata.csv"))
UNC_patients = read.csv("/genomics/users/lilian/Comparing_Cohorts/clinicalData_compareGroups/UNC-108/UNC-108_clinicalData_clean_totalTMB.csv")
UNC_patients = UNC_patients%>% select(Tumor_Sample_Barcode, Response) %>% drop_na()
names(UNC_patients) = c("patient","Response")
UNC_patients$patient = gsub("_.*","", UNC_patients$patient)
UNC_patients = UNC_patients %>% unique()
# UNC_patients = UNC_patients %>% mutate(Response = case_when(response == "PD" ~ "NR",
#                                                             response == "CR" ~ "R",
#                                                             response == "PR" ~ "R")) %>% drop_na() %>% select(-response)

table(UNC_patients$Response)
write.csv(UNC_patients, "/users/genomics/marta/BLCA/UNC-108/patients_response.csv", quote=F, row.names = F)
UNC_patients$Dataset = "UNC-108" 

patients = rbind(UNC_patients, HM_patients, MHS_patients, UCG_patients,SNY_patients)
patients
write.csv(patients, "/users/genomics/marta/BLCA/RNASeq_patients.csv", quote=F, row.names = F)
table(patients$Dataset)
