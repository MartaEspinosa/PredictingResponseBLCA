#!/usr/bin/env Rscript
#Marta Espinosa
#17/10/2022
#Script to convert read counts into FPKM values

# Run like this:
#Rscript counts2FPKM.R $path_until_project $project_name

################################################
################ Load libraries ################
################################################
#list.of.packages <- c("dplyr")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library("dplyr")
# BiocManager::install("edgeR")
library("edgeR")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Please provide the path until the project name and the project name, in that specific order", call.=FALSE)
}

# Read arguments
GENERAL <- args[1]
PROJECT <- args[2]

################################################
### Read input
################################################

DIR = file.path(GENERAL,PROJECT,"analysis/08_tumor_specific/gffcompare_1to1/healthytissues")
fc = read.csv(file.path(DIR,"healthytissues_featureCounts.header.txt"),sep="\t",header=T) # featureCounts output

################################################

featureLength = fc %>% select(transcript_id, Length)
## From counts to FPKM
y = DGEList(counts=fc[,c(7:ncol(fc))], genes=data.frame(Length=fc[,c(6)]))
y = calcNormFactors(y)
RPKM_fc = as.data.frame(rpkm(y))
RPKM_fc$transcript_id = fc$transcript_id
RPKM_fc$Length = fc$Length
RPKM_fc = RPKM_fc %>% select(transcript_id, Length, everything())

write.csv(RPKM_fc, file=file.path(GENERAL,PROJECT,"analysis/08_tumor_specific/gffcompare_1to1/healthytissues/table_of_counts_TPM.csv"), row.names = FALSE, quote = FALSE) #specify the path and filename for the output to be saved

### is there any novel genes not quantified at all?
# novels = fc %>% subset(grepl("TCONS", transcript_id))
# novels$total_counts = rowSums(novels[,c(8:ncol(novels))])
# length(which(novels$total_counts == 0))
