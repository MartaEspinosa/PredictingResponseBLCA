#!/bin/bash

#SBATCH -p bigmem,long            # Partition to submit to
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu 15Gb           # Memory in MB
#SBATCH -J STAR           # job name
#SBATCH -o /users/genomics/marta/logs/STAR_%A_%a.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/STAR_%A_%a.err    # File to which standard err will be written

####Change output and output error paths to redirect the log files generated

## BUILD VARS AND CREATE FOLDERS
PROJECT=$1
DIR=$2
OUTDIR=$DIR/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files
mkdir -p $CLUSTERDIR

ANNOTGENE=/users/genomics/sergiov/annotations_and_indexes
GNMIDX=/users/genomics/sergiov/annotations_and_indexes/index_50

FASTQDIR=/datasets/marta/TCGA/$PROJECT/fastq_files


#patients_pairs=${GENERAL}/$PROJECT/results/patients.csv
patients_pairs=$2/$PROJECT/results/paired_patients.csv # TCGA DATA


# Read the file, skip the header, and extract the first column
column=$(tail -n +2 "$patients_pairs" | cut -d ',' -f 1)

# Convert the extracted column to an array
IFS=$'\n' read -rd '' -a patients_array <<<"$column"


i=$(($SLURM_ARRAY_TASK_ID -1))

PATIENT=${patients_array[i]}
SAMPLE=${PATIENT}_tumor

echo $PATIENT

module load STAR/2.7.8a-GCC-10.2.0

if [ ! -f  ${OUTDIR}/{$SAMPLE}Aligned.sortedByCoord.out.bai ]; then
  R1=_r1.fastq.gz
  R2=_r2.fastq.gz #this should be commented if the dataset is single-end


  #####################################TUMOR########################################################

  ######################################################################################################
  #####################################ALIGNMENT########################################################

  #gzip fastq files are considered in the code, as well as paired-end reads samples.
  #two pass mode is activated
  #output with uniquely mapped reads ONLY
  STAR --runThreadN $SLURM_CPUS_PER_TASK --limitBAMsortRAM 50000000000 \
   --genomeDir $GNMIDX --readFilesCommand zcat --readFilesIn ${FASTQDIR}/$SAMPLE$R1 ${FASTQDIR}/$SAMPLE$R2 --outFileNamePrefix\
   ${OUTDIR}/$SAMPLE --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFilterType BySJout\
   --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
   --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
   --sjdbGTFfile $ANNOTGENE/gencode.v41.primary_assembly.annotation.gtf --twopassMode Basic

  ######################################################################################################
  #########################################INDEX########################################################


  module purge
  module load SAMtools/1.12-GCC-10.2.0


  samtools index ${OUTDIR}/${SAMPLE}Aligned.sortedByCoord.out.bam ${OUTDIR}/${SAMPLE}Aligned.sortedByCoord.out.bai
fi

#####################################NORMAL########################################################
SAMPLE=${PATIENT}_normal
module load STAR/2.7.8a-GCC-10.2.0
######################################################################################################
#####################################ALIGNMENT########################################################
if [ ! -f  ${OUTDIR}/{$SAMPLE}Aligned.sortedByCoord.out.bai ]; then

  #gzip fastq files are considered in the code, as well as paired-end reads samples.
  #two pass mode is activated
  #output with uniquely mapped reads ONLY
  STAR --runThreadN $SLURM_CPUS_PER_TASK --limitBAMsortRAM 50000000000 \
   --genomeDir $GNMIDX --readFilesCommand zcat --readFilesIn ${FASTQDIR}/$SAMPLE$R1 ${FASTQDIR}/$SAMPLE$R2 --outFileNamePrefix\
   ${OUTDIR}/$SAMPLE --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFilterType BySJout\
   --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
   --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
   --sjdbGTFfile $ANNOTGENE/gencode.v41.primary_assembly.annotation.gtf --twopassMode Basic

  ######################################################################################################
  #########################################INDEX########################################################


  module purge
  module load SAMtools/1.12-GCC-10.2.0


  samtools index ${OUTDIR}/${SAMPLE}Aligned.sortedByCoord.out.bam ${OUTDIR}/${SAMPLE}Aligned.sortedByCoord.out.bai

fi

#scp -rp ${CLUSTERDIR}/${PATIENT}* marta@hydra:$GENERAL/$PROJECT/analysis/05_STAR/uniquely_mapped_2pass_BAM_files
#rm -r ${CLUSTER}/${PATIENT}*

