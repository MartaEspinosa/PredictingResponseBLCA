#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu 16Gb           # Memory in MB
#SBATCH -J featCounts_gff           # job name
#SBATCH -o /users/genomics/marta/logs/featCounts_gff.%j.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/featCounts_gff.%j.err    # File to which standard err will be written

####Change output and output error paths to redirect the log files generated

###PREPARING NEEDED DATA
PROJECT=$1
DIR=$2/$PROJECT
p=$3 #single-end or paired-end dataset
strand=$4
CLUSTERDIR=$2/$PROJECT/analysis/07_quantification/straightforward # neither stringtie nor gffcompare has been run
mkdir -p $CLUSTERDIR
AnnotGTF=/users/genomics/sergiov/annotations_and_indexes/gencode.v41.primary_assembly.annotation.gtf

module load Subread/2.0.3
########################

# countReadPairs may need to be removed in case of single-end reads

if [ $p == "paired-end" ]; then
    if [ $strand == "firststrand" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -p -s 2 -g gene_name -O --countReadPairs -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    elif [ $strand == "secondstrand" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -p -s 1 -g gene_name -O --countReadPairs -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    elif [ $strand == "unstranded" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -p -s 0 -g gene_name -O --countReadPairs -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    fi
fi
if [ $p == "single-end" ]; then
    if [ $strand == "firststrand" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -s 2 -g gene_name -O -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    elif [ $strand == "secondstrand" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -s 1 -g gene_name -O -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    elif [ $strand == "unstranded" ]; then
        featureCounts -T $SLURM_CPUS_PER_TASK -s 0 -g gene_name -O -a $AnnotGTF -o ${CLUSTERDIR}/CountsTable_genename.txt $DIR/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*Aligned.sortedByCoord.out.bam
    fi
fi



