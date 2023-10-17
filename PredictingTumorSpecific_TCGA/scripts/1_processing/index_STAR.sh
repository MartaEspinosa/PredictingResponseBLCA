#!/bin/bash
#SBATCH -p short,normal,long		# Partition to submit to
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu 14Gb	# Memory in MB
#SBATCH -J Idx_Gencode_v38_hg38_readlength75		# job name
#SBATCH -o /users/genomics/marta/logs/Idx_Gencode_v38_hg38_readlength75.%j.out	# File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/Idx_Gencode_v38_hg38_readlength75.%j.err	# File to which standard err will be written

####Change output and output error paths to redirect the log files generated

DIR=$1
GENOMEDIR=$2

ANNOTGENE=$GENOMEDIR/Annot_files_GTF
FASTADIR=$GENOMEDIR/GRCh38
GNMIDX=$CLUSTERDIR/Index_Genomes_STAR_2.8/
mkdir $GNMIDX
GNMIDX=$GNMIDX/Idx_Gencode_v38_hg38_readlength75
mkdir $GNMIDX

######################################################################################################
#####################################ALIGNMENT########################################################
module load STAR/2.7.8a-GCC-10.2.0

# Generate index with
#sjdbOverhang (def 100) length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)

STAR --runThreadN $SLURM_CPUS_PER_TASK --runMode genomeGenerate --genomeDir $GNMIDX --genomeFastaFiles $FASTADIR/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $ANNOTGENE/gencode.v38.primary_assembly.annotation.gtf --sjdbOverhang 74

#scp -rp $GNMDIX/* marta@hydra:$GENOMEDIR/Index_Genomes_STAR/Idx_Gencode_v38_hg38_readlength75
#rm -r $GNMDIX/*

