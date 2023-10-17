# Tumorspecific_neoags

Conventional neoantigens are considered those that come from mutated peptides or from tumor-specific proteins. The second class has relied on cancer/testis antigens. However, we want to focus not only on those but also on lncRNA and novel transcripts that are specifically expressed in the tumor samples, since a great proportion of them contain potential translatable open reading frames.

The approach goes from the raw fastq files to the generation of a list of peptides with binding affinity for the major histocompatibility complex.

Both Python3 and R (v4.1.2) languages are used, as well as bash.

## Requirements

`Requirements.txt` file contains all programas required to run the pipeline.

## Reference genome

GRCh38 human genome version is used. Other version should work too, as long as paths are modified.

## Pipeline Summary

### 0. Preprocessing

Fastq files are obtained either from GEO datasets or from TCGA. On the second option, bams are downloaded and reverted to fastq files.

It is essential the last step of this first part, which creates a comma-separated file containing the identifiers of the patients and corresponding tumor and normal samples.

### 1. Processing

Some quality checks are done (FastQScreen, FastQC) and depending on the results, cutadapt can be used to remove the adaptors. However, since later the alignment is done using STAR software, it takes into account the possibility of having adaptors, so not removing them would not be a big deal in this case.
STAR is run to only keep uniquely mapped reads and it uses the two-pass mode.

Having stranded datasets is crucial for the assembly of novel transcripts, so RSeQC is used to know whether strand information is provided or not and, if so, in which orientation.

### 2. Tumor-Expressed

With the intention of simplyfing potential complications that may arise because of isoforms, we only keep the longest transcript per gene, having a relation 1:1 between transcripts and genes.

FeatureCounts quantifies and counts are converted to TPM.

We consider a minimum expression of 1 TPM (minimum expression per novel transcripts may need to be adapted, depending on the case). These transcripts will constitute the **tumor-expres subset**.

### 3. Tumor-Specific

Transcripts expressed in tumor samples but also in normal samples would not be relevant for neoantigens identification, since the body will be able to recognize those peptides as self. So, here we use the expression data from the normal samples to only keep those transcripts that are expressed (> 3 TPM) in the tumor sample but not in the normal one (< 0.1 TPM). This is done per patient, obtaining in the end a non-redundant list of all transcripts that are tumor-specific in at least one patient.

Another issue is whether those that fulfill tumor-specific expression conditions are expressed in healthy tissues or not. If so, we are not interested in them any more. For this, we only allow a median of < 0.5 TPM according to GTEx project with exception of testis and ovaries.

### 4. Tumor-Specific Criteria

Once the definition of the main subset of tumor-specific genes is done, we can change the expression cut-off to consider a smaller subset of genes.
Moreover, we want to increase tumor-specificity by removing genes that are tumor-specific in some patients but they are expressed more than 1 TPM in other patients.
Another applied criteria is to select only those genes tumor-specific in > 10% of the analyzed peptides.
