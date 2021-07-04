#code to run to dowload all the singularity images necessary to run the pipeline on the HPC
cd singularity

## tabix image
singularity pull docker://biocontainers/tabix:v1.9-11-deb_cv1

## angsd image
singularity pull docker://lifebitai/angsd:latest

## picard image
singularity pull docker://biocontainers/picard:v1.139_cv3

## bwa and samtools image
singularity pull docker://dukegcb/bwa-samtools:latest

## trimmomatic image
singularity pull shub://jlboat/BioinfoContainers:trimmomatic

## MBD toolbox for BamUtil, fastqc, samtools
singularity pull docker://mateongenaert/mbdtoolbox:latest

## vcflib image
singularity pull docker://shollizeck/vcflib:latest

## vcftools image
singularity pull docker://biocontainers/vcftools:v0.1.16-1-deb_cv1

## bcftools image
singularity pull docker://biocontainers/bcftools:v1.9-1-deb_cv1

## htslib image
singularity pull docker://clinicalgenomics/htslib:1.12

## beagle image
singularity pull docker://lifebitai/beagle4:latest

## plink image
singularity pull docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1