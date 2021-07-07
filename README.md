# nf-lcWGS-mapping-and-imputation

[![Nextflow version: >=0.30.2](https://img.shields.io/badge/nextflow-%E2%89%A50.30.2-brightgreen.svg)](http://nextflow.io)
[![Singularity version: >=3.x](https://img.shields.io/badge/singularity-3.x-blue.svg)](http://sylabs.io/singularity)

This pipeline allows mapping and imputation of non-human samples sequenced at low-coverage.

It uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operationg system (Linux, MacOS, etc)
- Java 8
- [Singularity](http://singularity.lbl.gov) engine

## Quickstart

1. Install Nextflow by using the following command:

    ```
    curl -s https://get.nextflow.io | bash
    ```

2. Download singularity images on /singularity directory:

    ```
    sh set-up-images.sh

    ```
3. Add the input parameters in the /data directory. For more details look at the dedicated section below. 

**NOTE**: the very first time you execute it, it will take a few minutes to download the pipeline from this GitHub repository and the associated Singularity images needed to execute the pipeline.

### Using Singularity

Singularity is the preferred container engine for running the pipeline in an HPC environment. In order to minimize the amount of issues that could arise we recommend the use of Singularity version 3.0 or higher. In my case I used Singularity version 3.6.1. The singularity images used were obtained from docker hub and singularity hub.

## Pipeline input
Before you deploy the mapping pipeline, make sure you have the raw reds file in fq.gz format in the data/raw_reads/ directory and the reference genome sequence of your species in the ref_genome/ directory. Edit the name of your files as appropriate in the sample_map.nf script. 

Before you deploy the imputation pipeline, the input files you need to have in your code directory are the following:
* the mapped bam (.bam) file and its index file (.bai) in the data/merged/ directory (this is the output of the mapping pipeline, but you can add your own files and change the name accordingly in the sample_imputation.nf script)
* the reference genome sequence of your species in the ref_genome/ directory.
* a list of SNP calls in each chromosome of the reference panel (one file per chromosome). This can be obtained with the following line of code and the bioinformatic tool vcftools given you have all the reference's phased vcf.gz for each chromosome in your current directory. 

```
vcftools --gzvcf chr$NUM.out.gt.vcf.gz --out chr$NUM.filtered --site-quality

##generate SNP sites list in a format compatible with ANGSD
awk 'NF{NF-=1};1' <chr$NUM.filtered.lqual >chr$NUM.sites
sed '1d' chr$NUM.sites > chr$NUM.txt

```
* the phased genotype calls of the reference panel divided by chromosome in the data/chr_phased_ref/ directory;
* a list of mendelian and QTL traits of interest in the format: chr start end if you are interested in finding whether your sample's genotype for some traits of interest. If you are not interested in this analysis, please comment out the process named "pheno_search" before running the pipeline;
* the name you want to use to identify your sample in a text file in the data/samplename directory.

## Pipeline results
### Output files

The mapping pipeline produces two outputs: 
* the mapped reads in the data/mapped directory
* a merged and processed file (deduplicated and overlapping read pairs clipped) in the data/merged directory. If you want to intermediate files you can add the "publishDir" command as I show in the following example:

```
process	merge { 
    label 'samtools'
    publishDir params.merged, mode:"copy" //this line publishes the output of the "merge" process in the data/merge directory

    input:
        path(input) from quality_filtered2.collect()
        path(index) from sam_indexed

    output:
        path "merged.bam" into sam2_merged     

    script:
    """
    samtools merge -f merged.bam $input                                                                                      
    """
}
```

The imputation pipeline produces:
* a combined vcf file (with variants from both the reference panel and the horse sample) in the results/combined directory;
* this file converted in a Plink and Admixture compatible format: .bed, .bim and .fam files in the results/ancestry directory; 
* a vcf file with the genotypes for the SNPs of interest found in the horse sample. 

## Pipeline configuration

### Executors

Nextflow provides different `Executors` run the processes on the local machine, on a computational cluster or different cloud providers without the need to change the pipeline code.

By default the local executor is used, but it can be changed by using  the [executor](https://www.nextflow.io/docs/latest/config.html#scope-executor) configuration scope.

For example, to run the pipeline in a computational cluster using Sun Grid Engine you can create a `nextflow.config` file in your current working directory with something like:

```
process {
    executor = 'sge'
    queue    = 'my-queue'
    penv     = 'smp'
}
```


## Run the pipeline




##  Tools versions

The versions of the tools that have been tested with the pipeline are the following:

- [vcftools v0.1.16](https://github.com/vcftools/vcftools/releases/tag/v0.1.16)
- [samtools v1.2.1](https://github.com/samtools/samtools/releases/tag/1.2.1)
- [bamtools v2.5.1](https://github.com/pezmaster31/bamtools/releases/tag/v2.5.1)
- [bcftools v1.9](https://github.com/samtools/bcftools/releases/tag/1.9)
- [htslib v1.12](https://github.com/samtools/htslib/releases/tag/1.12)
- [plink v1.9](https://github.com/singemanator/MGL804-PLINK1.9)
- [tabix v1.9](https://github.com/samtools/tabix)
- [vcflib v1.0.2](https://github.com/vcflib/vcflib/tree/v1.0.2)
- [picard v1.139](https://github.com/broadinstitute/picard/releases/tag/1.139)
- [bwa v0.7.12](https://github.com/lh3/bwa/releases/tag/0.7.12)
- [angsd v1.7](https://github.com/ANGSD/angsd)
- [beagle v4.1](https://faculty.washington.edu/browning/beagle/b4_1.html)
- [fastqc v.0.11.4](https://github.com/s-andrews/FastQC)
- [bamUtil v1.0.13](https://github.com/statgen/bamUtil/releases/tag/v1.0.13)