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

The mapping pipeline produces two outputs: the mapped reads in the data/mapped directory and a merged and processed file (deduplicated and overlapping read pairs clipped) in the data/merged directory. If you want to intermediate files you can add the "publishDir" command as I show in the following example:

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

### Pipeline profiles

The Grape pipeline can be run using different configuration profiles. The profiles essentially allow the user to run the analyses using
different tools and configurations. To specify a profile you can use the [`-profiles` Nextflow option](http://www.nextflow.io/docs/latest/config.html#config-profiles).

The following profiles are available at present:


profile | description
|-|-|
 `gemflux`  | uses `GEMtools` for mapping pipeline and `Flux Capacitor` for isoform expression quantification
 `starrsem` | uses `STAR` for mapping and bigwig and `RSEM` for isoform expression quantification
 `starflux` | uses `STAR` for mapping and `Flux Capacitor` for isoform expression quantification

The default profile is `starrsem`.

## Run the pipeline

Here is a simple example of how you can run the pipeline:

```
nextflow -bg run grape-nf --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG -resume > pipeline.log
```

By default the pipeline execution will stop as far as one of the processes fails. This behaviour can be changed using the [errorStrategy](http://www.nextflow.io/docs/latest/process.html#errorstrategy) process directive, which can also be specified on the command line. For example, to ignore errors and keep processing you can use:

`-process.errorStrategy=ignore`.

It is also possible to run a subset of the pipeline steps using the option ``--steps``. For example, the following command will only run ``mapping`` and ``quantification``:

```
nextflow -bg run grape-nf --steps mapping,quantification --index input-files.tsv --genome refs/hg38.AXYM.fa --annotation refs/gencode.v21.annotation.AXYM.gtf --rg-platform ILLUMINA --rg-center-name CRG > pipeline.log
```

##  Tools versions

The pipeline can be also run natively by installing the required software on the local system or by using [Environment Modules](http://www.nextflow.io/docs/latest/process.html?#module).

The versions of the tools that have been tested so far with the `standard` pipeline profile are the following:

- [bamstats v0.3.4](https://github.com/guigolab/bamstats/releases/tag/v0.3.4)
- [bedtools v2.19.1](https://github.com/arq5x/bedtools2/releases/tag/v2.19.1)
- [KentUtils v308](https://github.com/ucscGenomeBrowser/kent/releases/tag/v308_base)
- [RSEM v1.2.21](https://github.com/deweylab/RSEM/releases/tag/v1.2.21)
- [RSeQC v2.6.4](http://rseqc.sourceforge.net/)
- [sambamba v0.7.1](https://github.com/biod/sambamba/releases/tag/v0.7.1)
- [samtools v1.3](https://github.com/samtools/samtools/releases/tag/1.3.1)
- [STAR v2.4.0j](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0j)
