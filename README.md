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


**NOTE**: Fastq files from paired-end data will be grouped together by `runID`.

Here is an example from the test run:

```
sample1  test1   data/test1_1.fastq.gz   fastq   FqRd1
sample1  test1   data/test1_2.fastq.gz   fastq   FqRd2
```

Sample and id can be the same in case you don't have/know sample identifiers:

```
run1  run1   data/test1_1.fastq.gz   fastq   FqRd1
run1  run1   data/test1_2.fastq.gz   fastq   FqRd2
```

## Pipeline results

The paths of the resulting output files and the corresponding metadata are stored into the `pipeline.db` file (`TSV` formatted) which sits inside the current working folder. The format of this file is the same as the index file with few more columns:

||||
|-|-|-|
1 | `sampleID` | the sample identifier, used to merge bam files in case multiple runs for the same sample are present
2 | `runID`          | the run identifier (e.g. `test1`)
3 | `path`        | the path to the fastq file
4 | `type`        | the type (e.g. `bam`)
5 | `view`        | an attribute that specifies the content of the file (e.g. `GenomeAlignments`)
6 | `readType`    | the input data type (either `Single-End` or `Paired-End`)
7 | `readStrand`  | the inferred experiment strandedness if any (it can be `NONE` for unstranded data, `SENSE` or `ANTISENSE` for single-end data, `MATE1_SENSE` or `MATE2_SENSE` for paired-end data.)

Here is an example from the test run:

```
sample1   test1   /path/to/results/sample1.contigs.bed    bed      Contigs                     Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.isoforms.gtf   gtf      TranscriptQuantifications   Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.plusRaw.bw     bigWig   PlusRawSignal               Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.genes.gff      gtf      GeneQuantifications         Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/test1_m4_n10.bam       bam      GenomeAlignments            Paired-End   MATE2_SENSE
sample1   test1   /path/to/results/sample1.minusRaw.bw    bigWig   MinusRawSignal              Paired-End   MATE2_SENSE
```

### Output files

The pipeline produces several output files during the workflow execution. Many files are to be considered temporary and can be removed once the pipeline completes. The following files are the ones reported in the `pipeline.db` file and are to be considered as the pipeline final output.

#### Alignments to the reference genome

|views|
|-|
|`GenomeAlignments`|

This BAM file contains information on the alignments to the reference genome. It includes all the reads from the FASTQ input. Reads that do not align to the reference are set as unmapped in the bam file. The file can be the product of several steps of the pipeline depending on the given input parameters. It is initially produced by the `mapping` step, then it can be the result of merging of different runs from the same experiment and finally it can run through a marking duplicates process that can eventually remove reads that are marked as duplicates.

#### Alignments to the reference transcriptome

|views|
|-|
|`TranscriptomeAlignments`|

This BAM file contains information on the alignments to the reference transcriptome. It is generally used only for expression abundance estimation, as input in the `quantification` process. The file is generally produced in the `mapping` process and can be the result of merging of different runs from the same experiment.

#### Alignments statistics

|views|
|-|
|`BamStats`|

This JSON file contains alignment statistics computed with the [bamstats](https://github.com/guigolab/bamstats) program. It also reports RNA-Seq quality check metrics agreed within the IHEC consortium.

#### Signal tracks

|views|
|-|
|`RawSignal`|
|`MultipleRawSignal`|
|`MinusRawSignal`|
|`PlusRawSignal`|
|`MultipleMinusRawSignal`|
|`MultiplePlusRawSignal`|

These BigWig files (one or two, depending on the strandedness of the input data) represent the RNA-Seq signal.

#### Contigs

|views|
|-|
|`Contigs`|

This BED file reports RNA-seq contigs computed from the pooled signal tracks.

#### Quantifications

|views|
|-|
|`GeneQuantifications`|
|`TranscriptQuantifications`|

These two files report abundances for genes and transcripts in the processed RNA-seq samples. The format can be either GFF or TSV depending on the tool used to perform the quantification.

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
