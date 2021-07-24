#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.reads = "$baseDir/data/raw_reads/*_{1,2}.fq.gz"
params.reference = "$baseDir/ref_genome/*.fna"
params.log = "fastqc_logs"
params.trim = "data/trimmed"
params.index = "ref_genome"
params.mapped = "data/mapped"
params.merged = "data/merged"

log.info """\
         SAMPLE MAPPING
         ===================================
         sample reads : ${params.reads}
         reference    : ${params.reference}
         output dir   : ${params.merged}
         """
         .stripIndent()


/*
 * Pipeline for mapping of a low-coverage sample
 * from raw reads to merged bam file.
 */

Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set{ read_pairs_ch }

process unzip {

    input:
        tuple sample_id, path(reads) from read_pairs_ch

    output:
        tuple sample_id, path('*') into unzip_ch, unzip2_ch
  
    script:
    """
    gzip -d --force $reads
    """
}

process fastqc {
    label 'toolbox'
    publishDir "fastqc_logs", mode:'copy'

    input:
        tuple sample_id, path(unzipped) from unzip_ch

    output:
        path '*' into fastqc_ch

    script:
    """
     	fastqc -f fastq $unzipped
    """
}

process trim {
    publishDir params.trim, mode:'copy'

    input:
        tuple name, path(unzipped) from unzip2_ch
    
    output:
        tuple name, path("${name}_{1,2}_trimmed.fq") into trimmed_ch
        
    script:
    """
    java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE $unzipped \
    ${name}_{1,2}_trimmed.fq ${name}_{1,2}_trimmed_unpaired.fq \
    SLIDINGWINDOW:4:15
    """
}

process index {
    label 'index'
    publishDir params.index, mode:'copy'
    
    input:
        path reference from params.reference
     
    output:
        file "*.{amb,ann,bwt,pac,sa}" into index_ch

    script:       
    """
    bwa index $reference 
    """
}

process bwa_mem {
    
    input:
        tuple name, path(reads) from trimmed_ch
        path(bwa_index) from index_ch
        path reference from params.reference

    output:
        path "${name}_bwa.bam" into sorted_ch

    script:
    """
    bwa mem -t 29 $reference $reads | samtools view -bS - > \
	${name}_bwa.bam
    """
}

process sam_sort {
    label 'samtools'

    input:
        path(reads) from sorted_ch

    output:
        path "${reads.simpleName}_bwa.bam" into sorted2_ch    
  
    script:    
    """
    samtools sort $reads -o  \
    ${reads.simpleName}_bwa.bam
    """
}

process sam_qual_filter {
    label 'samtools'

    input:
        path(reads) from sorted2_ch

    output:
        path "${reads.simpleName}_minq10_sorted.bam" into quality_filtered, quality_filtered2

    script:
    """
    samtools view -h -q 10 $reads | samtools view -buS - | samtools sort -o ${reads.simpleName}_minq10_sorted.bam   
    """

}

process sam_index {
    label 'samtools'

    input:
        path(reads) from quality_filtered

    output:
        path "*" into sam_indexed

    script:
    """
    samtools index $reads 
    """

}

process	merge { 
    label 'samtools'

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

process remove_duplicates {
    
    input:
        path input from sam2_merged.collectFile(name: 'merged.bam')

    output:
        path "merged_dedup.bam" into dedup_ch
    
    script:
    """
    samtools rmdup $input merged_dedup.bam
    """
}

process clip_overlap {
    label 'toolbox'

    publishDir params.merged, mode:'copy'

    input:
        path(input) from dedup_ch

    output:
        path "merged_dedup_overlapclipped.bam" into clipoverlap_ch

    script:
    """
    bam clipOverlap --in $input --out merged_dedup_overlapclipped.bam
    """
}

process clip_overlap_idx {
    label 'toolbox'
    
    publishDir params.merged, mode:'copy'

    input:
        path(input) from clipoverlap_ch

    output:
        path "merged_dedup_overlapclipped.bam.bai" into clipoverlap_idx

    script:
    """
    samtools index $input
    """
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}