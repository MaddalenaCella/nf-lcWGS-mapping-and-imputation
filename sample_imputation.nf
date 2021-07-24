#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.bamfile = "$baseDir/data/merged/merged_dedup_overlapclipped.bam" //bam file generated from mapping nextflow script
params.indexed = "$baseDir/data/merged/merged_dedup_overlapclipped.bam.bai" //bai file generated from mapping nextflow script
params.snps = "$baseDir/data/snps/chr*.txt" //this is a list of SNPs found in each of the chromosomes in the reference panel
params.phased = "$baseDir/data/chr_phased_ref/chr*.out.gt.vcf.gz" //phased chromosomes of the reference panel
params.phenopos = "$baseDir/data/phenopos/angsd_final.file" //list of mendelian and QTLs of interest
params.reference = "$baseDir/ref_genome/*.fna"
params.samplename = "$baseDir/data/samplename/name.txt" //name you want to give to sample vcf
params.ref = "$baseDir/results/phased_ref"
params.imputed = "$baseDir/results/imputed"
params.combined = "$baseDir/results/combined"
params.ancestry = "$baseDir/results/ancestry"
params.phenotype = "$baseDir/results/phenotype"


log.info """\
         SAMPLE IMPUTATION
         ===================================
         sample bam   : ${params.bamfile}
         reference    : ${params.reference}
         phased panel : ${params.phased}
         output plink : ${params.ancestry}
         """
         .stripIndent()

/*
 * Pipeline for imputing variants in a low-coverage sample 
 * and coverting it to plink format for further analysis
 * starting from a mapped file in bam format and phased vcf
 * files for each chromosome of the reference panel of individuals
 */

Channel
    .fromPath( params.snps, checkIfExists: true )
    .set{ snp_ch }

Channel
    .fromPath( params.phased, checkIfExists: true )
    .into{ phased_ch; phased2_ch }

Channel
    .fromPath( params.phenopos, checkIfExists: true )
    .into{ phenopos_ch; phenopos2_ch }

process decompress {
    
    input:
        path phased_gz from phased_ch
    output:
        path ('*') into decompressed_ch
    script:
    """
    bgzip -d $phased_gz
    """
}

process SNP_calling {
    label 'angsd'

    input:
	path snps from snp_ch //you can remove if it doesn't work
        path bam from params.bamfile
        path bam_idx from params.indexed
        path ref from params.reference

    output:
	path "${snps.simpleName}.bcf" into SNP_call_ch

    script:
    """
    angsd sites index $snps

    angsd -i $bam -ref $ref -P 5 \
    -out ${snps.simpleName} -sites $snps -r ${snps.simpleName}: \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -checkBamHeaders 0 \
    -dobcf 1 -GL 1 -doGlf 2 -doPost 1 -doMajorMinor 1 -SNP_pval 1e-1 -doMaf 1 \
    -doGeno 1 -dumpCounts 2 -doDepth 1 -doCounts 1 --ignore-RG 0
    """
}

process bcf_to_vcf {
    label 'bcftools'
    
    input:
        path bcf from SNP_call_ch.flatten()
        path name from params.samplename
    
    output:
        path "${bcf.simpleName}.reheaded.vcf.gz" into VCF_ch
    
    script:
    """
    bcftools view -O z -o ${bcf.simpleName}.vcf.gz $bcf

    bcftools reheader -s $name ${bcf.simpleName}.vcf.gz \
    -o ${bcf.simpleName}.reheaded.vcf.gz
    """
}

process imputation {

    input:
        path snps from VCF_ch
        path deco_ref from decompressed_ch

    output:
        path "${snps.simpleName}*.gz" into imputed_ch

    script:
    """
    beagle \
    ref=$deco_ref \
    gt=$snps \
    out=${snps.simpleName}.out window=100000000
    """
}

process index_and_concat {

    label 'bcftools'
    publishDir params.imputed, mode:'copy'

    input:
        path snps_beagle from imputed_ch
    
    output:
        path ('Merged.vcf.gz') into concat_ch, concat2_ch, concat3_ch
        
    script:
    """
    bcftools index -f $snps_beagle
    bcftools concat $snps_beagle -O z -o Merged.vcf.gz
    """
}

process index_concatenated {

    label 'tabix'
    publishDir params.imputed, mode:'copy'

    input:
        path concatenated from concat_ch.collectFile(name: 'Merged.vcf.gz')
    
    output:
        path ('Merged.vcf.gz.tbi') into idx_concat_ch, idx_concatenated_ch
        
    script:
    """
    tabix -f $concatenated
    """
}

process annotate {

    label 'bcftools'
    publishDir params.imputed, mode:'copy'

    input:
        path concatenated from concat2_ch.unique()
        path idx from idx_concat_ch

    output:
        path ('MergedNoFI.vcf.gz') into annotated_ch, annotated2_ch
        
    script:
    """
    bcftools annotate -x 'FORMAT','INFO' $concatenated \
    -O z -o MergedNoFI.vcf.gz
    """
}

process index_annotated {
    
    label 'tabix'

    input:
        path annotated from annotated_ch
    
    output:
        path ('*') into idx_annotated_ch
        
    script:
    """
    tabix -f -p vcf $annotated
    """
}


process concat_phased {

    label 'bcftools'

    input:
        path phased from phased2_ch

    output:
        path ('RefPanel.vcf.gz') into concat_ref_ch
        
    script:
    """
    bcftools index -f $phased
    bcftools concat $phased -O z -o RefPanel.vcf.gz
    """

}

process annotate_ref {

    label 'bcftools'
    publishDir params.ref, mode:'copy'

    input:
        path concatenated from concat_ref_ch.unique()
       
    output:
        path ('*') into annotated_ref_ch, annotated2_ref_ch
        
    script:
    """
    bcftools annotate -x 'FORMAT','INFO' $concatenated \
    -O z -o RefPanelNoFI.vcf.gz

    """
}

process index_annotated_ref {

    label 'tabix'

    input:
        path annotated from annotated_ref_ch.collectFile(name: 'RefPanelNoFI.vcf.gz')
    
    output:
        path ('*') into idx_annotated_ref_ch
        
    script:
    """
    tabix -f -p vcf $annotated
    """
}

process combine {

    publishDir params.combined, mode:'copy'

    input:
        path annotated from annotated2_ch
        path annotated_ref from annotated2_ref_ch.collectFile(name: 'RefPanelNoFI.vcf.gz')
        path idx_ann from idx_annotated_ch
        path idx_ann_ref from idx_annotated_ref_ch
    
    output:
        path ('combinedNoFI.vcf') into combined_ch, combined2_ch
        
    script:
    """
    vcfcombine $annotated_ref $annotated > combinedNoFI.vcf
    """
}

process singletons_list {

    label 'vcftools'

    input:
        path combined from combined_ch
        
    output:
        path ('combined*') into singletons_ch
        
    script:
    
    """
    vcftools --vcf $combined --singletons --out combined
    """
}

process remove_singletons {

    label 'vcftools'
    publishDir params.ancestry, mode:'copy'

    input:
        path combined from combined2_ch 
        path singletons from singletons_ch
        
    output:
        path ('final_vcf*') into nosingletons_ch, nosingletons2_ch
        
    script:
    
    """
    vcftools --vcf $combined \
    --exclude-positions $singletons \
    --recode --recode-INFO-all --out final_vcf
    """
}

process LD_pruning_positions {

    label 'plink'

    input:
        path nosingletons from nosingletons_ch

    output:
    	path ('final-noLD.prune.in') into to_keep_ch
	    path ('final-noLD.prune.out') into not_to_keep_ch

    script:

    """
    plink1.9 --vcf $nosingletons --const-fid 0 --maf 0.01 --indep-pairwise 50 5 0.2 --horse --out final-noLD
    """
}

process LD_pruning { //plink is compatible with both Admixture and Plink --pca. Two commonly used tools 

    label 'plink'
    publishDir params.ancestry, mode: 'copy'

    input:
	    path prune_in from to_keep_ch
	    path nosingletons from nosingletons2_ch

    output:
       path ('PCA_ready.bed') into bed_ch
       path ('PCA_ready.bim') into bim_ch
       path ('PCA_ready.fam') into fam_ch
       path ('PCA_ready.nosex') into nosex_ch
       path ('PCA_ready.log') into log_ch
       
    script:

    """
    plink1.9 --vcf $nosingletons --extract $prune_in --make-bed --horse --const-fid 0 --out PCA_ready
    """
}

process pheno_sites_idx {
    label 'angsd'
    
    input:
        path phenofile from phenopos2_ch

    output:
        path ("*") into idx_pheno_file_ch

    script:
    """
    angsd sites index $phenofile
    """
}

process pheno_search {

    label 'vcftools'
    publishDir params.phenotype, mode:'copy'

    input:
        path concatenated from concat3_ch
        path pheno_pos from phenopos_ch
	path idx_pos from idx_pheno_file_ch	
	path idx from idx_concatenated_ch
        
    output:
        path ('pheno_sample*') into phenotype_ch
        
    script:
    
    """
    vcftools --gzvcf $concatenated --recode --out pheno_sample --positions $pheno_pos
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
