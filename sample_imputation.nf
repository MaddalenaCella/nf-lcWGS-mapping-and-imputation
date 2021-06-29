#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.bamfile = "$baseDir/data/merged/merged_dedup_overlapclipped.bam" //bam file generated from mapping nextflow script
params.indexed = "$baseDir/data/merged/merged_dedup_overlapclipped.bam.bai" //bai file generated from mapping nextflow script
params.snps = "$baseDir/data/snps/chr*.txt" //this is a list of SNPs found in each of the chromosomes in the reference panel
params.phased = "$baseDir/data/chr_phased_ref/chr*.out.gt.vcf.gz" //phased chromosomes of the reference panel
params.phenopos = "$baseDir/data/phenopos/angsd_final.file" //list of mendelian and QTLs of interest
params.reference = "$baseDir/ref_genome/EquCab3.fna"
params.samplename = "$baseDir/data/sample_name/name.txt" //name you want to give to sample vcf
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
    .into{ snp_ch; snp2_ch }

Channel
    .fromPath( params.bamfile, checkIfExists: true )
    .into{ bam_ch; bam2_ch }

Channel
    .fromPath( params.phased, checkIfExists: true )
    .into{ phased_ch; phased2_ch; phased3_ch }

Channel
    .fromPath( params.phenopos, checkIfExists: true )
    .into{ phenopos_ch; phenopos2_ch }


process site_index {
    label 'angsd'

    input:
        path sites_files from snp_ch
    output:    
        path "*.{bin,idx}" into snp_idx_ch
    script:
    """
    /opt/view/bin/angsd sites index $sites_files
    """
}

//before SNP calling I need to split global bam file into chromosomes with samtools and then do SNP calling on those files
//https://genome.sph.umich.edu/wiki/BamUtil:_splitChromosome this could be an alternative if bamtools does not work
process split_bam {
    label 'bamtools'

    input:
        path bam from bam_ch

    output:    
        path ('*.REF_chr{1,2}.bam') into chr_bam_ch, chr_bam2_ch

    script:
    """
    bamtools split -in $bam -reference
    """
}

process clip_overlap_idx {
    label 'toolbox'

   //publishDir params.merged, mode:'copy'

    input:
	path input from chr_bam2_ch.flatten()

    output:
	path "*.bai" into clipoverlap_idx_ch

    script:
    """
    samtools index $input
    """
}

process SNP_calling {
    label 'angsd'

    input:
	path snps from snp2_ch
        path snp_idx from snp_idx_ch
        path bam from chr_bam_ch.flatten()
        path bam_idx from clipoverlap_idx_ch
        path index from params.indexed
        path ref from params.reference

    output:
	path "${snps.simpleName}.vcf.*" into SNP_ch, SNP2_ch
    script:
    """
    /opt/view/bin/angsd -i $bam -ref $ref -P 5 \
    -out ${snps.simpleName} -sites $snps -r ${snps.simpleName}: \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
    -checkBamHeaders 0 \
    -dovcf 1 -GL 1 -doGlf 2 -doPost 1 -doMajorMinor 1 -SNP_pval 1e-3 -doMaf 1 \
    -doGeno 8 -dumpCounts 2 -doDepth 1 -doCounts 1
    """
}

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
//do I need an index before???
//need to add this before imputation!!!
process reheader {

    label 'bcftools'

    input:
        path name from params.samplename
        path calls from SNP2_ch
    
    output:
        path ('*') into reheaded_ch
        
    script:
    """
    bcftools reheader -s $name $calls \
    -o ${calls.simpleName}.reheaded.vcf.gz
    """
}

process imputation {

    publishDir params.imputed, mode:'copy'

    input:
        path snps from reheaded_ch
        path deco_ref from decompressed_ch

    output:
        path ('*') into imputed_ch, imputed2_ch

    script:
    """
    beagle \
    ref=$deco_ref \
    gt=$snps \
    out=${snps.simpleName}.out window=100000000
    """
}

process index_imputed {

    label 'tabix'

    input:
        path snps_beagle from imputed_ch
    
    output:
        path ('*') into idx_imputed_ch
        
    script:
    """
    tabix -f -p vcf $snps_beagle
    """
}

process concat_imputed {

    label 'bcftools'

    input:
        tuple sample_id, path(snps_beagle) from imputed2_ch
        path idx from idx_imputed_ch

    output:
        path ('Merged.vcf.gz') into concat_ch, concat2_ch, concat3_ch
        
    script:
    """
    bcftools concat $snps_beagle -O z -o Merged.vcf.gz
    """

}

process index_concat {

    label 'tabix'

    input:
        path concatenated from concat_ch
    
    output:
        path ('*') into idx_concatenated_ch
        
    script:
    """
    tabix -f -p vcf $concatenated
    """
}


process annotate {

    label 'bcftools'
    publishDir params.imputed, mode:'copy'

    input:
        path reheaded from reheaded_ch
       
    output:
        path ('*') into annotated_ch, annotated2_ch
        
    script:
    """
    bcftools annotate -x 'FORMAT','INFO' $reheaded \
    -o MergedNoFI.vcf.gz

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

process index_phased_ref {

    label 'tabix'

    input:
        path phased from phased2_ch
    
    output:
        path ('*') into idx_phased_ch
        
    script:
    """
    tabix -f -p vcf $phased
    """
}

process concat_phased {

    label 'bcftools'

    input:
        tuple sample_id, path(phased) from phased3_ch
        path idx from idx_phased_ch

    output:
        path ('RefPanel.vcf.gz') into concat_ref_ch, concat2_ref_ch, concat3_ref_ch
        
    script:
    """
    bcftools concat $phased -O z -o RefPanel.vcf.gz
    """

}

process index_concat_ref {

    label 'tabix'

    input:
        path concatenated from concat_ref_ch
    
    output:
        path ('*') into idx_concatenated_ref_ch
        
    script:
    """
    tabix -f -p vcf $concatenated
    """
}

process annotate_ref {

    label 'bcftools'
    publishDir params.ref, mode:'copy'

    input:
        path concatenated from concat2_ref_ch
       
    output:
        path ('*') into annotated_ref_ch, annotated2_ref_ch
        
    script:
    """
    bcftools annotate -x 'FORMAT','INFO' $concatenated \
    -o RefPanelNoFI.vcf.gz

    """
}

process index_annotated_ref {

    label 'tabix'

    input:
        path annotated from annotated_ref_ch
    
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
        path annotated_ref from annotated2_ref_ch
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
        path ('final_vcf*') into nosingletons_ch
        
    script:
    
    """
    vcftools --vcf $combined \
    --exclude-positions $singletons \
    --recode --recode-INFO-all --out final_vcf
    """
}

process plink_convert_admixture {

    publishDir params.ancestry, mode:'copy'

    input:
        path no_singletons from nosingletons_ch
        
    output:
        path ('admixture12*') into admixture_ch
        
    script:
    
    """
    plink1.9 --vcf $no_singleton --recode12 --horse --out admixture12
    """    

}

process pheno_search {

    label 'vcftools'
    publishDir params.phenotype, mode:'copy'

    input:
        path concatenated from concat3_ch
        path pheno_pos from phenopos_ch
        
    output:
        path ('pheno_sample*') into phenotype_ch
        
    script:
    
    """
    vcftools --gzvcf $concatenated --recode --out pheno_sample --positions $pheno_pos
    """    

}

process pheno_search_ref {

    label 'vcftools'
    publishDir params.phenotype, mode:'copy'

    input:
        path concatenated from concat3_ref_ch
        path pheno_pos from phenopos2_ch
        
    output:
        path ('pheno_pos_Ref*') into phenotype_ch
        
    script:
    
    """
    vcftools --gzvcf $concatenated --recode --out pheno_pos_Ref --positions $pheno_pos
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

