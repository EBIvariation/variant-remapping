#!/usr/bin/env nextflow


// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    Process a VCF file containing SNV or short indels associated with the old genome to remap the variant's coordinates to the new genome.

    Inputs:
            --vcffile      original vcf file containing variants to be remapped [required]
            --oldgenome    source genome file used to discover variants from the file provided with --vcffile [required]
            --newgenome    new genome which the original vcf file will be mapped against [required]
    Outputs:
            --outfile      final vcf containing variants mapped to new genome [required]
    Parameters:
            --flankingseq  The length of the flanking sequences that generates reads (default 50)
            --scorecutoff  Percentage of the flanking sequences that should be used as the alignment score cut-off threshold (default 0.6)
            --diffcutoff   Percentage of the flanking sequences that should be used as the AS-XS difference cut-off threshold (default 0.04)
    """
}

params.flankingseq = 50
params.scorecutoff = 0.6
params.diffcutoff =  0.04
params.outfile = null
params.oldgenome = null
params.newgenome = null
params.vcffile = null
params.help = null

// Show help message
if (params.help) exit 0, helpMessage()

// Test input files
if (!params.vcffile || !params.oldgenome || !params.outfile || !params.newgenome) { 
    if (!params.vcffile)    log.warn('Provide a input vcf file using --vcffile')
    if (!params.oldgenome)  log.warn('Provide a fasta file for the old genome --oldgenome') 
    if (!params.outfile)    log.warn('Provide a path to the output vcf file using --outfile') 
    if (!params.newgenome)  log.warn('Provide a fasta file for the new genome file using --newgenome') 
    exit 1, helpMessage()
}
 
// Index files for both old and new genomes 
oldgenome_dir = file(params.oldgenome).getParent()
oldgenome_basename = file(params.oldgenome).getName()
oldgenome_fai = file("${params.oldgenome}.fai")
newgenome_dir = file(params.newgenome).getParent()
newgenome_basename = file(params.newgenome).getName()
newgenome_fai = file("${params.newgenome}.fai")

// basename and basedir of the output file to know how to name the output file
outfile_basename = file(params.outfile).getName()
outfile_dir = file(params.outfile).getParent()

/*
 * Uncompress VCF file
 */
process uncompressInputVCF {

    input:
        path "source.vcf"

    output:
        path "uncompressed.vcf", emit: vcf_file

    script:
        if ( file(params.vcffile).getExtension() == 'gz' )
            """
            gunzip -c source.vcf > uncompressed.vcf
            """
        else
            """
            ln -nfs source.vcf uncompressed.vcf
            """
}

/*
 * Store the original VCF header for later use
 */
process storeVCFHeader {

    input:
        path "source.vcf"

    output:
        path "vcf_header.txt", emit: vcf_header

    """
    bcftools view --header-only source.vcf | grep -v '^##FORMAT' >  vcf_header.txt
    """
}

include { prepare_old_genome; prepare_new_genome; prepare_new_genome_bowtie } from './prepare_genome.nf'
include { process_split_reads; process_split_reads_with_bowtie } from './variant_to_realignment.nf'


/*
 * Create the header for the output VCF
 */
process buildHeader {

    input:
        path "variants_remapped_sorted.vcf"
        path "vcf_header.txt"

    output:
        path "final_header.txt", emit: final_header

    """
    # Create list of contigs/chromosomes to be added to the header
    cut -f 1 variants_remapped_sorted.vcf | sort -u > contig_names.txt
    while read CHR; do echo "##contig=<ID=\${CHR}>"; done < contig_names.txt > contigs.txt
    # Add the reference assembly
    echo "##reference=${params.newgenome}" >> contigs.txt

    # Copy everything that isn't #CHROM (the column titles), ##contig or ##reference from the old header to a temp
    awk '(\$1 !~ /^##contig/ && \$1 !~ /^##reference/ && \$1 !~ /^#CHROM/) {print \$0}' vcf_header.txt > temp_header.txt

    # Add the two headers together and add the column names
    cat temp_header.txt contigs.txt > final_header.txt
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> final_header.txt
    """
}

/*
 * Add header to output VCF and merge with reference allele file
 */
process mergeHeaderAndContent {

    input:
        path "final_header.txt"
        path "variants_remapped_sorted.vcf"

    output:
        path "vcf_out_with_header.vcf", emit: final_vcf_with_header

    """
    # Add header to the vcf file:
    cat final_header.txt variants_remapped_sorted.vcf > vcf_out_with_header.vcf
    """
}

/*
 * Sort VCF file
 */
process sortVCF {

    input:
        path "variants_remapped.vcf"

    output:
        path "variants_remapped_sorted.vcf.gz", emit: variants_remapped_sorted_gz

    """
    bgzip variants_remapped.vcf
    tabix variants_remapped.vcf.gz
    bcftools sort -o variants_remapped_sorted.vcf.gz -Oz variants_remapped.vcf.gz
    """
}

/*
 * Run bcftools norm to swap the REF and ALT alleles if the REF doesn't match the new assembly
 */
process normalise {

    publishDir outfile_dir,
        overwrite: true,
        mode: "copy"

    input:
        path "variants_remapped_sorted.vcf.gz"
        path "genome.fa"

    output:
        path "${outfile_basename}", emit: final_output_vcf

    """
    bcftools norm -c ws -f genome.fa -N variants_remapped_sorted.vcf.gz -o ${outfile_basename} -O v
    """
}

/*
 * Create file containing remapping stats
 */
process calculateStats {

    publishDir outfile_dir,
        overwrite: true,
        mode: "copy"

    input:
        path outfile_basename

    output:
        path "${outfile_basename}.stats"

    """
    bcftools stats ${outfile_basename} > ${outfile_basename}.stats
    """
}

// Take variants remapped to the new genome and merge them back with the original header and sort the output
workflow finalise {
    take:
        variants_remapped
        vcf_header

    main:
        buildHeader(variants_remapped, vcf_header)
        mergeHeaderAndContent(buildHeader.out.final_header, variants_remapped)
        sortVCF(mergeHeaderAndContent.out.final_vcf_with_header)
        normalise(sortVCF.out.variants_remapped_sorted_gz, params.newgenome)
        calculateStats(normalise.out.final_output_vcf)
}


//process_with_minimap
// Workflow without a name is the default workflow that gets executed when the file is run through nextflow
workflow {
    main:
        prepare_old_genome(params.oldgenome)
        prepare_new_genome(params.newgenome)
        uncompressInputVCF(params.vcffile)
        storeVCFHeader(uncompressInputVCF.out.vcf_file)
        process_split_reads(
            uncompressInputVCF.out.vcf_file,
            params.oldgenome,
            prepare_old_genome.out.genome_fai,
            prepare_old_genome.out.genome_chrom_sizes,
            params.newgenome,
            prepare_new_genome.out.genome_fai
        )
        finalise(process_split_reads.out.variants_remapped, storeVCFHeader.out.vcf_header)
}

//process_with_bowtie
workflow process_with_bowtie {
    main:
        prepare_old_genome(params.oldgenome)
        prepare_new_genome_bowtie(params.newgenome)
        uncompressInputVCF(params.vcffile)
        storeVCFHeader(uncompressInputVCF.out.vcf_file)
        process_split_reads_with_bowtie(
            uncompressInputVCF.out.vcf_file,
            params.oldgenome,
            prepare_old_genome.out.genome_fai,
            prepare_old_genome.out.genome_chrom_sizes,
            params.newgenome,
            prepare_new_genome_bowtie.out.genome_fai,
            prepare_new_genome_bowtie.out.bowtie_indexes
        )
        finalise(process_split_reads_with_bowtie.out.variants_remapped, storeVCFHeader.out.vcf_header)
}
