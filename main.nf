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
    """
}

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

// basename and basedir of the output file to know how to name the output files
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
    bcftools view --header-only source.vcf > vcf_header.txt
    """
}

include { prepare_old_genome; prepare_new_genome; prepare_new_genome_bowtie } from './prepare_genome.nf'
include { process_split_reads; process_split_reads_mid; process_split_reads_long; process_split_reads_with_bowtie } from './variant_to_realignment.nf'


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

    # Add variant remapping INFO definition to the header
    echo -e '##INFO=<ID=st,Number=1,Type=String,Description="Strand change observed in the alignment.">' >> temp_header.txt
    echo -e '##INFO=<ID=rac,Number=1,Type=String,Description="Reference allele change during the alignment.">' >> temp_header.txt
    echo -e '##INFO=<ID=nra,Number=0,Type=Flag,Description="Novel reference allele that was not observed in the previous set of alternate">' >> temp_header.txt
    echo -e '##INFO=<ID=zlr,Number=0,Type=Flag,Description="Zero length allele. Had to be expanded from the reference.">' >> temp_header.txt
    # Add the two headers together and add the column names
    cat temp_header.txt contigs.txt > final_header.txt
    tail -n 1 vcf_header.txt >> final_header.txt
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
    bcftools norm --check-ref e -f genome.fa  variants_remapped_sorted.vcf.gz -o ${outfile_basename} -O v
    """
}

/*
 * Create file containing remapping stats
 */
process outputStats {

    publishDir outfile_dir,
        overwrite: true,
        mode: "copy"

    input:
        path "summary"

    output:
        path "${outfile_basename}.yml"

    """
    ln -s summary "${outfile_basename}.yml"
    """
}

process combineVCF {
    input:
        path "variants1.vcf"
        path "variants2.vcf"
        path "variants3.vcf"
    output:
        path "merge.vcf", emit: merge_vcf

    """
    cat variants1.vcf variants2.vcf variants3.vcf > merge.vcf
    """
}

process combineYaml {
    input:
        path "round1.yml"
        path "round2.yml"
        path "round3.yml"

    output:
        path "merge.yml", emit: merge_yml

    """
    cat round1.yml round2.yml round3.yml > merge.yml
    """
}

// Take variants remapped to the new genome and merge them back with the original header and sort the output
workflow finalise {
    take:
        variants_remapped
        vcf_header
        genome
        summary

    main:
        buildHeader(variants_remapped, vcf_header)
        mergeHeaderAndContent(buildHeader.out.final_header, variants_remapped)
        sortVCF(mergeHeaderAndContent.out.final_vcf_with_header)
        normalise(sortVCF.out.variants_remapped_sorted_gz, genome)
        outputStats(summary)
}


// process_with_minimap
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
        process_split_reads_mid(
            process_split_reads.out.variants_unmapped,
            params.oldgenome,
            prepare_old_genome.out.genome_fai,
            prepare_old_genome.out.genome_chrom_sizes,
            params.newgenome,
            prepare_new_genome.out.genome_fai
        )
        process_split_reads_long(
            process_split_reads_mid.out.variants_unmapped,
            params.oldgenome,
            prepare_old_genome.out.genome_fai,
            prepare_old_genome.out.genome_chrom_sizes,
            params.newgenome,
            prepare_new_genome.out.genome_fai
        )
        combineVCF(
            process_split_reads.out.variants_remapped,
            process_split_reads_mid.out.variants_remapped,
            process_split_reads_long.out.variants_remapped
        )
        combineYaml(
            process_split_reads.out.summary_yml,
            process_split_reads_mid.out.summary_yml,
            process_split_reads_long.out.summary_yml,
        )

        finalise(combineVCF.out.merge_vcf, storeVCFHeader.out.vcf_header, params.newgenome, combineYaml.out.merge_yml)
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
        finalise(process_split_reads_with_bowtie.out.variants_remapped, storeVCFHeader.out.vcf_header,
                 params.newgenome, process_split_reads_with_bowtie.out.summary_yml)
}
