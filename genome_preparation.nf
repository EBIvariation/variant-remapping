#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Index the provided genome with bowtie and samtools
    Usage:
            --fasta                     Fasta file containing the reference genome.
            --outdir                    Directory where the index will be deployed.
    """
}

// Show help message
if (params.help) exit 0, helpMessage()

params.fasta  = "$baseDir/resources/genome.fa"
params.outdir = "$baseDir/resources/"

// The basename variable contains the name and extension of the input file, e. g. "genome.fa"
basename = file(params.fasta).getName()
/*
 * Index the provided reference genome using bowtie build
 */
process bowtieGenomeIndex {

    // Memory required is 10 times the size of the fasta in Bytes or at least 1GB
    memory Math.max(file(params.fasta).size() * 10, 1073741824) + ' B'

    publishDir params.outdir,
        overwrite: false,
        mode: "move"

    input:
        path "genome_fasta" from params.fasta

    output:
        path "$basename.*.bt2"

    """
    bowtie2-build genome_fasta $basename
    """
}
/*
 * Index the provided reference genome using samtools faidx
 */
process samtoolsFaidx {

    publishDir params.outdir,
        overwrite: false,
        mode: "copy"

    input:
        path "${basename}" from params.fasta

    output:
        path "${basename}.fai" into genome_fai

    """
    samtools faidx ${basename}
    """
}

/*
 * Extract chomosome/contig sizes
 */
process chromSizes {

    publishDir params.outdir,
        overwrite: false,
        mode: "copy"

    input:
        path "${basename}.fai" from genome_fai

    output:
        path "${basename}.chrom.sizes" into chrom_sizes

    """
    cut -f1,2 ${basename}.fai > ${basename}.chrom.sizes
    """
} 
