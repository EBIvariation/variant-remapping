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

/*
 * Index the provided reference genome using bowtie build
 */
process bowtieGenomeIndex {

    publishDir params.outdir,
        overwrite: false,
        mode: "move",
        saveAs: { filename -> file(params.fasta).getName() + filename.replaceFirst(/index/, "") }

    input:
        path 'genome_fasta' from params.fasta

    output:
        path "index.*.bt2"

    """
    bowtie2-build genome_fasta index
    """
}


/*
 * Index the provided reference genome using samtools faidx
 */
process samtoolsFaidx {

    publishDir params.outdir,
        overwrite: false,
        mode: "copy",
        saveAs: { filename -> file(params.fasta).getName() + filename.replaceFirst(/reference/, "") }

    input:
        path 'genome_fasta' from params.fasta

    output:
        path "genome_fasta.fai" into genome_fai

    """
    samtools faidx genome_fasta
    """
}

/*
 * Extract chomosome/contig sizes
 */
process chromSizes {

    publishDir params.outdir,
        overwrite: false,
        mode: "copy",
        saveAs: { filename -> file(params.fasta).getName() + filename.replaceFirst(/reference/, "") }

    input:
        path 'genome_fasta.fai' from genome_fai

    output:
        path "genome_fasta.chrom.sizes" into chrom_sizes

    """
    cut -f1,2 genome_fasta.fai > genome_fasta.chrom.sizes
    """
}
