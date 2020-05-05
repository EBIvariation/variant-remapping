#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Index the provided genome with bowtie
    Usage:
            --fasta                     Fasta file containing the reference genome.
            --outdir                    Directory where th index will be deployed.
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
        mode: "move"

    input:
    path 'reference' from params.fasta

    output:
        file("${reference}.*") into bowtieIndexes

    """
    bowtie2-build reference reference
    """
}
