#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Index the provided genome with bowtie and samtools
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
        mode: "move",
        saveAs: { filename -> file(params.fasta).getName() + filename.replaceFirst(/reference/, "") }

    input:
        path 'reference' from params.fasta

    output:
        path "reference.*.bt2"

    """
    bowtie2-build reference reference
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
        path 'reference' from params.fasta

    output:
        path "reference.fai" into genome_fai

    """
    samtools faidx reference
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
        path 'reference.fai' from genome_fai

    output:
        path "reference.chrom.sizes" into chrom_sizes

    """
    cut -f1,2 reference.fai > reference.chrom.sizes
    """
}
