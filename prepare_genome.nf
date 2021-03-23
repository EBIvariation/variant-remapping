#!/usr/bin/env nextflow

// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2

process samtoolsFaidx {

    input:
        path "genome_basename"

    output:
        path "genome_basename.fai", emit: genome_fai

    """
    samtools faidx genome_basename
    """
}

/*
 * Extract chomosome/contig sizes
 */
process chromSizes {

    input:
        path "genome.fa.fai"

    output:
        path "genome.fa.chrom.sizes", emit: genome_chrom_sizes

    """
    cut -f1,2 genome.fa.fai > genome.fa.chrom.sizes 
    """
} 


workflow prepare_old_genome {
    take:
        genome
    main:
        samtoolsFaidx(genome)
        chromSizes(samtoolsFaidx.out.genome_fai)
    emit:
        genome_fai = samtoolsFaidx.out.genome_fai
        genome_chrom_sizes = chromSizes.out.genome_chrom_sizes
}

workflow prepare_new_genome {
    take:
        genome
    main:
        samtoolsFaidx(genome)
    emit:
        genome_fai = samtoolsFaidx.out.genome_fai

}
