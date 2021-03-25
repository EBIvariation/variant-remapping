#!/usr/bin/env nextflow

// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2



/*
* Index the new reference genome using bowtie_build
*/
process bowtieGenomeIndex {
    // Memory required is 10 times the size of the fasta in Bytes or at least 1GB
    memory Math.max(file(params.newgenome).size() * 10, 1073741824) + ' B'

    input:
        path "genome_fasta"

    output:
        path "bowtie_index.*.bt2", emit: bowtie_indexes

    """
    bowtie2-build genome_fasta bowtie_index
    """
}


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


workflow prepare_new_genome_bowtie {
    take:
        genome
    main:
        samtoolsFaidx(genome)
        bowtieGenomeIndex(genome)
    emit:
        genome_fai = samtoolsFaidx.out.genome_fai
        bowtie_indexes = bowtieGenomeIndex.out.bowtie_indexes
}