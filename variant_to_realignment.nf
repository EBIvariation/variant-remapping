#!/usr/bin/env nextflow


// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2


/*
 * Convert the VCF file to BED format.
 */
process ConvertVCFToBed {

    input:
        path "source.vcf"

    output:
        path "variants.bed", emit: variants_bed

    """
    # TODO: Change vcf2bed so it does not split the alternates into two lines
    vcf2bed < source.vcf > variants.bed
    """
}

/*
 * Based on variants BED, generate the flanking regions BED.
 */
process flankingRegionBed {

    input:
        path "variants.bed"
        path "genome.chrom.sizes"

    output:
        path "flanking_r1.bed", emit: flanking_r1_bed
        path "flanking_r2.bed", emit: flanking_r2_bed

    script:
    """
    awk 'BEGIN{OFS="\t"}{\$2=\$2-1;\$3=\$3-1; print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l $params.flankingseq -r 0  > flanking_r1.bed

    awk 'BEGIN{OFS="\t"}{\$2=\$2+length(\$6);\$3=\$3+length(\$6); print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l 0 -r $params.flankingseq  > flanking_r2.bed
    """
}

/*
 * Extract the actual flanking region in fasta format.
 */
process flankingRegionFasta {

    memory '4 GB'

    input:  
        path "flanking_r1.bed"
        path "flanking_r2.bed"
        path "genome.fa"
        path "genome.fa.fai"
    
    output:
        path "variants_read1.fa", emit: variants_read1
        path "variants_read2.fa", emit: variants_read2

    '''
    # Get the fasta sequences for these intervals
    bedtools getfasta -fi genome.fa -bed flanking_r1.bed -fo variants_read1.fa
    bedtools getfasta -fi genome.fa -bed flanking_r2.bed -fo variants_read2.fa

    # Replace the colon separators with "|":
    # Storing this information for later on in the script when we split the name by "|" to extract the relevant 
    # information (ALT allele, QUAL, FILT, INFO)
    # This is done because the INFO column can contain ":", which means we wouldn't be able to split by ":", so "|" was 
    # chosen
    sed -i 's/:/|/' variants_read1.fa
    sed -i 's/:/|/' variants_read2.fa
    '''
}

/*
 * Extract information about the original variants and put it in the fasta header
 */
process extractVariantInfoToFastaHeader {

    input:  
        path "flanking_r1.bed"
        path "flanking_r2.bed"
        path "variants_read1.fa"
        path "variants_read2.fa"

    output:
        path "variant_read1.out.fa", emit: variant_read1_with_info
        path "variant_read2.out.fa", emit: variant_read2_with_info

    // Disable the string interpolation using single quotes
    // https://www.nextflow.io/docs/latest/script.html#string-interpolation
    '''
    # Store variant positions
    cut -f 1,3 flanking_r1.bed > position.txt

    # Store ref bases
    cut -f 6 flanking_r1.bed > old_ref_bases.txt

    # Store rsIDs
    cut -f 4 flanking_r1.bed > rsIDs.txt

    # Store variant bases
    cut -f 7  flanking_r1.bed > variant_bases.txt

    # Store the other vcf columns
    cut -f 5,8,9 flanking_r1.bed > qual_filt_info.txt

    # Paste the names, variant bases, then fasta sequences into a new file
    paste position.txt old_ref_bases.txt variant_bases.txt rsIDs.txt qual_filt_info.txt \
    <(grep -v '^>' variants_read1.fa) | awk '{print ">"$0}' > temp1.txt

    paste position.txt old_ref_bases.txt variant_bases.txt rsIDs.txt qual_filt_info.txt \
    <(grep -v '^>' variants_read2.fa) | awk '{print ">"\$0}' > temp2.txt

    # Reformat the fasta ID: inconsistencies in the separators, and no new line before the sequence mean that this next 
    # command is a bit ugly
    # Input:
    # >[chr]|[pos interval]   [REF]       [ALT]       [rsID]    [QUAL] [FILT] [INFO]      [seq]
    # Replace all spaces and tabs with "|":
    # >[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO]|[seq]
    # Replace the last "|" with a space (this is between the header and the sequence):
    # >[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO] [seq]
    # And finally replace the space with a newline:
    # Output:
    # >[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO]
    # [seq]
    sed 's/\\t/|/g; s/ /|/g; s/\\(.*\\)|/\\1 /' temp1.txt | tr ' ' '\\n' > variant_read1.out.fa
    sed 's/\\t/|/g; s/ /|/g; s/\\(.*\\)|/\\1 /' temp2.txt | tr ' ' '\\n' > variant_read2.out.fa
    '''
}

/*
 * Align sequence with minimap2
 */
process alignWithMinimap {

    // Memory required is 5 times the size of the fasta in Bytes or at least 1GB
    memory Math.max(file(params.newgenome).size() * 5, 1073741824) + ' B'

    input:
        path "variant_read1.fa"
        path "variant_read2.fa"
        // indexing is done on the fly so get the genome directly
        path 'genome.fa'

    output:
        path "reads_aligned.bam", emit: reads_aligned_bam

    """
    # Options used by the 'sr' preset but allowing secondary alignments
    minimap2 -k21 -w11 --sr --frag=yes -A2 -B5 -O6,16 --end-bonus 20 -E2,1 -r50 -p.5 -z 800,200\
             -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=yes -N 2 \
             -a genome.fa variant_read1.fa variant_read2.fa | samtools view -bS - > reads_aligned.bam
    """
}

/*
 * Align sequence with bowtie2
 */
process alignWithBowtie {

    // Memory required is 5 times the size of the fasta in Bytes or at least 1GB
    memory Math.max(file(params.newgenome).size() * 5, 1073741824) + ' B'

    input:
        path "variant_read1.fa"
        path "variant_read2.fa"
        // This will get all the index files named as they were placed in a subdirectory
        file "bowtie_index/*"

    output:
        path "reads_aligned.bam", emit: reads_aligned_bam


    """
    bowtie2 -k 2 --end-to-end --np 0 -f -x bowtie_index/bowtie_index -1 variant_read1.fa -2 variant_read2.fa | samtools view -bS - > reads_aligned.bam
    """
}


/*
 * Take the reads and process them to get the remapped variants
 *
 */
process readsToRemappedVariants {

    input:
        path "reads_aligned.bam"
        path "genome.fa"

    output:
        path "variants_remapped.vcf", emit: variants_remapped
        path "variants_unmapped.vcf", emit: variants_unmapped

    """
    # Ensure that we will use the reads_to_remapped_variants.py from this repo
    ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads_aligned.bam \
        -o variants_remapped.vcf  --filter_align_with_secondary \
        --newgenome genome.fa --out_failed_file variants_unmapped.vcf
    """
}

workflow process_split_reads {
    take:
        source_vcf
        old_genome_fa
        old_genome_fa_fai
        old_genome_chrom_sizes
        new_genome_fa
        new_genome_fa_fai

    main:
        ConvertVCFToBed(source_vcf)
        flankingRegionBed(ConvertVCFToBed.out.variants_bed, old_genome_chrom_sizes)
        flankingRegionFasta(
            flankingRegionBed.out.flanking_r1_bed, flankingRegionBed.out.flanking_r2_bed,
            old_genome_fa, old_genome_fa_fai
        )
        extractVariantInfoToFastaHeader(
            flankingRegionBed.out.flanking_r1_bed, flankingRegionBed.out.flanking_r2_bed,
            flankingRegionFasta.out.variants_read1, flankingRegionFasta.out.variants_read2
        )
        alignWithMinimap(
            extractVariantInfoToFastaHeader.out.variant_read1_with_info,
            extractVariantInfoToFastaHeader.out.variant_read2_with_info,
            new_genome_fa
        )
        readsToRemappedVariants(alignWithMinimap.out.reads_aligned_bam, new_genome_fa)

    emit:
        variants_remapped = readsToRemappedVariants.out.variants_remapped
        variants_unmapped = readsToRemappedVariants.out.variants_unmapped
}

workflow process_split_reads_with_bowtie {
    take:
        source_vcf
        old_genome_fa
        old_genome_fa_fai
        old_genome_chrom_sizes
        new_genome_fa
        new_genome_fa_fai
        new_genome_bowtie_index

    main:
        ConvertVCFToBed(source_vcf)
        flankingRegionBed(ConvertVCFToBed.out.variants_bed, old_genome_chrom_sizes)
        flankingRegionFasta(
            flankingRegionBed.out.flanking_r1_bed, flankingRegionBed.out.flanking_r2_bed,
            old_genome_fa, old_genome_fa_fai
        )
        extractVariantInfoToFastaHeader(
            flankingRegionBed.out.flanking_r1_bed, flankingRegionBed.out.flanking_r2_bed,
            flankingRegionFasta.out.variants_read1, flankingRegionFasta.out.variants_read2
        )
        alignWithBowtie(
            extractVariantInfoToFastaHeader.out.variant_read1_with_info,
            extractVariantInfoToFastaHeader.out.variant_read2_with_info,
            new_genome_bowtie_index
        )
        readsToRemappedVariants(alignWithBowtie.out.reads_aligned_bam, new_genome_fa)

    emit:
        variants_remapped = readsToRemappedVariants.out.variants_remapped
        variants_unmapped = readsToRemappedVariants.out.variants_unmapped
}


