#!/usr/bin/env nextflow


// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2


/*
 * Convert the VCF file to BED format.
 */
process convertVCFToBed {

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
 * Based on variants BED, generate the flanking regions BED files.
 */
process flankingRegionBed {

    input:
        path "variants.bed"
        path "genome.chrom.sizes"
        val flankingseq

    output:
        path "flanking_r1.bed", emit: flanking_r1_bed
        path "flanking_r2.bed", emit: flanking_r2_bed

    script:
    // The Flanking sequence will start/end one base up/downstream  of the variant.
    // We need to add only (flankingseq - 1) to that base to have the correct flank length
    flankingseq = flankingseq - 1
    """
    awk 'BEGIN{OFS="\t"}{\$2=\$2-1;\$3=\$3-1; print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l $flankingseq -r 0  > flanking_r1.bed

    awk 'BEGIN{OFS="\t"}{\$2=\$2+length(\$6);\$3=\$3+length(\$6); print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l 0 -r $flankingseq  > flanking_r2.bed
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
    # Store variant positions (add + 1 to revert to one based position)
    awk '{print $1"\t"$3 + 1}' flanking_r1.bed > position.txt

    # Store ref bases
    cut -f 6 flanking_r1.bed > old_ref_bases.txt

    # Store rsIDs
    cut -f 4 flanking_r1.bed > rsIDs.txt

    # Store variant bases
    cut -f 7 flanking_r1.bed > variant_bases.txt

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
        path "genome.fa"
        val flanklength

    output:
        path "reads_aligned.bam", emit: reads_aligned_bam


    script:
    if (flanklength < 500)
        """
        # Options used by the 'sr' preset with some modifications:
        # -O6,16 instead of -O12,32 --> reduce indel cost
        # -B5 instead of -B10 --> reduce mismatch cost
        # --end-bonus 20 --> bonus score when the end of the read aligns to mimic global alignment.
        # --secondary=yes -N 2 --> allow up to 2 secondary alignments
        minimap2 -k21 -w11 --sr --frag=yes -A2 -B5 -O6,16 --end-bonus 20 -E2,1 -r50 -p.5 -z 800,200\
                 -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=yes -N 2 \
                 -a genome.fa variant_read1.fa variant_read2.fa | samtools view -bS - > reads_aligned.bam
        """
    else
        """
        minimap2 -k19 -w19 -A2 -B5 -O6,16 --end-bonus 20 -E3,1 -s200 -z200 -N50 --min-occ-floor=100 \
                 --secondary=yes -N 2 \
                 -a genome.fa variant_read1.fa variant_read2.fa | samtools view -bS - > reads_aligned.bam
        """
}

/*
 * Sort BAM file by name
 */
process sortByName {

    input:
        path "reads_aligned.bam"

    output:
        path "reads_aligned_name_sorted.bam", emit: reads_aligned_sorted_bam

    """
    samtools sort -n -O BAM -o reads_aligned_name_sorted.bam reads_aligned.bam
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
        val flank_length
        val filter_align_with_secondary

    output:
        path "variants_remapped.vcf", emit: variants_remapped
        path "variants_unmapped.vcf", emit: variants_unmapped
        path "summary.yml", emit: summary_yml

    script:
        if (filter_align_with_secondary)
            """
            # Ensure that we will use the reads_to_remapped_variants.py from this repo
            ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads_aligned.bam \
                -o variants_remapped.vcf  --newgenome genome.fa --out_failed_file variants_unmapped.vcf \
                --flank_length $flank_length --summary summary.yml --filter_align_with_secondary
            """
        else
            """
            # Ensure that we will use the reads_to_remapped_variants.py from this repo
            ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads_aligned.bam \
                -o variants_remapped.vcf  --newgenome genome.fa --out_failed_file variants_unmapped.vcf \
                --flank_length $flank_length --summary summary.yml
           """

}

workflow process_split_reads_generic {
    take:
        source_vcf
        old_genome_fa
        old_genome_fa_fai
        old_genome_chrom_sizes
        new_genome_fa
        new_genome_fa_fai
        flank_length
        filter_align_with_secondary

    main:
        convertVCFToBed(source_vcf)
        flankingRegionBed(convertVCFToBed.out.variants_bed, old_genome_chrom_sizes, flank_length)
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
            new_genome_fa,
            flank_length
        )
        sortByName(alignWithMinimap.out.reads_aligned_bam)
        readsToRemappedVariants(
            sortByName.out.reads_aligned_sorted_bam, new_genome_fa,
            flank_length, filter_align_with_secondary
        )

    emit:
        variants_remapped = readsToRemappedVariants.out.variants_remapped
        variants_unmapped = readsToRemappedVariants.out.variants_unmapped
        summary_yml = readsToRemappedVariants.out.summary_yml
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
        flank_length = 50
        filter_align_with_secondary = true
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary
        )
    emit:
        variants_remapped = process_split_reads_generic.out.variants_remapped
        variants_unmapped = process_split_reads_generic.out.variants_unmapped
        summary_yml = process_split_reads_generic.out.summary_yml
}


workflow process_split_reads_mid {
    take:
        source_vcf
        old_genome_fa
        old_genome_fa_fai
        old_genome_chrom_sizes
        new_genome_fa
        new_genome_fa_fai

    main:
        flank_length = 2000
        filter_align_with_secondary = true
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary
        )
    emit:
        variants_remapped = process_split_reads_generic.out.variants_remapped
        variants_unmapped = process_split_reads_generic.out.variants_unmapped
        summary_yml = process_split_reads_generic.out.summary_yml
}


workflow process_split_reads_long {
    take:
        source_vcf
        old_genome_fa
        old_genome_fa_fai
        old_genome_chrom_sizes
        new_genome_fa
        new_genome_fa_fai

    main:
        flank_length = 50000
        filter_align_with_secondary = false
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary
        )
    emit:
        variants_remapped = process_split_reads_generic.out.variants_remapped
        variants_unmapped = process_split_reads_generic.out.variants_unmapped
        summary_yml = process_split_reads_generic.out.summary_yml
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
        flank_length = 50
        convertVCFToBed(source_vcf)
        flankingRegionBed(convertVCFToBed.out.variants_bed, old_genome_chrom_sizes, flank_length)
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
        readsToRemappedVariants(alignWithBowtie.out.reads_aligned_bam, new_genome_fa, flank_length, true)

    emit:
        variants_remapped = readsToRemappedVariants.out.variants_remapped
        variants_unmapped = readsToRemappedVariants.out.variants_unmapped
        summary_yml = readsToRemappedVariants.out.summary_yml
}
