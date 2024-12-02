#!/usr/bin/env nextflow

// Enable syntax extension
// See https://www.nextflow.io/docs/latest/dsl2.html
nextflow.enable.dsl=2


/*
 * Convert the VCF file to BED format storing the VCF line in the "name" column and the reference allele in the
 * "strand" column.
 */
process convertVCFToBed {
    label 'default_time', 'med_mem'

    input:
        path "source.vcf"

    output:
        path "variants.bed", emit: variants_bed

    '''
    # Convert the vcf file to bed format:
    #  - Remove all headers
    #  - Switch to 0 based coordinates system
    #  - Add the reference allele so it can be used in flankingRegionBed to adjust the position of the right flank
    #  - add all VCF fields separated by 2 characters pipe and caret (|^) to avoid impacting existing formatting of
    #    the VCF line. The sub replacing percent is to protect the % character that would be interpreted by printf
    #    otherwise. the sub replacing space is to prevent bedtools from using them as a field separator
    awk -F '\\t' '{ if (!/^#/){ \
                    printf $1"\\t"$2-1"\\t"$2"\\t"$1; \
                    for (i=2; i<=NF; i++){ gsub(/%/, "%%", $i); gsub(/ /, "£€", $i); printf "|^"$i }; print "\\t"$4}; \
                  }' source.vcf \
                  > variants.bed
    '''
}

/*
 * Based on variants BED, generate the BED file for each flank.
 */
process flankingRegionBed {
    label 'default_time', 'med_mem'

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
    # Adjust the end position of the flank to be one base upstream of the variant
    awk 'BEGIN{OFS="\\t"}{\$2=\$2-1; \$3=\$3-1; print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l $flankingseq -r 0  > flanking_r1.bed

    # Adjust the start position of the flank to be one base downstream of the end of variant (\$5 is the reference allele)
    awk 'BEGIN{OFS="\\t"}{ \$2=\$2+length(\$5); \$3=\$3+length(\$5); print \$0}' variants.bed \
        | bedtools slop  -g genome.chrom.sizes -l 0 -r $flankingseq  > flanking_r2.bed
    """
}

/*
 * Extract the actual flanking region in fasta format.
 */
process flankingRegionFasta {
    label 'default_time', 'med_mem'

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
    '''
}

/*
 * Extract information about the original variants and put it in the fasta header
 */
process extractVariantInfoToFastaHeader {
    label 'default_time', 'med_mem'

    input:  
        path "flanking_r1.bed"
        path "flanking_r2.bed"
        path "variants_read1.fa"
        path "variants_read2.fa"

    output:
        path "interleaved.fa", emit: interleaved_fasta

    // Disable Nextflow string interpolation using single quotes
    // https://www.nextflow.io/docs/latest/script.html#string-interpolation
    '''
    # Store variant position in the file to have a unique name
    awk '{print ">" NR }' flanking_r1.bed > position.txt

    # Store position of the variant in the file and replace '£€' with the original whitespace from convertVCFToBed
    cut -f 4 flanking_r1.bed | sed 's/£€/ /g' > vcf_fields.txt

    # Paste the names, variant bases, then fasta sequences into a new file
    # A space will be inserted between the position and the vcf fields
    # Then a newline is inserted between the vcf fields and the sequence
    # The vcf fields are regarded as comment to the fasta entry.
    paste -d ' \\n' position.txt vcf_fields.txt <(grep -v '^>' variants_read1.fa) > variant_read1.out.fa
    paste -d '\\n' position.txt <(grep -v '^>' variants_read2.fa) > variant_read2.out.fa

    paste variant_read1.out.fa variant_read2.out.fa | paste - - | awk -F "\\t" 'BEGIN {OFS="\\n"} {print $1,$3,$2,$4}' > interleaved.fa
    '''
}

/*
 * Split fasta entries into multiple chunks
 */
process split_fasta {
    label 'short_time', 'small_mem'

    input:
        path interleaved_fasta
        val chunk_size

    output:
        path("read_chunk-*", emit: read_split)

    script:
    if (interleaved_fasta.size() > 0)
        """
        split -a 5 -d -l ${chunk_size * 4} ${interleaved_fasta} read_chunk-
        """
    else
        """
        ln -s ${interleaved_fasta} read_chunk-00001
        """
}

/*
 * Align sequence with minimap2
 */
process alignWithMinimap {
    label 'med_time'

    // Memory required is 10 times the size of the fasta in Bytes or at least 2GB
    // Overwrite base_memory so that the standard retry strategy is used
    ext base_memory: { Math.max(file(params.newgenome).size() * 10, 2000000000) }

    input:
        // reads contains paired interleaved (first and second read in the same file)
        each path(reads)
        // indexing is done on the fly so get the genome directly
        path "genome.fa"
        val flanklength

    output:
        path "reads_aligned.bam", emit: reads_aligned_bam

    script:
    index_size = ""
    if (file(params.newgenome).size() > 4000000000){
        index_size = " -I " + file(params.newgenome).size() * 1.1
    }
    if (flanklength < 500)
        """
        set -o pipefail
        # Options used by the 'sr' preset with some modifications:
        # -O6,16 instead of -O12,32 --> reduce indel cost
        # -B5 instead of -B10 --> reduce mismatch cost
        # --end-bonus 20 --> bonus score when the end of the read aligns to mimic global alignment.
        # --secondary=yes -N 2 --> allow up to 2 secondary alignments
        # -y option will take the comment from the fasta entry and output it
        # the awk script will convert this comment in valid SAM tag
        minimap2 $index_size -k21 -w11 --sr --frag=yes -A2 -B5 -O6,16 --end-bonus 20 -E2,1 -r50 -p.5 -z 800,200\
                 -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=yes -N 2 -y \
                 -a genome.fa ${reads} | \
                 awk -F '\\t' 'BEGIN{OFS="\\t"}{if(!/^@/){\$NF="vr:Z:"\$NF}; print \$0;}' | \
                 samtools view -bS - > reads_aligned.bam
        """
    else
        """
        set -o pipefail
        minimap2 $index_size -k19 -w19 -A2 -B5 -O6,16 --end-bonus 20 -E3,1 -s200 -z200 -N50 --min-occ-floor=100 \
                 --secondary=yes -N 2 -y \
                 -a genome.fa ${reads} | \
                 awk -F '\\t' 'BEGIN{OFS="\\t"}{if(!/^@/){\$NF="vr:Z:"\$NF}; print \$0;}' | \
                 samtools view -bS - > reads_aligned.bam
        """
}

/*
 * Sort BAM file by name
 */
process sortByName {
    label 'default_time', 'med_mem'

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
    label 'med_time'

    // Memory required is 5 times the size of the fasta in Bytes or at least 1GB
    // Overwrite base_memory so that the standard retry strategy is used
    ext base_memory: { Math.max(file(params.newgenome).size() * 5, 1073741824) }

    input:
        path "variant_read1.fa"
        path "variant_read2.fa"
        // This will get all the index files named as they were placed in a subdirectory
        file "bowtie_index/*"

    output:
        path "reads_aligned.bam", emit: reads_aligned_bam


    """
    bowtie2 -k 2 --end-to-end --np 0 --sam-append-comment -f -x bowtie_index/bowtie_index \
      -1 variant_read1.fa -2 variant_read2.fa \
      | awk -F '\\t' 'BEGIN{OFS="\\t"}{if(!/^@/){\$NF="vr:Z:"\$NF}; print \$0;}' \
      | samtools view -bS - > reads_aligned.bam
    """
}


/*
 * Take the reads and process them to get the remapped variants
 */
process readsToRemappedVariants {
    label 'default_time', 'med_mem'

    input:
        path "reads.bam"
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
            ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads.bam \
                -o variants_remapped.vcf  --newgenome genome.fa --out_failed_file variants_unmapped.vcf \
                --flank_length $flank_length --summary summary.yml --filter_align_with_secondary
            """
        else
            """
            # Ensure that we will use the reads_to_remapped_variants.py from this repo
            ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads*.bam \
                -o variants_remapped.vcf  --newgenome genome.fa --out_failed_file variants_unmapped.vcf \
                --flank_length $flank_length --summary summary.yml
           """
}

/*
 * Gather step for remapped and unmapped variants and the summary yaml file
 *
 */
process merge_variants {
    label 'short_time', 'small_mem'

    input:
        path "remapped*.vcf"
        path "unmapped*.vcf"
        path "summary*.yml"

    output:
       path "variants_remapped.vcf", emit: variants_remapped
       path "variants_unmapped.vcf", emit: variants_unmapped
       path "output_summary.yml", emit: summary_yml

    """
    cat remapped*.vcf > variants_remapped.vcf
    cat unmapped*.vcf > variants_unmapped.vcf
    ${baseDir}/variant_remapping_tools/merge_yaml.py --input summary*.yml --output output_summary.yml
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
        chunk_size

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

        split_fasta(
            extractVariantInfoToFastaHeader.out.interleaved_fasta,
            chunk_size
        )


        alignWithMinimap(
            split_fasta.out.read_split,
            new_genome_fa,
            flank_length
        )
        sortByName(alignWithMinimap.out.reads_aligned_bam)
        // Collect all the bam files in the next step
        readsToRemappedVariants(
            sortByName.out.reads_aligned_sorted_bam, new_genome_fa,
            flank_length, filter_align_with_secondary
        )
        merge_variants(
            readsToRemappedVariants.out.variants_remapped.collect(),
            readsToRemappedVariants.out.variants_unmapped.collect(),
            readsToRemappedVariants.out.summary_yml.collect()
        )

    emit:
        variants_remapped = merge_variants.out.variants_remapped
        variants_unmapped = merge_variants.out.variants_unmapped
        summary_yml = merge_variants.out.summary_yml
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
        chunk_size = 5000000
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary,
            chunk_size
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
        chunk_size = 500000
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary,
            chunk_size
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
        chunk_size = 50000
        process_split_reads_generic(
            source_vcf, old_genome_fa, old_genome_fa_fai, old_genome_chrom_sizes,
            new_genome_fa, new_genome_fa_fai, flank_length, filter_align_with_secondary,
            chunk_size
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
