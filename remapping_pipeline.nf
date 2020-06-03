#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Please ensure all input files are available
    Inputs:
            VCF file = $baseDir/resources/source.vcf"
            Original genome = "$baseDir/resources/genome.fa
            New genome = "$baseDir/resources/genome.fa"
    Parameters:
            params.vcffile = original vcf file containing variants to be remapped
            params.oldgenome = source genome file used to discover variants from params.vcffile
            params.outfile = final vcf containing variants mapped to new genome
            params.newgenome = new genome which the original vcf file will be mapped against
            params.flankingseq = The length of the flanking sequences that generates reads
            params.scorecutoff = Percentage of the flanking sequences that should be used as the alignment Score cut-off threshold
            params.diffcutoff = Percentage of the flanking sequences that should be used as the AS-XS difference cut-off threshold
    """
}
// Show help message
if (params.help) exit 0, helpMessage()

params.flankingseq = 50
params.scorecutoff = 0.6
params.diffcutoff =  0.04
params.outfile = "$baseDir/resources/remap_vcf.vcf"
params.oldgenome = "$baseDir/resources/genome.fa"
params.newgenome = "$baseDir/resources/genome.fa"
params.vcffile = "$baseDir/resources/source.vcf"

// Index files for both old and new genomes 
oldgenome_chrom_sizes = file("${params.oldgenome}.chrom.sizes")
oldgenome_fai = file("${params.oldgenome}.fai")
newgenome_dir = file(params.newgenome).getParent()
newgenome_basename = file(params.newgenome).getName()
newgenome_fai = file("${params.newgenome}.fai")

// basename and basedir of the output file to know how to name the output file
outfile_basename = file(params.outfile).getName()
outfile_dir = file(params.outfile).getParent()


/*
 * Store the original VCF header for later use
 */
process StoreVCFHeader {

    input:
        path "source.vcf" from params.vcffile

    output:
        path "vcf_header.txt" into vcf_header

    """
    bcftools view --header-only source.vcf -o vcf_header.txt
    """
}

/*
 * Convert VCF file to bed format
 */
process ConvertVCFToBed {

    input:
        path "source.vcf" from params.vcffile

    output:
        path "variants.bed" into variants_bed

    """
    vcf2bed < source.vcf > variants.bed
    """
}

/*
 * Convert get the flanking region in bed format.
 */
process flankingRegionBed {

    input:
        path "variants.bed" from variants_bed
        path "genome.chrom.sizes" from oldgenome_chrom_sizes

    output:
        path "flanking.filtered.bed" into flanking_bed

    script:
    readlength = params.flankingseq * 2 + 1
    """
    bedtools slop -i variants.bed -g genome.chrom.sizes -b ${params.flankingseq} > flanking.bed
    # Check if the intervals are less than the required length (2*flanking length + 1): this filters out reads that aren't 
    # the correct number of bases long (origin: variant too close to the start or end of the chromosome)
    awk  '(\$3 - \$2 == $readlength){print \$0}' flanking.bed > flanking.filtered.bed
    """
}

/*
 * Extract the actual flanking region in fasta format.
 */
process flankingRegionFasta {

    input:  
        path "flanking.bed" from flanking_bed
        path "genome.fa" from params.oldgenome
        path "genome.fa.fai" from oldgenome_fai
    
    output:
        path "variants_reads.fa" into variants_reads

    '''
    # Get the fasta sequences for these intervals
    bedtools getfasta -fi genome.fa -bed flanking.bed -fo variants_reads.fa

    # Replace the colon separators with "|":
    # Storing this information for later on in the script when we split the name by "|" to extract the relevant 
    # information (ALT allele, QUAL, FILT, INFO)
    # This is done because the INFO column can contain ":", which means we wouldn't be able to split by ":", so "|" was 
    # chosen
    sed -i 's/:/|/' variants_reads.fa
    '''
}

/*
 * Extract information about the original variants and put it in the fasta header
 */
process extractVariantInfoToFastaHeader {

    input:  
        path "flanking.bed" from flanking_bed
        path "variants_reads.fa" from variants_reads
    
    output:
        path "variant_reads.out.fa" into variant_reads_with_info

    // Disable the string interpolation using single quotes
    // https://www.nextflow.io/docs/latest/script.html#string-interpolation
    '''
    # Store ref bases
    cut -f 6 flanking.bed > old_ref_bases.txt

    # Store rsIDs
    cut -f 4 flanking.bed > rsIDs.txt

    # Store variant bases
    cut -f 7  flanking.bed > variant_bases.txt

    # Store the other vcf columns
    cut -f 5,8,9 flanking.bed > qual_filt_info.txt

    # Paste the names, variant bases, then fasta sequences into a new file
    paste <(grep '^>' variants_reads.fa) old_ref_bases.txt variant_bases.txt rsIDs.txt qual_filt_info.txt \
    <(grep -v '^>' variants_reads.fa) > temp.txt

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
    sed 's/\\t/|/g; s/ /|/g; s/\\(.*\\)|/\\1 /' temp.txt | tr ' ' '\\n' > variant_reads.out.fa
    '''
}

/*
 * Align sequence with bowtie2
 */
process alignWithBowtie {

    // Memory required is 5 times the size of the fasta in Bytes or at least 1GB
    memory  Math.max(file(params.newgenome).size() * 5, 1073741824) + ' B'

    input:  
        path "variant_reads.fa" from variant_reads_with_info
        // This will get the directory containing the bowtie index linked in the directory
        file 'bowtie_index' from newgenome_dir


    output:
        path "reads_aligned.bam" into reads_aligned_bam

    """
    bowtie2 -k 2 --np 0 -f -x bowtie_index/${newgenome_basename} variant_reads.fa | samtools view -bS - > reads_aligned.bam
    """
}

/*
 * Sort the bam file with samtools
 */
process sortBam {

    input:
        path "reads_aligned.bam" from reads_aligned_bam

    output:
        path "reads_aligned.sorted.bam" into reads_sorted_bam
        path "reads_aligned.sorted.bam.bai" into bam_index

    """
    samtools sort reads_aligned.bam -o reads_aligned.sorted.bam
    samtools index reads_aligned.sorted.bam
    """
}

/*
 * Correct reverse strand alleles
 */
process reverseStrand {

    input:
        path "reads_aligned.sorted.bam" from reads_sorted_bam

    output:
        path "old_ref_alleles.txt" into old_ref_alleles
        path "variants_remapped.vcf" into variants_remapped

    """
    # Ensure that we will use the reverse_strand.py from this repo
    ${baseDir}/reverse_strand.py -i reads_aligned.sorted.bam -p old_ref_alleles.txt -o variants_remapped.vcf -f $params.flankingseq -s $params.scorecutoff -d $params.diffcutoff 
    """
}

/*
 * Insert the reference base from the reference genome. 
 */
process insertReferenceAllele {

    input:
        path "variants_remapped.vcf" from variants_remapped
        path "genome.fa" from params.newgenome
        path "genome.fa.fai" from newgenome_fai

    output:
        path "var_pre_final.vcf" into var_pre_final

    '''
    # Add interval for variant position in bed format (required by getfasta)
    # Reprints all the columns, adding an extra column before the pos column as the pos-1, as bedtools getfasta requires a 
    # bed file
    # Input: 
    # [chr]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
    # Output:
    # [chr]	[pos-1]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
    awk '{print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7}' variants_remapped.vcf > variants_remapped.bed
    
    # Getfasta to get the REF genome alleles and then extract the genome alleles without the positions
    bedtools getfasta -tab -fi genome.fa -bed variants_remapped.bed | awk '{print $2}'> genome_alleles.txt
    
    # Make all alleles caps (sometimes there are lower-case bases)
    tr a-z A-Z < genome_alleles.txt > genome_alleles.fixed.txt
    
    # Paste all the information into the correct columns
    cut -f1-3 variants_remapped.vcf > temp_first_3_columns.txt
    cut -f4-7 variants_remapped.vcf > temp_last_4_columns.txt
    paste temp_first_3_columns.txt genome_alleles.fixed.txt temp_last_4_columns.txt > var_pre_final.vcf
    '''
}

/*
 * Create the header for the output VCF
 */
process buildHeader {

    input:
        path "var_pre_final.vcf" from var_pre_final
        path "vcf_header.txt" from vcf_header

    output:
        path "final_header.txt" into final_header

    """
    # Create list of contigs/chromosomes to be added to the header
    cut -f 1 var_pre_final.vcf | sort -u > contig_names.txt
    while read CHR; do echo "##contig=<ID=\${CHR}>"; done < contig_names.txt > contigs.txt
    # Add the reference assembly
    echo "##reference=${params.newgenome}" >> contigs.txt
    
    # Copy everything that isn't #CHROM (the column titles), ##contig or ##reference from the old header to a temp
    awk '(\$1 !~ /^##contig/ && \$1 !~ /^##reference/ && \$1 !~ /^#CHROM/) {print \$0}' vcf_header.txt > temp_header.txt

    # Add the two headers together and add the column names
    cat temp_header.txt contigs.txt > final_header.txt
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> final_header.txt
    """
}

/*
 * Add header to output VCF and merge with reference allele file
 */
process mergeHeaderAndContent {

    input:
        path "final_header.txt" from final_header
        path "var_pre_final.vcf" from var_pre_final
        path "old_ref_alleles.txt" from old_ref_alleles

    output:
        path "vcf_out_with_header.vcf" into final_vcf_with_header
        path "old_ref_alleles_with_header.txt" into old_ref_alleles_with_header
        
    """
    # Add header to the vcf file:
    cat final_header.txt var_pre_final.vcf > vcf_out_with_header.vcf
    # Also add it to the ref alleles file so that the line numbers match:
    cat final_header.txt old_ref_alleles.txt > old_ref_alleles_with_header.txt
    """
}

/*
 * Fix any variants with matching REF & ALT alleles
 */
process fixRefAllele {

    input:
        path "vcf_out_with_header.vcf" from final_vcf_with_header
        path "old_ref_alleles_with_header.txt" from old_ref_alleles_with_header

    output:
        path "pre_final_vcf.vcf" into pre_final_vcf
        
    """
    # Test each variant to see if REF = ALT, if so, replace the current REF with the corresponding old REF (this is to 
    # deal with variants such as G > G)
    ${baseDir}/replace_refs.py -i vcf_out_with_header.vcf -r old_ref_alleles_with_header.txt -o pre_final_vcf.vcf
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
        path "genome.fa" from params.newgenome
        path "pre_final_vcf.vcf" from pre_final_vcf

    output:
        path "${outfile_basename}" into final_output_vcf        
    """
    bgzip -c pre_final_vcf.vcf > pre_final_vcf.vcf.gz
    bcftools norm -c ws -f genome.fa -N pre_final_vcf.vcf.gz -o ${outfile_basename} -O v
    """
}

/*
 * Create file containing remapping stats
 */
process calculateStats {

    publishDir outfile_dir,
        overwrite: true,
        mode: "copy"

    input:
        path outfile_basename from final_output_vcf

    output:
        path "${outfile_basename}.stats"
        
    """
    bcftools stats ${outfile_basename} > ${outfile_basename}.stats
    """
}
