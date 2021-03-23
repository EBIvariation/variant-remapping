#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Process a VCF file containing SNV or short indels associated with the old genome to remap the variant's coordinates to the new genome.

    Inputs:
            --vcffile      original vcf file containing variants to be remapped [required]
            --oldgenome    source genome file used to discover variants from the file provided with --vcffile [required]
            --newgenome    new genome which the original vcf file will be mapped against [required]
    Outputs:
            --outfile      final vcf containing variants mapped to new genome [required]
    Parameters:
            --flankingseq  The length of the flanking sequences that generates reads (default 50)
            --scorecutoff  Percentage of the flanking sequences that should be used as the alignment score cut-off threshold (default 0.6)
            --diffcutoff   Percentage of the flanking sequences that should be used as the AS-XS difference cut-off threshold (default 0.04)
    """
}

params.flankingseq = 50
params.scorecutoff = 0.6
params.diffcutoff =  0.04
params.outfile = null
params.oldgenome = null
params.newgenome = null
params.vcffile = null
params.help = null

// Show help message
if (params.help) exit 0, helpMessage()

// Test input files
if (!params.vcffile || !params.oldgenome || !params.outfile || !params.newgenome) { 
    if (!params.vcffile)    log.warn('Provide a input vcf file using --vcffile')
    if (!params.oldgenome)  log.warn('Provide a fasta file for the old genome --oldgenome') 
    if (!params.outfile)    log.warn('Provide a path to the output vcf file using --outfile') 
    if (!params.newgenome)  log.warn('Provide a fasta file for the new genome file using --newgenome') 
    exit 1, helpMessage()
}
 
// Index files for both old and new genomes 
oldgenome_dir = file(params.oldgenome).getParent()
oldgenome_basename = file(params.oldgenome).getName()
oldgenome_fai = file("${params.oldgenome}.fai")
newgenome_dir = file(params.newgenome).getParent()
newgenome_basename = file(params.newgenome).getName()
newgenome_fai = file("${params.newgenome}.fai")

// basename and basedir of the output file to know how to name the output file
outfile_basename = file(params.outfile).getName()
outfile_dir = file(params.outfile).getParent()


/*
* Check that the fai index file for old genome exists and if it does not create it.
* Once created the publishDir directive will place it in the location described by the oldgenome_fai variable
*/
if (!oldgenome_fai.exists()){
    process samtoolsFaidxOld {

        publishDir oldgenome_dir,
            overwrite: false,
            mode: "copy"

        input:
            path "${oldgenome_basename}" from params.oldgenome

        output:
            path "${oldgenome_basename}.fai" into oldgenome_fai

        """
        samtools faidx ${oldgenome_basename}
        """
    }
}

/*
* Check that the fai index file for new genome exists and if it does not create it.
* Once created the publishDir directive will place it in the location described by the newgenome_fai variable
*/
if (!newgenome_fai.exists()){
    process samtoolsFaidxNew {

        publishDir newgenome_dir,
            overwrite: false,
            mode: "copy"

        input:
            path "${newgenome_basename}" from params.newgenome

        output:
            path "${newgenome_basename}.fai" into newgenome_fai

        """
        samtools faidx ${newgenome_basename}
        """
    }
}

/*
 * Extract chomosome/contig sizes
 */
process chromSizes {

    input:
        path "genome.fa.fai" from oldgenome_fai


    output:
        path "genome.fa.chrom.sizes" into oldgenome_chrom_sizes

    """
    cut -f1,2 genome.fa.fai > genome.fa.chrom.sizes 
    """
} 


/*
 * Uncompress VCF file
 */
process uncompressInputVCF {

    input:
        path "source.vcf" from params.vcffile

    output:
        path "uncompressed.vcf" into vcf_file

    script:
        if ( file(params.vcffile).getExtension() == 'gz' )
            """
            gunzip -c source.vcf > uncompressed.vcf
            """
        else
            """
            ln -nfs source.vcf uncompressed.vcf
            """
}
/*
 * Store the original VCF header for later use
 */
process StoreVCFHeader {

    input:
        path "source.vcf" from vcf_file

    output:
        path "vcf_header.txt" into vcf_header

    """
    bcftools view --header-only source.vcf | grep -v '^##FORMAT' >  vcf_header.txt
    """
}

/*
 * Convert VCF file to bed format
 */
process ConvertVCFToBed {

    input:
        path "source.vcf" from vcf_file

    output:
        path "variants.bed" into variants_bed

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
        path "variants.bed" from variants_bed
        path "genome.chrom.sizes" from oldgenome_chrom_sizes

    output:
        path "flanking_r1.bed" into flanking_r1_bed
        path "flanking_r2.bed" into flanking_r2_bed

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
        path "flanking_r1.bed" from flanking_r1_bed
        path "flanking_r2.bed" from flanking_r2_bed
        path "genome.fa" from params.oldgenome
        path "genome.fa.fai" from oldgenome_fai
    
    output:
        path "variants_read1.fa" into variants_read1
        path "variants_read2.fa" into variants_read2

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
        path "flanking_r1.bed" from flanking_r1_bed
        path "flanking_r2.bed" from flanking_r2_bed
        path "variants_read1.fa" from variants_read1
        path "variants_read2.fa" from variants_read2

    output:
        path "variant_read1.out.fa" into variant_read1_with_info
        path "variant_read2.out.fa" into variant_read2_with_info

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
        path "variant_read1.fa" from variant_read1_with_info
        path "variant_read2.fa" from variant_read2_with_info
        // indexing is done on the fly so get the genome directly
        path 'genome.fa' from params.newgenome

    output:
        path "reads_aligned.bam" into reads_aligned_bam

    """
    # Options used by the 'sr' preset but allowing secondary alignments
    minimap2 -k21 -w11 --sr --frag=yes -A2 -B5 --end-bonus 10 -O12,32 -E2,1 -r50 -p.5 -z 800,200\
             -f1000,5000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes --secondary=yes -N 2 \
             -a genome.fa variant_read1.fa variant_read2.fa | samtools view -bS - > reads_aligned.bam
    """
}



/*
 * Take the reads and process them to get the remapped variants
 *
 */
process readsToRemappedVariants {

    input:
        path "reads_aligned.bam" from reads_aligned_bam
        path "genome.fa" from params.newgenome

    output:
        path "variants_remapped.vcf" into variants_remapped

    """
    # Ensure that we will use the reads_to_remapped_variants.py from this repo
    ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads_aligned.bam \
        -o variants_remapped.vcf  --filter_align_with_secondary \
        --newgenome genome.fa
    """
}

/*
 * Sort VCF file
 */
process sortVCF {

    input:
        path "variants_remapped.vcf" from variants_remapped

    output:
        path "variants_remapped_sorted.vcf" into variants_remapped_sorted

    """
    cat variants_remapped.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > variants_remapped_sorted.vcf
    """
}

/*
 * Create the header for the output VCF
 */
process buildHeader {

    input:
        path "variants_remapped_sorted.vcf" from variants_remapped_sorted
        path "vcf_header.txt" from vcf_header

    output:
        path "final_header.txt" into final_header

    """
    # Create list of contigs/chromosomes to be added to the header
    cut -f 1 variants_remapped_sorted.vcf | sort -u > contig_names.txt
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
        path "variants_remapped_sorted.vcf" from variants_remapped_sorted

    output:
        path "vcf_out_with_header.vcf" into final_vcf_with_header

    """
    # Add header to the vcf file:
    cat final_header.txt variants_remapped_sorted.vcf > vcf_out_with_header.vcf
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
        path "vcf_out_with_header.vcf" from final_vcf_with_header

    output:
        path "${outfile_basename}" into final_output_vcf        

    """
    bgzip -c vcf_out_with_header.vcf > vcf_out_with_header.vcf.gz
    bcftools norm -c ws -f genome.fa -N vcf_out_with_header.vcf.gz -o ${outfile_basename} -O v
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
