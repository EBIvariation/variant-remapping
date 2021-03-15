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
* Index the new reference genome using bowtie_build
*/
if (!file("${newgenome_dir}/${newgenome_basename}.1.bt2").exists()){

    process bowtieGenomeIndex {
        // Memory required is 10 times the size of the fasta in Bytes or at least 1GB
        memory Math.max(file(params.newgenome).size() * 10, 1073741824) + ' B'

        publishDir newgenome_dir,
            overwrite: false,
            mode: "move"

        input:
            path "genome_fasta" from params.newgenome

        output:
            path "$newgenome_basename.*.bt2" 
            val "done" into bt2_done

        """
        bowtie2-build genome_fasta $newgenome_basename
        """
    }
} else {
    bt2_done = Channel.value("done")
}

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
    memory Math.max(file(params.newgenome).size() * 5, 1073741824) + ' B'

    input:
        path "variant_reads.fa" from variant_reads_with_info
        // This will get the directory containing the bowtie index linked in the directory
        file 'bowtie_index' from newgenome_dir
        // This ensures that bowtie index exists
        val "unused" from bt2_done

    output:
        path "reads_aligned.bam" into reads_aligned_bam

    """
    bowtie2 -k 2 --np 0 -f -x bowtie_index/${newgenome_basename} variant_reads.fa | samtools view -bS - > reads_aligned.bam
    """
}

/*
 * Sort the bam file with samtools and index the result.
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
        path "genome.fa" from params.newgenome

    output:
        path "variants_remapped.vcf" into variants_remapped

    """
    # Ensure that we will use the reverse_strand.py from this repo
    ${baseDir}/variant_remapping_tools/reads_to_remapped_variants.py -i reads_aligned.sorted.bam \
        -o variants_remapped.vcf -f $params.flankingseq \
        -s $params.scorecutoff -d $params.diffcutoff --newgenome genome.fa
    """
}


/*
 * Create the header for the output VCF
 */
process buildHeader {

    input:
        path "variants_remapped.vcf" from variants_remapped
        path "vcf_header.txt" from vcf_header

    output:
        path "final_header.txt" into final_header

    """
    # Create list of contigs/chromosomes to be added to the header
    cut -f 1 variants_remapped.vcf | sort -u > contig_names.txt
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
        path "variants_remapped.vcf" from variants_remapped

    output:
        path "vcf_out_with_header.vcf" into final_vcf_with_header

    """
    # Add header to the vcf file:
    cat final_header.txt variants_remapped.vcf > vcf_out_with_header.vcf
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
