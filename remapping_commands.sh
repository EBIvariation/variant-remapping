#!/bin/bash

oldgenome=''
newgenome=''
newgenomeaccession=''
vcffile=''
flankingseq=''
outfile=''

while getopts 'g:n:a:v:f:o:' flag; do
	case "${flag}" in
		g) oldgenome="${OPTARG}";;
		n) newgenome="${OPTARG}" ;;
		a) newgenomeaccession="${OPTARG}" ;;
		v) vcffile="${OPTARG}" ;;
		f) flankingseq="${OPTARG}" ;;
		o) outfile="${OPTARG}" ;;
	esac
done

TMPDIR="$(mktemp -d -t --tmpdir=$PWD remapping_tmp_dir_XXX)"
cd "$TMPDIR"
echo "-----------------------1) Flanking sequence generation-----------------------"
# Generate chrom sizes:
samtools faidx ../"$oldgenome"
cut -f1,2 ../"$oldgenome".fai > "$oldgenome".chrom.sizes

# Filter SNPs only
bcftools filter -i 'TYPE="snp"' ../"$vcffile" -o snps_only.vcf

# Store header
bgzip snps_only.vcf
bcftools view --header-only snps_only.vcf.gz -o vcf_header.txt

# Convert vcf to bed:
gzip -d snps_only.vcf.gz
vcf2bed < snps_only.vcf > variants.bed
# The actual position of the variant is the second coordinate

# Generate the flanking sequence intervals
bedtools slop -i variants.bed -g "$oldgenome".chrom.sizes -b "$flankingseq" > flanking.bed

# Check if the intervals are less than the required length (2*flanking length + 1): this filters out reads that aren't 
# the correct number of bases long (origin: variant too close to the start or end of the chromosome)
readlength=$(bc <<< "scale=1; ($flankingseq*2)+1")
awk -v rlength="$readlength" '($3 - $2 == rlength){print $0}' flanking.bed > flanking.filtered.bed

# Get the fasta sequences for these intervals
bedtools getfasta -fi ../"$oldgenome" -bed flanking.filtered.bed -fo variants_reads.fa

# Replace the colon separators with "|":
# Storing this information for later on in the script when we split the name by "|" to extract the relevant 
# information (ALT allele, QUAL, FILT, INFO)
# This is done because the INFO column can contain ":", which means we wouldn't be able to split by ":", so "|" was 
# chosen
sed -i 's/:/|/' variants_reads.fa

# Store ref bases
awk '{print $6}' flanking.filtered.bed > old_ref_bases.txt

# Store rsIDs
awk '{print $4}' flanking.filtered.bed > rsIDs.txt

# Store variant bases
awk '{print $7}' flanking.filtered.bed > variant_bases.txt

# Store the other vcf columns
awk '{print $5, $8, $9}' flanking.filtered.bed > qual_filt_info.txt

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
sed 's/\t/|/g; s/ /|/g; s/\(.*\)|/\1 /' temp.txt | tr ' ' '\n' > variant_reads.out.fa


echo "---------------------------2) Mapping with bowtie2---------------------------"
# "index" is the prefix of the output filenames, not a directory
echo "#### Building index: may take a while ####"

# Indexes are kept in the working directory so you can keep them if you want, as they take a while to generate
cd -
bowtie2-build "$newgenome" index
echo "#### Mapping results: ####"

# Bowtie2 arguments:
# -k 1 means that bowtie2 will only report the best alignment for each variant (default behaviour)
# -k 2 means that bowtie2 will report the best alignment and the second best
# -f specifies that the reads are in fasta format
# -x [index] is the prefix of the index files

# Samtools view arguments:
# -b: output in BAM format
# -S: used to specifiy input is SAM format, but only in older versions of samtools

bowtie2 -k 2 -f -x index "$TMPDIR"/variant_reads.out.fa | samtools view -bS - > "$TMPDIR"/reads_aligned.bam

# Example output:
# >1094 reads; of these:
# >  1094 (100.00%) were unpaired; of these:
# >   7 (0.64%) aligned 0 times
# >    48 (4.39%) aligned exactly 1 time
# >    1039 (94.97%) aligned >1 times
# >99.36% overall alignment rate


echo '------------------------------3) Data extraction-----------------------------'
# Extract chromosome, rsID, position on new genome, original variant allele, qual, filter and info fields, also 
# creates old_ref_alleles.txt
cd "$TMPDIR"
samtools sort reads_aligned.bam -o reads_aligned.sorted.bam
samtools index reads_aligned.sorted.bam
../reverse_strand.py -i reads_aligned.sorted.bam -p old_ref_alleles.txt -o variants_remapped.vcf -f "$flankingseq"

# Add interval for variant position in bed format (required by getfasta)
# Reprints all the columns, adding an extra column before the pos column as the pos-1, as bedtools getfasta requires a 
# bed file
# Input: 
# [chr]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
# Output:
# [chr]	[pos-1]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' variants_remapped.vcf > variants_remapped.bed

# Getfasta to get the REF genome alleles and then extract the genome alleles without the positions
bedtools getfasta -tab -fi ../"$newgenome" -bed variants_remapped.bed | awk '{print $2}'> genome_alleles.txt
# Make all alleles caps (sometimes there are lower-case bases)
tr a-z A-Z < genome_alleles.txt > genome_alleles.fixed.txt

# Paste all the information into the correct columns
cut -f1-3 variants_remapped.vcf > temp_first_3_columns.txt
cut -f4-7 variants_remapped.vcf > temp_last_4_columns.txt
paste temp_first_3_columns.txt genome_alleles.fixed.txt temp_last_4_columns.txt > var_pre_final.vcf

# Header:
# Create list of contigs/chromosomes to be added to the header
awk '{print $1}' var_pre_final.vcf | sort | uniq > contig_names.txt
while read CHR; do echo "##contig=<ID=$CHR>"; done < contig_names.txt > contigs.txt

# Add the reference assembly
echo "##reference=$newgenomeaccession" >> contigs.txt

# Copy everything that isn't #CHROM (the column titles), ##contig or ##reference from the old header to a temp
awk '($1 !~ /^##contig/ && $1 !~ /^##reference/ && $1 !~ /^#CHROM/) {print $0}' vcf_header.txt > temp_header.txt

# Add the two headers together and add the column names
cat temp_header.txt contigs.txt > final_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> final_header.txt

# Add header to the vcf file:
cat final_header.txt | cat - var_pre_final.vcf > temp && mv temp vcf_out_with_header.vcf
# Also add it to the ref alleles file so that the line numbers match:
cat final_header.txt | cat - old_ref_alleles.txt > temp && mv temp old_ref_alleles_with_header.txt

# Test each variant to see if REF = ALT, if so, replace the current REF with the corresponding old REF (this is to 
# deal with variants such as G > G)
../replace_refs.py -i vcf_out_with_header.vcf -r old_ref_alleles_with_header.txt -o pre_final_vcf.vcf

# Run bcftools norm to swap the REF and ALT alleles if the REF doesn't match the new assembly
bgzip pre_final_vcf.vcf
bcftools norm -c ws -f ../"$newgenome" -N pre_final_vcf.vcf.gz -o ../"$outfile" -O v

cd -
bgzip "$outfile"
bcftools stats "$outfile".gz > stats.txt
gzip -d "$outfile".gz

cd "$TMPDIR"
echo "-----------------------------------Results-----------------------------------"
echo "Total number of input variants:"
grep "^[^#]" ../"$vcffile" | wc -l
echo "Total number of variants processed (after filters):"
total=$(grep "^[^>]" variant_reads.out.fa | wc -l)
echo $total
echo "Number of remapped variants:"
remapped=$(grep "^[^#]" ../"$outfile" | wc -l)
echo $remapped
echo "Percentage of remapped variants:"
percentage=$(bc <<< "scale=3;($remapped/$total)*100")
echo $percentage"%"

cd -
rm -R "$TMPDIR"