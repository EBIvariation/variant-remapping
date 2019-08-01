#!/bin/bash
#requires vcf2bed, bowtie2, samtools, bedtools and reverse_strand.py (uses bamnostic)

#This basic pipeline allows you to input two assemblies and a variant vcf file, 
#and remap the variants to the new genome


#example command:
#time ./remapping_commands.sh -g droso_dm3.fasta -n GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna -a GCA_000001215.4 -v droso_variants_renamed.vcf -o test.vcf


oldgenome=''
newgenome=''
newgenomeaccession=''
vcffile=''
outfile=''

while getopts 'g:n:a:v:o:' flag; do
	case "${flag}" in
		g) oldgenome="${OPTARG}";;
		n) newgenome="${OPTARG}" ;;
		a) newgenomeaccession="${OPTARG}" ;;
		v) vcffile="${OPTARG}" ;;
		o) outfile="${OPTARG}" ;;
	esac
done

TMPDIR="$(mktemp -d -t --tmpdir=$PWD remapping_tmp_dir_XXX)"
cd "$TMPDIR"
echo "-----------------------1) Flanking sequence generation-----------------------"
#generate chrom sizes:
samtools faidx ../"$oldgenome"
cut -f1,2 ../"$oldgenome".fai > "$oldgenome".chrom.sizes

#filter SNPs only
bcftools filter -i 'TYPE="snp"' ../"$vcffile" -o snps_only.vcf

#store header
bgzip snps_only.vcf
bcftools view --header-only snps_only.vcf.gz -o vcf_header.txt

#convert vcf to bed:
gzip -d snps_only.vcf.gz
vcf2bed < snps_only.vcf > variants.bed
#the actual position of the variant is the second coordinate

#generate the flanking sequence intervals
bedtools slop -i variants.bed -g "$oldgenome".chrom.sizes -b 50 > flanking.bed


#get the fasta sequences for these intervals
bedtools getfasta -fi ../"$oldgenome" -bed flanking.bed -fo variants_reads.fa

#replace the colon separators with "|":
#storing this information for later on in the script when we split the name by "|" to extract the relevant information (ALT allele, QUAL, FILT, INFO)
#this is done because the INFO column can contain ":", which means we wouldn't be able to split by ":", so "|" was chosen
sed -i 's/:/|/' variants_reads.fa

#store ref bases
awk '{print $6}' flanking.bed > old_ref_bases.txt

#store rsIDs
awk '{print $4}' flanking.bed > rsIDs.txt

#store variant bases
awk '{print $7}' flanking.bed > variant_bases.txt

#store the other vcf columns
awk '{print $5, $8, $9}' flanking.bed > qual_filt_info.txt

#paste the names, variant bases, then fasta sequences into a new file
paste <(grep '^>' variants_reads.fa) old_ref_bases.txt variant_bases.txt rsIDs.txt qual_filt_info.txt <(grep -v '^>' variants_reads.fa) > temp.txt

#reformat the fasta ID: inconsistencies in the separators, and no new line before the sequence mean that this next command is a bit ugly
#Input:
#>[chr]|[pos interval]   [REF]       [ALT]       [rsID]    [QUAL] [FILT] [INFO]      [seq]
#replace all spaces and tabs with "|":
#>[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO]|[seq]
#replace the last "|" with a space (this is between the header and the sequence):
#>[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO] [seq]
#and finally replace the space with a newline:
#Output:
#>[chr]|[pos interval]|[REF]|[ALT|[rsID]|[QUAL|[FILT]|[INFO]
#[seq]
sed 's/\t/|/g; s/ /|/g; s/\(.*\)|/\1 /' temp.txt | tr ' ' '\n' > variant_reads.out.fa


echo "---------------------------2) Mapping with bowtie2---------------------------"
#"index" is the prefix of the output filenames, not a directory
echo "#### Building index: may take a while ####"
#indexes are kept in the working directory so you can keep them if you want, as they take a while to generate
cd ../
bowtie2-build "$newgenome" index
echo "#### Mapping results: ####"
#bowtie2 arguments:
#-k 1 means that bowtie2 will only report the best alignment for each variant (default behaviour)
#-f specifies that the reads are in fasta format
#-x [index] is the prefix of the index files
#samtools view arguments:
#-b: output in BAM format
#-S: used to specifiy input is SAM format, but only in older versions of samtools
bowtie2 -k 1 -f -x index "$TMPDIR"/variant_reads.out.fa | samtools view -bS - > "$TMPDIR"/reads_aligned.bam
#example output:
#>1094 reads; of these:
#>  1094 (100.00%) were unpaired; of these:
#>   7 (0.64%) aligned 0 times
#>    48 (4.39%) aligned exactly 1 time
#>    1039 (94.97%) aligned >1 times
#>99.36% overall alignment rate


echo '------------------------------3) Data extraction-----------------------------'
#extract chromosome, rsID, position on new genome, original variant allele, qual, filter and info fields, also creates old_ref_alleles.txt
cd "$TMPDIR"
samtools sort reads_aligned.bam -o reads_aligned.sorted.bam
samtools index reads_aligned.sorted.bam
.././reverse_strand.py -i reads_aligned.sorted.bam -p old_ref_alleles.txt -o variants_remapped.vcf

#add interval for variant position in bed format (required by getfasta)
#reprints all the columns, adding an extra column before the pos column as the pos-1, as bedtools getfasta requires a bed file
#Input: 
#[chr]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
#Output:
#[chr]	[pos-1]	[pos]	[rsID]	[ALT]	[QUAL]	[FILT]	[INFO]
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' variants_remapped.vcf > variants_remapped.bed

#getfasta to get the REF genome alleles and then extract the genome alleles without the positions
bedtools getfasta -tab -fi ../"$newgenome" -bed variants_remapped.bed | awk '{print $2}'> genome_alleles.txt
#make all alleles caps (sometimes there are lower-case bases)
tr a-z A-Z < genome_alleles.txt > genome_alleles.fixed.txt

#paste all the information into the correct columns
cut -f1-3 variants_remapped.vcf > temp_first_3_columns.txt
cut -f4-7 variants_remapped.vcf > temp_last_4_columns.txt
paste temp_first_3_columns.txt genome_alleles.fixed.txt temp_last_4_columns.txt > var_pre_final.vcf

#header:
#create list of contigs/chromosomes to be added to the header
awk '{print $1}' var_pre_final.vcf | sort | uniq > contig_names.txt
while read CHR; do echo "##contig=<ID=$CHR>"; done < contig_names.txt > contigs.txt

#add the reference assembly
echo "##reference=$newgenomeaccession" >> contigs.txt

#copy everything that isn't #CHROM (the column titles), ##contig or ##reference from the old header to a temp
awk '($1 !~ /^##contig/ && $1 !~ /^##reference/ && $1 !~ /^#CHROM/) {print $0}' vcf_header.txt > temp_header.txt

#add the two headers together and add the column names
cat temp_header.txt contigs.txt > final_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> final_header.txt

#add header to the vcf file:
cat final_header.txt | cat - var_pre_final.vcf > temp && mv temp vcf_out_with_header.vcf
#also add it to the ref alleles file so that the line numbers match:
cat final_header.txt | cat - old_ref_alleles.txt > temp && mv temp old_ref_alleles_with_header.txt

#test each variant to see if REF = ALT, if so, replace the current REF with the corresponding old REF (this is to deal with variants like G > G)
.././replace_refs.py -i vcf_out_with_header.vcf -r old_ref_alleles_with_header.txt -o pre_final_vcf.vcf

#run bcftools norm to swap the REF and ALT alleles if the REF doesn't match the new assembly
bgzip pre_final_vcf.vcf
bcftools norm -c ws -f ../"$newgenome" -N pre_final_vcf.vcf.gz -o ../"$outfile" -O v

cd ../
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

cd ../
rm -R "$TMPDIR"