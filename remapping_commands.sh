#!/bin/bash
#requires vcf2bed, bowtie2, samtools, bedtools and reverse_strand.py (uses bamnostic)

#This basic pipeline allows you to input two assemblies and a variant vcf file, 
#and remap the variants to the new genome


#example command:
#time ./remapping_commands.sh -g droso_dm3.fasta -n GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna -a GCA_000001215.4 -v droso_variants_renamed.vcf -o test.vcf


#To be added:
# -support for indels
# -find a way to extract new assembly accession so it's not required as input? or just use new assembly file name?
# -decide what to do in the case of new REF allele being the same as the variant allele (was A > G on old assembly, now is G > G on new assembly)

#Completed:
# -filter only SNPs (not indels) using bcftools filter
# -copy filter etc columns from original vcf
# -copy vcf header to the new vcf file, but regenerate contigs from the list of contigs in new assembly, and change reference tag
# -take into account reverse strand
# -create a command line version where you input file names
# -% of remapped variants calculation


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


#1) Flanking sequence generation
echo "-----------------------1) Flanking sequence generation-----------------------"
#generate chrom sizes:
samtools faidx $oldgenome 
cut -f1,2 "$oldgenome".fai > "$oldgenome".chrom.sizes

#filter SNPs only
bcftools filter -i 'TYPE="snp"' $vcffile > snps_only.vcf

#store header
awk '$1 ~ /^#/ {print $0}' snps_only.vcf > vcf_header.txt

#convert vcf to bed:
vcf2bed < snps_only.vcf > variants.bed
#the actual position of the variant is the second coordinate
rm snps_only.vcf

#generate the flanking sequence intervals
bedtools flank -i variants.bed -g "$oldgenome".chrom.sizes -b 50 > flanking.bed
rm variants.bed "$oldgenome".chrom.sizes

#this outputs a flanking.bed file: for 1 variant we have the following lines:
#[chr] [-50 pos] [variant pos-1] [rsID] [qual] [REF allele] [ALT allele] [filter] [info]
#[chr] [variant pos] [+50 pos] [rsID] [qual] [REF allele] [ALT allele] [filter] [info]

#we want 1 line per variant, with 1 interval to create the read:
#[chr] [-50 pos] [+50 pos] [rsID] [qual] [REF allele] [ALT allele] [filter] [info]

#fix rows: put both flanking sequences on 1 row so each read is 101 bases: [50 bases]-[variant]-[50 bases]
awk 'NR == 1 {prev = $2} NR>1 && NR%2==0 {print $1"\t"prev"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} NR%2 == 1 { prev = $2}' flanking.bed > flanking.corrected.bed
rm flanking.bed

#get the fasta sequences for these intervals
bedtools getfasta -fi $oldgenome -bed flanking.corrected.bed -fo variants_reads.fa

#replace the colon separators with "|"
sed -i 's/:/|/' variants_reads.fa

#store rsIDs
awk '{print $4}' flanking.corrected.bed > rsIDs.txt

#store variant bases
awk '{print $7}' flanking.corrected.bed > variant_bases.txt

#store the other vcf columns
awk '{print $5, $8, $9}' flanking.corrected.bed > qual_filt_info.txt

#paste the names, variant bases, then fasta sequences into a new file
paste <(grep '^>' variants_reads.fa) variant_bases.txt rsIDs.txt qual_filt_info.txt <(grep -v '^>' variants_reads.fa) > temp.txt
rm flanking.corrected.bed variant_bases.txt rsIDs.txt qual_filt_info.txt  variants_reads.fa

#reformat the fasta ID
sed 's/\t/|/g; s/ /|/g; s/\(.*\)|/\1 /' temp.txt | tr ' ' '\n' > variant_reads.out.fa
rm temp.txt






#2) Mapping with bowtie2
echo "---------------------------2) Mapping with bowtie2---------------------------"
#"index" is the prefix of the output filenames, not a directory
echo "#### Building index: may take a while ####"
bowtie2-build $newgenome index
echo "#### Mapping results: ####"
bowtie2 -f -x index variant_reads.out.fa | samtools view -bS - > reads_aligned.bam
#example output:
#>1094 reads; of these:
#>  1094 (100.00%) were unpaired; of these:
#>   7 (0.64%) aligned 0 times
#>    48 (4.39%) aligned exactly 1 time
#>    1039 (94.97%) aligned >1 times
#>99.36% overall alignment rate





#3) Data extraction
echo '------------------------------3) Data extraction-----------------------------'
#extract chromosome, rsID, position on new genome, original variant allele, qual, filter and info fields
samtools sort reads_aligned.bam -o reads_aligned.sorted.bam
samtools index reads_aligned.sorted.bam
./reverse_strand.py -i reads_aligned.sorted.bam -o variants_remapped.vcf
rm reads_aligned.bam reads_aligned.sorted.bam reads_aligned.sorted.bam.bai

#add interval for variant position in bed format (required by getfasta)
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' variants_remapped.vcf > variants_remapped.bed

#getfasta to get the REF genome alleles and then extract the genome alleles without the positions
bedtools getfasta -tab -fi $newgenome -bed variants_remapped.bed | awk '{print $2}'> genome_alleles.txt
rm variants_remapped.bed
#make all alleles caps (sometimes there are lower-case bases)
tr a-z A-Z < genome_alleles.txt > genome_alleles.fixed.txt

#paste all the information into the correct columns
awk '{print $1"\t"$2"\t"$3}' variants_remapped.vcf > temp_first_3_columns.txt
awk '{print $4"\t"$5"\t"$6"\t"$7}' variants_remapped.vcf > temp_last_4_columns.txt
paste temp_first_3_columns.txt genome_alleles.fixed.txt temp_last_4_columns.txt > var_pre_final.vcf
rm temp_first_3_columns.txt temp_last_4_columns.txt genome_alleles.txt variants_remapped.vcf genome_alleles.fixed.txt


#header:
#create list of contigs/chromosomes to be added to the header
awk '{print $1}' var_pre_final.vcf | sort | uniq > contig_names.txt
while read CHR; do echo "##contig=<ID=$CHR>"; done < contig_names.txt > contigs.txt
rm contig_names.txt

#add the reference assembly
echo "##reference="$newgenomeaccession >> contigs.txt

#copy everything that isn't #CHROM (the column titles), ##contig or ##reference from the old header to a temp
awk '($1 !~ /^##contig/ && $1 !~ /^##reference/ && $1 !~ /^#CHROM/) {print $0}' vcf_header.txt > temp_header.txt
rm vcf_header.txt

#add the two headers together and add the column names
cat temp_header.txt contigs.txt > final_header.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> final_header.txt
rm temp_header.txt contigs.txt

#add header to the vcf file:
cat final_header.txt | cat - var_pre_final.vcf > temp && mv temp $outfile
rm final_header.txt var_pre_final.vcf


echo "-----------------------------------Results-----------------------------------"
echo "Total number of input variants:"
grep "^[^#]" $vcffile | wc -l
echo "Total number of variants processed (after filters):"
total=$(grep "^[^>]" variant_reads.out.fa | wc -l)
rm variant_reads.out.fa
echo $total
echo "Number of remapped variants:"
remapped=$(grep "^[^#]" $outfile | wc -l)
echo $remapped
echo "Percentage of remapped variants:"
percentage=$(bc <<< "scale=3;($remapped/$total)*100")
echo $percentage"%"