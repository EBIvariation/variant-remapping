#requires vcf2bed, bowtie2, samtools, bedtools, and write_header.py

#This basic pipeline allows you to input two assemblies and a variant vcf file, 
#and remap the variants to the new genome

#It's in a very basic state:
#To be added:
# -take into account reverse strand
# -create a command line version where you input file names

#Completed:
# -filter only SNPs (not indels) using bcftools filter
# -copy filter etc columns from original vcf
# -copy vcf header to the new vcf file, but regenerate contigs from the list of contigs in new assembly, and change reference tag




#environment variables:
#input vcf file:
export vcfFile="droso_variants_renamed.vcf"
#input old genome for which you have variants:
export oldgenome="droso_dm3.fasta"
#input new genome that you want to remap the variants to:
export newgenome="droso_dm6.fasta"
#assembly accession for the new genome:
export newgenomeaccession="GCA_000001215.4"
#species name for file naming:
export species="droso"





#1) Flanking sequence generation

#generate chrom sizes:
samtools faidx $oldgenome 
cut -f1,2 "$oldgenome".fai > "$species".chrom.sizes

#filter SNPs only
bcftools filter -i 'TYPE="snp"' $vcfFile > snps_only.vcf

#store header
awk '$1 ~ /^#/ {print $0}' snps_only.vcf > vcf_header.txt

#convert vcf to bed:
vcf2bed < snps_only.vcf > "$species"_variants.bed
#the actual position of the variant is the second coordinate
rm snps_only.vcf

#generate the flanking sequence intervals
bedtools flank -i "$species"_variants.bed -g "$species".chrom.sizes -b 50 > flanking.bed
rm "$species"_variants.bed "$species".chrom.sizes

#fix rows: put both flanking sequences on 1 row = 1 read and add 1 to start position so the read is 100 bases:
awk 'NR == 1 {prev = $2} NR>1 && NR%2==0 {print $1"\t"prev+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} NR%2 == 1 { prev = $2}' flanking.bed > flanking.corrected.bed
rm flanking.bed

#get the fasta sequences for these intervals
bedtools getfasta -fi $oldgenome -bed flanking.corrected.bed -fo "$species"_variants_reads.fa

#replace the colon separators with "|"
sed 's/:/|/' "$species"_variants_reads.fa > "$species"_variants_reads2.fa

#store rsIDs
awk '{print $4}' flanking.corrected.bed > rsIDs.txt

#store variant bases
awk '{print $7}' flanking.corrected.bed > variant_bases.txt

#store the other vcf columns
awk '{print $5, $8, $9}' flanking.corrected.bed > qual_filt_info.txt

#paste the names, variant bases, then fasta sequences into a new file
paste <(grep '^>' "$species"_variants_reads2.fa) variant_bases.txt rsIDs.txt qual_filt_info.txt <(grep '^[A|C|G|T]' "$species"_variants_reads.fa) > temp.txt
rm flanking.corrected.bed variant_bases.txt rsIDs.txt qual_filt_info.txt "$species"_variants_reads2.fa "$species"_variants_reads.fa

#reformat the fasta ID
sed 's/\t/|/g; s/ /|/g; s/\(.*\)|/\1 /' temp.txt | tr ' ' '\n' > "$species"_variant_reads.out.fa
rm temp.txt






#2) Mapping with bowtie2

# index the target assembly
#"index" is the name of the output files, not a directory
mkdir "$newgenome"_index
bowtie2-build $newgenome "$newgenome"_index"/$newgenome"_index

bowtie2 -f -x "$newgenome"_index"/$newgenome"_index "$species"_variant_reads.out.fa | samtools view -bS - > "$species"_reads_aligned.bam
#example output:
#>1094 reads; of these:
#>  1094 (100.00%) were unpaired; of these:
#>   7 (0.64%) aligned 0 times
#>    48 (4.39%) aligned exactly 1 time
#>    1039 (94.97%) aligned >1 times
#>99.36% overall alignment rate





#3) Data extraction

#extract chromosome, rsID, position on new genome, original variant allele, qual, filter and info fields
samtools view "$species"_reads_aligned.bam | awk '($1 ~ "|") {split($1,info,"|"); print $3"\t"$4+50"\t"info[4]"\t"info[3]"\t"info[5]"\t"info[6]"\t"info[7]}' > variants_remapped.vcf
rm "$species"_reads_aligned.bam

#delete missing chromosomes (*)
awk '($1 != "*") {print $0}' variants_remapped.vcf > variants_remapped.fixed.vcf
rm variants_remapped.vcf

#add interval for variant position in bed format (required by getfasta)
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' variants_remapped.fixed.vcf > variants_remapped.bed

#getfasta to get the genome alleles and then extract the genome alleles without the positions
bedtools getfasta -tab -fi droso_dm6.fasta -bed variants_remapped.bed | awk '{print $2}'> genome_alleles.txt
rm variants_remapped.bed

#paste all the information into the correct columns
awk '{print $1"\t"$2"\t"$3}' variants_remapped.fixed.vcf > temp_var1.txt
awk '{print $4"\t"$5"\t"$6"\t"$7}' variants_remapped.fixed.vcf > temp_var2.txt
paste temp_var1.txt genome_alleles.txt temp_var2.txt > var_pre_final.vcf
rm temp_var1.txt temp_var2.txt genome_alleles.txt variants_remapped.fixed.vcf

#header:
#create list of contigs/chromosomes to be added to the header
awk '{print $1}' var_pre_final.vcf | sort | uniq > contig_names.txt
while read CHR; do echo "##contig=<ID>=$CHR"; done < contig_names.txt > contigs.txt
rm contig_names.txt

#add the reference assembly
echo "##reference="$newgenomeaccession >> contigs.txt

#copy everything that isn't ##contig or ##reference from the old header to a temp
awk '($1 !~ /^##contig/ && $1 !~ /^##reference/) {print $0}' vcf_header.txt > temp_header.txt
rm vcf_header.txt

#find the end of the INFO tags and print the contigs (using custom python script)
chmod a+x write_header.py
./write_header.py -i temp_header.txt -c contigs.txt -o final_header.txt
rm temp_header.txt contigs.txt

#add header to the vcf file:
cat final_header.txt | cat - var_pre_final.vcf > temp && mv temp "$species"_variants_remapped.vcf
rm final_header.txt var_pre_final.vcf




#bam:
#$10 = sequence
#$2 = flag

#test:
#rs80391166
#TGCTGGTGAAGAATATATATATTCTGTATAGGTAAAGGTTCGGCAACGCCTTCTTCGAACTCCCACAGACTTTTTAAAGAATTTGGTATACTCAATTACT
#echo ${x:49:1}
#>T


