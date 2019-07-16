#requires vcf2bed, bowtie2, samtools, bedtools

#This basic pipeline allows you to input two assemblies and a variant vcf file, 
#and remap the variants to the new genome

#It's in a very basic state:
#To be added:
# -filter only SNPs (not indels) using bcftools filter
# -copy vcf header to the new vcf file, but regenerate contigs from the list of contigs in new assembly, and change reference tag
# -copy filter etc columns from original vcf



#environment variables:
#input vcf file:
export vcfFile="droso_variants_renamed.vcf"
#input old genome for which you have variants:
export oldgenome="droso_dm3.fasta"
#input new genome that you want to remap the variants to:
export newgenome="droso_dm6.fasta"
#species name for file naming:
export species="droso"





#1) Flanking sequence generation

#generate chrom sizes:
samtools faidx $oldgenome 
cut -f1,2 "$oldgenome".fai > "$species".chrom.sizes

#convert vcf to bed:
vcf2bed < $vcfFile > "$species"_variants.bed
#the actual position of the variant is the second coordinate

#generate the flanking sequence intervals
bedtools flank -i "$species"_variants.bed -g "$species".chrom.sizes -b 50 > flanking.bed

#fix rows: put both flanking sequences on 1 row = 1 read and add 1 to start position so the read is 100 bases:
awk 'NR == 1 {prev = $2} NR>1 && NR%2==0 {print $1"\t"prev+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} NR%2 == 1 { prev = $2}' flanking.bed > flanking.corrected.bed

#get the fasta sequences for these intervals
bedtools getfasta -name -fi $oldgenome -bed flanking.corrected.bed -fo "$species"_variants_reads.fa

#store variant bases
awk '{print $7}' flanking.corrected.bed > variant_bases.txt

paste <(grep '^>' "$species"_variants_reads.fa) variant_bases.txt <(grep '^[A|C|G|T]' droso_variants_reads.fa) > temp.txt

#store the variant allele in the fasta ID
sed $'s/\t/:/' temp.txt | tr "\t" "\n" > "$species"_variant_reads.out.fa

#remove intermediate files:
rm temp.txt variant_bases.txt flanking.corrected.bed "$species"_variants_reads.fa flanking.bed "$species"_variants.bed "$species".chrom.sizes "$oldgenome".fai





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

#extract chromosome, rsID, position on new genome, and original variant allele
samtools view "$species"_reads_aligned.bam | awk '($1 ~ ":") {split($1,info,":"); print $3"\t"$4+50"\t"info[1]"\t"info[5]}' > variants_remapped.vcf

#delete missing chromosomes (*)
awk '($1 != "*") {print $0}' variants_remapped.vcf > variants_remapped.fixed.vcf

#add interval for variant position in bed format (required by getfasta)
awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' variants_remapped.fixed.vcf > variants_remapped.bed

#getfasta to get the genome alleles
bedtools getfasta -tab -fi droso_dm6.fasta -bed variants_remapped.bed -fo genome_allele_pos.txt

#extract the genome alleles without the positions
awk '{print $2}' genome_allele_pos.txt > genome_alleles.txt

#paste all the information into a vcf format
paste variants_remapped.fixed.vcf genome_alleles.txt <(cut -f 4 variants_remapped.fixed.vcf) > var_pre_final.vcf

#remove the duplicate column
cut -f4 --complement var_pre_final.vcf > "$species"_variants_remapped.vcf

#delete intermediate files:
rm var_pre_final.vcf variants_remapped.fixed.vcf genome_alleles.txt genome_allele_pos.txt variants_remapped.bed variants_remapped.fixed.vcf variants_remapped.vcf "$species"_variant_reads.out.fa




#$10 = sequence
#$2 = flag

#test:
#rs80391166
#TGCTGGTGAAGAATATATATATTCTGTATAGGTAAAGGTTCGGCAACGCCTTCTTCGAACTCCCACAGACTTTTTAAAGAATTTGGTATACTCAATTACT
#echo ${x:49:1}
#>T

#To do:

