# variant-remapping
Pipeline for remapping VCF variants between two arbitrary assemblies in FASTA format. No chain file is required.

**Method**: creates reads from the flanking sequences of each variant, then maps them to the new assembly using 
bowtie2.

Currently, it only supports SNPs, and not indels.

**Prerequisites:**
- `reverse_strand.py`, included (uses [pysam](https://pysam.readthedocs.io/en/latest/api.html), allows reverse strand 
allele correction (see Note below for further explanation), and filtering based on bowtie2's alignment score)
- `replace_refs.py`, included (uses old REF allele in the case of new REF=ALT, these then get swapped (see Note 2))
- [vcf2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/vcf2bed.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [bcftools](http://www.htslib.org/download/)

**Note: reverse strand correction:**  
When a variant is mapped to the reverse strand, the corresponding allele in the output VCF file is reversed, aka 
converted to the forward strand, as alleles in VCF files are always described on the forward strand. For example:
INPUT VCF:
Old genome: `G (REF) > A (ALT)`
This maps onto the new genome on the reverse strand:
OUTPUT VCF:
New genome: `C (REF) > T (ALT)`

**Note 2: what happens when the new REF=ALT?**
Example: Original variant on the old reference sequence: `G (REF) > T (ALT)` 
Once remapped to the new reference sequence, we get: `T (REF) > T (ALT)`  
What happens is `replace_refs.py` will search for the old REF allele, and replace the new REF with it. So we get:  
`G (REF) > T (ALT)`  
And then `bcftools norm` will check to see if any REF alleles are different to the new reference genome, and if so, it 
will swap them:
`T (REF) > G (ALT)`  
This makes sense because we know that T is now the REF allele, which means the original G allele was technically a 
variant of this new reference genome.  

## Input
- Old genome assembly file (FASTA format): the genome you have variants for.
- New genome assembly file (FASTA format): the genome you want to remap the variants to.
- Variants file (VCF format): contains the list of variants you want to remap.
- The length of the flanking sequences that generate the reads.
- Percentage of the flanking sequences that should be used as Alignment Score cut-off threshold.

## Output
A VCF file containing:
- remapped coordinates (position on the new assembly)
- rsIDs
- the correct chromosome/contig names
- the new REF alleles (reverse strand mapping taken into account)
- the ALT, QUAL, FILT and INFO columns of the input VCF

Console output: 
```
"Total number of input variants:"
[number]
"Total number of variants processed (after filters):"
[number]
"Number of remapped variants:"
[number]
"Percentage of remapped variants:"
[percentage]
```

## Usage
First, make sure the scripts are executable:
`chmod a+x remapping_commands.sh reverse_strand.py replace_refs.py`

Command usage: (see Input for details)
```
./remapping_commands.sh \
  -g [old genome] \
  -n [new genome] \
  -v [vcf file] \
  -f [flanking sequence length] \
  -s [score percentage cut-off] \
  -o [outfile name]
```

You can use `time` at the beginning of the previous command to get the runtime and CPU time.

**Example:**
"I want to remap the variants in `droso_variants_renamed.vcf` from `droso_dm3.fasta` to `droso_dm6.fasta` (its 
accession is: `GCA_000001215.4`), with flanking sequences of 50 bases, which will create 101-base reads. The alignment score cut-off will be -(50 x 0.6) = -30, meaning that reads with alignment scores lower than -30 will not be kept. The remapped variants will be in `test.vcf`."  

**Command:**
```
time ./remapping_commands.sh \
  -g droso_dm3.fasta \
  -n droso_dm6.fasta \
  -v droso_variants_renamed.vcf \
  -f 50 \
  -s 0.6 \
  -o test.vcf
```
**Output:**
```
-----------------------1) Flanking sequence generation-----------------------
---------------------------2) Mapping with bowtie2---------------------------
Mapping results:
1077 reads; of these:
  1077 (100.00%) were unpaired; of these:
    7 (0.65%) aligned 0 times
    48 (4.46%) aligned exactly 1 time
    1022 (94.89%) aligned >1 times
99.35% overall alignment rate
------------------------------3) Data extraction-----------------------------
-----------------------------------Results-----------------------------------
Total number of input variants:
1094
Total number of variants processed (after filters):
1077
Number of remapped variants:
1070
Percentage of remapped variants:
99.300%

real	0m1.363s
user	0m1.091s
sys	0m0.171s
```