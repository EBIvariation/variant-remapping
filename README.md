# variant-remapping
Pipeline for remapping VCF variants between two arbitrary assemblies in FASTA format. No chain file is required.

**Method**: creates reads from the flanking sequences of each variant, then maps them to the new assembly using bowtie2.

Currently, it only supports SNPs, and not indels.

Prerequisites:
- reverse_strand.py, included (reverse strand allele correction, uses bamnostic)
- write_header.py, included (necessary for VCF header generation)
- vcf2bed
- bowtie2
- samtools
- bedtools
- bcftools

## Input
- Old genome assembly file (FASTA format): the genome you have variants for.
- New genome assembly file (FASTA format): the genome you want to remap the variants to.
- Variants file (VCF format): contains the list of variants you want to remap.
- Accession name for the new assembly: required to recreate the header for the output VCF file.

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
`chmod a+x remapping_commands.sh reverse_strand.py write_header.py`

Command usage: (see Input for details)
`./remapping_commands.sh -g [oldgenome] -n [newgenome] -a [newgenomeaccession] -v [vcffile] -o [outfile name]`

You can use `time` at the beginning of the previous command to get the runtime and CPU time.

**Example command:**
`time ./remapping_commands.sh -g droso_dm3.fasta -n droso_dm6.fasta -a GCA_000001215.4 -v droso_variants_renamed.vcf -o test.vcf`
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

## Task list
- [x] Filter only SNPs
- [x] Copy VCF header and recreate necessary lines
- [x] Take into account reverse strand mapping (the REF alleles change)
- [x] Create a command line version
- [x] Calculate % of remapped variants
- [ ] Support for indels
- [ ] Find a way to extract new assembly accession so it's not required as input?
