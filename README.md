# variant-remapping
Pipeline for remapping VCF variants between two arbitrary assemblies in FASTA format. No chain file is required.

**Method**: creates reads from the flanking sequences of each variant, then maps them to the new assembly using 
bowtie2.

Currently, it only SNPs and short indels but has not been tested with larger or more complex variants

## Prerequisites
To run this pipeline you will need to install and configure [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation). 
The pipeline uses other software that needs to be downloaded and installed locally. You can obtain them manually or use [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### Installation using conda
```bash
git clone https://github.com/EBIvariation/variant-remapping.git
conda env create -f conda.yml
conda activate variant-remapping
pip install -r requirements.txt
```

### Installation without conda
Download, manually install the following program and make sure the executable are in your `PATH`
- [vcf2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/vcf2bed.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [bcftools](http://www.htslib.org/download/)
- [Python](https://www.python.org/downloads/)

Then run
```bash
git clone https://github.com/EBIvariation/variant-remapping.git
pip install -r requirements.txt
```

## Testing the installation
Run the test script to check that you have all the right dependencies installed properly 
```bash
tests/test_pipeline.sh
```

## Executing the pipeline
```bash
nextflow run main.nf 
    --oldgenome <genome.fa> 
    --newgenome <new_genome.fa> 
    --vcffile <source.vcf> 
    --outfile <remap.vcf>
    [--flankingseq 50]
    [--scorecutoff 0.6]
    [--diffcutoff 0.04]
```

### Input
- `--oldgenome`: Old genome assembly file (FASTA format): the genome you have variants for.
- `--newgenome`: New genome assembly file (FASTA format): the genome you want to remap the variants to.
- `--vcffile`: Variants file (VCF format): contains the list of variants you want to remap.
- `--flankingseq`: The length of the flanking sequences that generate the reads.
- `--scorecutoff`: Percentage of the flanking sequences that should be used as Alignment Score cut-off threshold.
- `--diffcutoff`: Percentage of the flanking sequences that should be used as AS-XS difference cut-off threshold.

### Output
`--outfile` specify a VCF file containing:
- remapped coordinates (position on the new assembly)
- rsIDs
- the correct chromosome/contig names
- the new REF alleles (reverse strand mapping taken into account)
- the ALT, QUAL, FILT and INFO columns of the input VCF

### Example:
"I want to remap the variants in `droso_variants_renamed.vcf` from `droso_dm3.fasta` to `droso_dm6.fasta` (its 
accession is: `GCA_000001215.4`), with flanking sequences of 50 bases, which will create 101-base reads. The alignment 
score cut-off will be -(50 x 0.6) = -30, meaning that reads with alignment scores lower than -30 will not be kept. The 
AS-XS difference threshold will be 0.04, meaning that AS-XS with a difference of less than 50 * 0.04 = 2 will not be 
kept. The remapped variants will be in `test.vcf`."  
