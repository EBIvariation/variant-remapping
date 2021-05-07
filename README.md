# variant-remapping
Pipeline for remapping VCF variants between two arbitrary assemblies in FASTA format. No chain file is required. 
It does however assume that the source and destination genomes are closely related 
and was designed with the explicit purpose of lifting over variants from one ersiong of the genome to another.

**Method**: creates reads from the flanking sequences of each variant, then maps them to the new assembly using 
minimap2.

Currently, it only SNPs and short indels but has not been tested with larger or more complex variants.

## Prerequisites
To run this pipeline you will need to install and configure [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) version 20.7 or later. 
The pipeline uses other software that needs to be downloaded and installed locally. You can obtain them manually or use [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

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
    --oldgenome <genome.fa> \
    --newgenome <new_genome.fa> \
    --vcffile <source.vcf> \
    --outfile <remap.vcf> \
```

### Input
- `--oldgenome`: Old genome assembly file (FASTA format): the genome you have variants for.
- `--newgenome`: New genome assembly file (FASTA format): the genome you want to remap the variants to.
- `--vcffile`: Variants file (VCF format): contains the list of variants you want to remap.

### Output
`--outfile` specify a VCF file containing:
- remapped coordinates (chromosome and position on the new assembly)
- the new REF alleles from the new assembly
- the ALT field possibly modified if the strand or REF has changed ID, QUAL, FILT and INFO columns of the input VCF
- Additional fields in the INFO column
- FORMAT and Sample columns if they were present in the input
