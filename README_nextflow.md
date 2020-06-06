# Variant remapping using Nextflow 

## Prerequisites
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

## Installation
```bash
git clone https://github.com/EBIvariation/variant-remapping.git
conda env create -f conda.yml
```

## Execution 
First activate the Conda environment:
```bash
conda activate variant-remapping
pip install -r requirements.txt
```

Run the test script to check that you have all the right dependencies installed properly 
```bash
tests/test_pipeline.sh
```

Or with real data:
```bash
nextflow run remapping_pipeline.nf \
  --oldgenome <genome.fa> \
  --newgenome <new_genome.fa> \
  --vcffile <source.vcf> \
  --outfile <remap.vcf>
```
