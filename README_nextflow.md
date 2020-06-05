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

Test data is located in the `resources` directory:
```bash
nextflow run genome_preparation.nf
```

Or with real data:
```bash
nextflow run genome_preparation.nf --fasta <genome.fa> --outdir <output directory>
```
