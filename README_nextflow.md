# Variant-remapping using nextflow 

#### Prerequisites:
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)


#### Installation
```bash
git clone https://github.com/EBIvariation/variant-remapping.git
conda env create -f conda.yml
```

#### Execution 
First activate the conda environment
```bash
conda activate variant-remapping
pip install -r requirements.txt
```

Test data located in `resources` directory
```bash
nextflow run genome_preparation.nf
```

Or with real data
```bash
nextflow run genome_preparation.nf --fasta <genome.fa> --outdir <output directory>
```
