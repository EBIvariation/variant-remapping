# Evaluation of the variant remapping pipeline

To assess the quality of the remapping pipeline we remap GiaB variants from GRCh37 to GRCh38 and compare them with GiaB variants from GRCh38. To ensure that the remapped and standard variants are comparable, we first filter the gold standard to keep only variants that are present in both datasets.

We then compare the remapped and standard dataset using [hap.py](https://github.com/Illumina/hap.py)

You need at least 30 GB RAM and 45 GB disk space to run this pipeline.

## Prerequisite and Environment variables

* `GIAB_DIR`: directory where the gold standard variants are stored
* `REMAPPING_DIR`: directory where the output of the remapping pipeline is stored
* `GENOME_DIR`: directory where the reference genomes are stored
* `SINGULARITY_DIR`: directory where the singularity images are stored

Some software required by the remapping pipeline are also assumed to be installed and in the `PATH`.

## Preparing GiaB datasets 

Set `GIAB_DIR` to points to the directory where the gold standard variants will be stored.

|:warning: This can be skipped if the datasets have been prepared before.|
|---|

The variant remapping pipeline was run on two variant datasets that were downloaded from the GiaB website. The variant data is here: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/

The reference genomes used for calling the GiaB variants are defined in [the FAQ point 5](https://www.nist.gov/programs-projects/faqs-genome-bottle)

GRCh37 reference with decoy:

```bash
mkdir -p ${GENOME_DIR}/hs37d5
cd ${GENOME_DIR}/hs37d5
wget --no-check-certificate https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -dk hs37d5.fa.gz
```

GRCh38 reference with no ALT loci: 

```bash
mkdir -p ${GENOME_DIR}/GRCh38_no_alt
cd ${GENOME_DIR}/GRCh38_no_alt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -dk GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Prepare directory structure:

```bash
mkdir -p ${GIAB_DIR}/NA12878/GRCh37
mkdir -p ${GIAB_DIR}/NA12878/GRCh38
```

Download the variant datasets:

```bash
cd ${GIAB_DIR}/NA12878/GRCh37
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed
```

```bash
cd ${GIAB_DIR}/NA12878/GRCh38
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed
```

|:warning: GRCh38 is not annotated -> Need to annotate it with dbsnp data|
|---|

Download dbSNPs vcf to get rsids:

```bash
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz
```

|:warning: dbSNP chromosome name are 1, 2, 3, ... where hs37d5 is chr1, chr2, chr3 ...  We'll prefix all chromosome name with chr|
|---|

```bash
zcat common_all_20180418.vcf.gz | awk '{if (/^#/){print } else{print "chr"$0}}' | bgzip -c > common_all_20180418_with_chr.vcf.gz
tabix -p vcf common_all_20180418_with_chr.vcf.gz
```

```bash 
bcftools annotate \
  -a common_all_20180418_with_chr.vcf.gz \
  -c ID HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
| bgzip -c > HG001_GRCh38_annotated.vcf.gz
```

|:warning: Hap.py compares the variants and the genotype provided. However the remapping currently does not carry over the genotypes. This means that the genotypes needs to be mocked in both the remapped and the standard datasets to allow the variants to be compares accurately|
|---|

```bash
cd ${GIAB_DIR}
gunzip -c NA12878/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
   | awk 'BEGIN{OFS="\t"} {if  (/^#/){print}else{$9="GT"; $10="1/1"; print $0}}' \
   | bgzip -c > NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes.vcf.gz

gunzip -c NA12878/GRCh38/HG001_GRCh38_annotated.vcf.gz \
   | awk 'BEGIN{OFS="\t"} {if  (/^#/){print}else{$9="GT"; $10="1/1"; print $0}}' \
   | bgzip -c > NA12878/GRCh38/HG001_GRCh38_annotated_Fake_Genotypes.vcf.gz 
```

Then we get all shared rsID from both files:

```bash
# Match the rsids that are in both GRCh37 and GRCh38
cat <(zcat NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes.vcf.gz | grep -v '^#' | cut -f 3 | sort | uniq ) \
    <(zcat NA12878/GRCh38/HG001_GRCh38_annotated_Fake_Genotypes.vcf.gz | grep -v '^#' | cut -f 3 | sort | uniq ) \
    | sort | uniq -c | awk '{if ($1>1){print $2}}' > valid_rsids.txt
```

Only keep the variants that are the same in both datasets:

```bash
cat <(zcat NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes.vcf.gz | grep  '^#') \
    <(zcat NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes.vcf.gz | grep -wFf valid_rsids.txt) \
    | bgzip -c > NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz 

cat <(zcat NA12878/GRCh38/HG001_GRCh38_annotated_Fake_Genotypes.vcf.gz | grep '^#') \
    <(zcat NA12878/GRCh38/HG001_GRCh38_annotated_Fake_Genotypes.vcf.gz | grep -wFf valid_rsids.txt) \
    | bgzip -c > NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz 
```

Then index the final files:

```bash 
tabix -p vcf NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz 
tabix -p vcf NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz 
```

Confirm that the number of variant is the same in both VCF:

```bash 
zcat NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz | grep -v '#' | cut -f 3 | sort | uniq | wc -l
3123249
zcat NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz | grep -v '#' | cut -f 3 | sort | uniq | wc -l
3123249
```

## Run the remapping pipeline 

See documentation in the [README](README.md) for how to perform the remapping.
* Complete Prerequisites and Testing the installation steps
* Activate variant-remapping Python environment
* Set `$GENOME_DIR` to the directory containing the genomes 
* Set `$REMAPPING_DIR` to the output directory

At the moment the pipeline is run this way: 

GRCh37 variants to GRCh38 genome:

```bash
mkdir ${REMAPPING_DIR}/GRCh37_to_GRCh38
nextflow run main.nf \
    --oldgenome ${GENOME_DIR}/hs37d5/hs37d5.fa \
    --newgenome ${GENOME_DIR}/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --vcffile ${GIAB_DIR}/NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz \
    --outfile ${REMAPPING_DIR}/GRCh37_to_GRCh38/NA12878.vcf
```

GRCh38 variants to GRCh37 genome:

```bash
mkdir ${REMAPPING_DIR}/GRCh38_to_GRCh37
nextflow run main.nf \
    --oldgenome ${GENOME_DIR}/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --newgenome ${GENOME_DIR}/hs37d5/hs37d5.fa \
    --vcffile ${GIAB_DIR}/NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz \
    --outfile ${REMAPPING_DIR}/GRCh38_to_GRCh37/NA12878.vcf
```

## Comparison between remapped and standard using Hap.py

### Run Hap.py with Singularity (needs access to EBI infrastructure)

Hap.py was installed using Singularity and can be accessed by setting `$SINGULARITY_DIR` to the location where singularity images are kept. See [this comment](https://www.ebi.ac.uk/panda/jira/browse/EVA-1835?focusedCommentId=312411&page=com.atlassian.jira.plugin.system.issuetabpanels%3Acomment-tabpanel#comment-312411) for details.

```bash
singularity exec \
    -B ${GIAB_DIR}:/giab \
    -B ${REMAPPING_DIR}:/remapping \
    -B ${GENOME_DIR}:/genomes \
    ${SINGULARITY_DIR}/hap.py.simg \
    /opt/hap.py/bin/hap.py  \
        /giab/NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz \
        /remapping/GRCh37_to_GRCh38/NA12878_f50_s0.6_d0.04_pre_final_vcf_fixed.vcf.gz \
        -f /giab/NA12878/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
        -r /genomes/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        -o /remapping/GRCh37_to_GRCh38/hap.py_assessement \
        -V
```

```bash
singularity exec \
    -B ${GIAB_DIR}:/giab \
    -B ${REMAPPING_DIR}:/remapping \
    -B ${GENOME_DIR}:/genomes \
    ${SINGULARITY_DIR}/hap.py.simg \
    /opt/hap.py/bin/hap.py  \
        /giab/NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz \
        /remapping/GRCh38_to_GRCh37/NA12878_f50_s0.6_d0.04_pre_final_vcf_fixed.vcf.gz \
        -f /giab/NA12878/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
        -r /genomes/hs37d5/hs37d5.fa \
        -o /remapping/GRCh38_to_GRCh37/hap.py_assessement \
        -V
```

### Run Hap.py with docker

Hap.py could be installed with docker through this [instruction](https://github.com/Illumina/hap.py#docker).

```bash
docker run -it \
    -v ${GIAB_DIR}:/giab \
    -v ${REMAPPING_DIR}:/remapping \
    -v ${GENOME_DIR}:/genomes \
    pkrusche/hap.py /opt/hap.py/bin/hap.py \
        /giab/NA12878/GRCh38/HG001_GRCh38_Fake_Genotypes_shared.vcf.gz \
        /remapping/GRCh37_to_GRCh38/NA12878.vcf \
        -f /giab/NA12878/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
        -r /genomes/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        -o /remapping/GRCh37_to_GRCh38/hap.py_assessement \
        -V
```

```bash
docker run -it \
    -v ${GIAB_DIR}:/giab \
    -v ${REMAPPING_DIR}:/remapping \
    -v ${GENOME_DIR}:/genomes \
    pkrusche/hap.py /opt/hap.py/bin/hap.py \
        /giab/NA12878/GRCh37/HG001_GRCh37_Fake_Genotypes_shared.vcf.gz \
        /remapping/GRCh38_to_GRCh37/NA12878.vcf \
        -f /giab/NA12878/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed \
        -r /genomes/hs37d5/hs37d5.fa \
        -o /remapping/GRCh38_to_GRCh37/hap.py_assessement \
        -V
```
