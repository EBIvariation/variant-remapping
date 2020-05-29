#!/bin/bash

set -e 

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR=$(dirname $SCRIPT_DIR)

nextflow run ${SOURCE_DIR}/genome_preparation.nf \
--fasta ${SOURCE_DIR}/resources/genome.fa \
--outdir ${SOURCE_DIR}/resources/


nextflow run ${SOURCE_DIR}/remapping_pipeline.nf \
--oldgenome ${SOURCE_DIR}/resources/genome.fa \
--newgenome ${SOURCE_DIR}/resources/genome.fa \
--vcffile ${SOURCE_DIR}/resources/source_vcf.vcf \
--outfile ${SOURCE_DIR}/resources/remap_vcf.vcf

ls ${SOURCE_DIR}/resources/remap_vcf.vcf \
   ${SOURCE_DIR}/resources/remap_vcf.vcf.stats \
   ${SOURCE_DIR}/resources/genome.fa.1.bt2 \
   ${SOURCE_DIR}/resources/genome.fa.2.bt2 \
   ${SOURCE_DIR}/resources/genome.fa.3.bt2 \
   ${SOURCE_DIR}/resources/genome.fa.4.bt2 \
   ${SOURCE_DIR}/resources/genome.fa.chrom.sizes \
   ${SOURCE_DIR}/resources/genome.fa.fai \
   ${SOURCE_DIR}/resources/genome.fa.rev.1.bt2 \
   ${SOURCE_DIR}/resources/genome.fa.rev.2.bt2

rm -rf work .nextflow* ${SOURCE_DIR}/resources/remap_vcf.vcf* ${SOURCE_DIR}/resources/genome.fa.*