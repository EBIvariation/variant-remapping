#!/bin/bash

set -Eeuxo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR=$(dirname $SCRIPT_DIR)

nextflow run ${SOURCE_DIR}/remapping_pipeline.nf \
--oldgenome ${SCRIPT_DIR}/resources/genome.fa \
--newgenome ${SCRIPT_DIR}/resources/new_genome.fa \
--vcffile ${SCRIPT_DIR}/resources/source.vcf \
--outfile ${SCRIPT_DIR}/resources/remap.vcf

# Check the presence of the output file
ls ${SCRIPT_DIR}/resources/remap.vcf \
   ${SCRIPT_DIR}/resources/remap.vcf.stats \
   ${SCRIPT_DIR}/resources/new_genome.fa.1.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.2.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.3.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.4.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.rev.1.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.rev.2.bt2 \
   ${SCRIPT_DIR}/resources/new_genome.fa.fai \
   ${SCRIPT_DIR}/resources/genome.fa.fai

# Clean up after the test
rm -rf work .nextflow* ${SCRIPT_DIR}/resources/remap.vcf* ${SCRIPT_DIR}/resources/new_genome.fa.* ${SCRIPT_DIR}/resources/genome.fa.fai
