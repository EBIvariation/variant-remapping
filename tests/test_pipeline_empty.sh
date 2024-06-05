#!/bin/bash

set -Eeuo pipefail

function asserteq() {
  if [[ ! "$1" -eq "$2" ]]
  then
    echo "Assertion Error: $1 not equal to $2"
    exit 1
  fi

}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR=$(dirname $SCRIPT_DIR)

# Build the Source VCF
cat << EOT > "${SCRIPT_DIR}/resources/source_empty.vcf"
##fileformat=VCFv4.3
##INFO=<ID=COMMENT,Number=1,Type=String,Description="Comment">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype across all datasets with called genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG001
EOT

nextflow run ${SOURCE_DIR}/main.nf \
-config ${SCRIPT_DIR}/resources/config.yml \
--oldgenome ${SCRIPT_DIR}/resources/genome.fa \
--newgenome ${SCRIPT_DIR}/resources/new_genome.fa \
--vcffile ${SCRIPT_DIR}/resources/source_empty.vcf \
--outfile ${SCRIPT_DIR}/resources/remap_empty.vcf

# Check the presence of the output file
ls ${SCRIPT_DIR}/resources/remap_empty.vcf \
   ${SCRIPT_DIR}/resources/remap_empty_unmapped.vcf \
   ${SCRIPT_DIR}/resources/remap_empty_counts.yml

# Build the expected VCF
cat << EOT > "${SCRIPT_DIR}/resources/expected_remap.vcf"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG001
EOT

# Compare vs the expected VCF
diff "${SCRIPT_DIR}/resources/expected_remap.vcf" <(grep -v '^##' "${SCRIPT_DIR}/resources/remap_empty.vcf")

asserteq `cat ${SCRIPT_DIR}/resources/remap_empty_counts.yml | grep 'all:' | cut -d ' ' -f 2`  0
asserteq `cat ${SCRIPT_DIR}/resources/remap_empty_counts.yml | grep 'filtered:' | cut -d ' ' -f 2`  0


# Clean up after the test
rm -rf work .nextflow* \
       ${SCRIPT_DIR}/resources/source_empty.vcf \
       ${SCRIPT_DIR}/resources/expected_remap.vcf \
       ${SCRIPT_DIR}/resources/remap_empty.vcf \
       ${SCRIPT_DIR}/resources/remap_empty_counts.yml \
       ${SCRIPT_DIR}/resources/remap_empty_unmapped.vcf \
       ${SCRIPT_DIR}/resources/new_genome.fa.* \
       ${SCRIPT_DIR}/resources/genome.fa.fai
