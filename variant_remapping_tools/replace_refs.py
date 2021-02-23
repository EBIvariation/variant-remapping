#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import linecache


def process_variant(line, line_num, oldrefalleles_file):
    variant = line.split('\t')
    # Gets corresponding REF allele from the oldrefalleles file.
    old_ref = linecache.getline(oldrefalleles_file, line_num + 1).rstrip()
    # The reference allele extracted from the new genome in insertReferenceAllele and stored in variant[3]
    # is only one base long. There test "variant[3] == variant[4]" can only work if the ref and alt have
    # the same length which won't be the case for indels.
    # For insertion, we're leaving the single base extracted from the new genome (regardless of the alt)
    # for deletion, we're using the old reference sequence (regardless of the alt)
    # TODO: Correct this behaviour to have a consistent treatment of all variants types
    if variant[3] == variant[4] or len(old_ref) != 1:
        return '\t'.join((str(e) for e in [variant[0], variant[1], variant[2], old_ref, variant[4], variant[5], variant[6], variant[7]]))
    else:
        return line

if __name__ == "__main__":
    description = "replaces the REF allele with the old REF allele if REF = ALT in the remapped VCF"

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--vcf", help="VCF file containing remapped variants")
    parser.add_argument("-r", "--oldrefalleles", help="file containing the same header as the VCF file + the old "
                                                      "REF alleles in the same order as the VCF")
    parser.add_argument("-o", "--outfile", help="name of new file")
    args = parser.parse_args()

    with open(args.vcf, 'r') as vcf, open(args.outfile, 'w') as out_vcf:
        for(i, line) in enumerate(vcf):
            if line.startswith('#'):
                out_vcf.write(line)
            else:
                out_vcf.write(process_variant(line, i, args.oldrefalleles))
