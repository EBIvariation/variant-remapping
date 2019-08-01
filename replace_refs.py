#! /usr/bin/python3
import argparse
from argparse import RawTextHelpFormatter
import linecache

description="replaces the REF allele with the old REF allele if REF = ALT in the remapped VCF"

parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--vcf", help="VCF file containing remapped variants")
parser.add_argument("-r", "--oldrefalleles", help="file containing the same header as the VCF file + the old REF alleles in the same order as the VCF")
parser.add_argument("-o", "--outfile", help="name of new file")
args = parser.parse_args()

out_vcf=open(args.outfile, 'w')

with open(args.vcf, 'r') as vcf:
        for(i, line) in enumerate(vcf):
                if line.startswith('#'):
                        out_vcf.write(line)
                else:
                        variant=line.split("\t")
                        if variant[3] == variant[4]:
                        	#gets corresponding REF allele from the oldrefalleles file
                                out_vcf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (variant[0],variant[1],variant[2],linecache.getline(args.oldrefalleles, i+1).rstrip(),variant[4],variant[5],variant[6],variant[7]))
                        else:
                                out_vcf.write(line)
out_vcf.close()
