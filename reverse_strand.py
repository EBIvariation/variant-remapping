#! /usr/bin/python3
import argparse
from argparse import RawTextHelpFormatter
import bamnostic as bs

description="Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything to a new file\nAlso filters out the non-mapped"

parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--bam", help="bam file containing remapped variants ")
parser.add_argument("-o", "--outfile", help="name of new file")
parser.add_argument("-p", "--old_ref_alleles", help="name of output old ref alleles")
args = parser.parse_args()


bamfile=bs.AlignmentFile(args.bam, 'rb')
outfile=open(args.outfile, 'w')
old_ref_alleles=open(args.old_ref_alleles, 'w')

#reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything into correct columns

for read in bamfile:
	name=read.read_name
	info=name.split("|")
	#mapped onto reverse strand:
	if read.flag & 16: #decoding bitwise flag with &; may need to be changed back to ==
		if info[2]=='A':
			rev_nucl="T"
		elif info[2]=='T':
			rev_nucl="A"
		elif info[2]=='C':
			rev_nucl="G"
		elif info[2]=='G':
			rev_nucl="C"
		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+51,info[4],rev_nucl,info[5],info[6],info[7]))
		#store old reference allele:
		old_ref_alleles.write(info[2])
	#only mapped onto direct strand:
	elif read.flag==0:
		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+51,info[4],info[3],info[5],info[6],info[7]))
		#store old reference allele:
		old_ref_alleles.write(info[2])

outfile.close()
old_ref_alleles.close()