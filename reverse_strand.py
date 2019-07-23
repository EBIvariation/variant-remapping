#! /usr/bin/python3
import argparse
from argparse import RawTextHelpFormatter
import bamnostic as bs

description="Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything to a new file\nAlso filters out the non-mapped"

parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--bam", help="bam file containing remapped variants ")
parser.add_argument("-o", "--outfile", help="name of new file")
args = parser.parse_args()

bamfile=bs.AlignmentFile(args.bam, 'rb')
outfile=open(args.outfile, 'w')

#reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything into correct columns

for read in bamfile:
	name=read.read_name
	info=name.split("|")
	#mapped onto reverse strand:
	if read.flag==16:
		if info[2]=='A':
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+50,info[3],"T",info[4],info[5],info[6]))
		elif info[2]=='T':
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+50,info[3],"A",info[4],info[5],info[6]))
		elif info[2]=='C':
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+50,info[3],"G",info[4],info[5],info[6]))
		elif info[2]=='G':
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+50,info[3],"C",info[4],info[5],info[6]))
	#only mapped onto direct strand:
	elif read.flag==0:
		outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+50,info[3],info[2],info[4],info[5],info[6]))

outfile.close()