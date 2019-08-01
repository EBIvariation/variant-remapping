#! /home/juliette/miniconda3/bin/python3
import argparse
from argparse import RawTextHelpFormatter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam

description="Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything to a new file\n"

parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--bam", help="bam file containing remapped variants ")
parser.add_argument("-o", "--outfile", help="name of new file")
parser.add_argument("-p", "--old_ref_alleles", help="name of output old ref alleles")
args = parser.parse_args()


bamfile=pysam.AlignmentFile(args.bam, 'rb')
outfile=open(args.outfile, 'w')
old_ref_alleles=open(args.old_ref_alleles, 'w')

#reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything into correct columns

for read in bamfile:
		name=read.query_name
		info=name.split("|")
		nucl=info[3]
		#mapped onto reverse strand:
		if read.flag & 16: #decoding bitwise flag with &; may need to be changed back to ==
			nucl = Seq(nucl, generic_dna).complement()
		#write it all to the file:
		if read.flag != 4: #4 means unmapped
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,read.pos+51,info[4],nucl,info[5],info[6],info[7]))
			#store old reference allele:
			old_ref_alleles.write("%s\n" % info[2])
outfile.close()
old_ref_alleles.close()