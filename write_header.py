#! /usr/bin/python3
import argparse
from argparse import RawTextHelpFormatter

description="Adds the first part of a VCF header up until the end of the ##INFOs, \nthen adds the newly generated ##configs and ##reference, and finally add the rest of the header."

parser = argparse.ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--old_header", help="old header containing ##fileformat, ##INFO etc")
parser.add_argument("-c", "--contigs", help="file containing ##contigs and ##reference")
parser.add_argument("-o", "--outfile", help="name of new renamed header file")
args = parser.parse_args()

newheader=open(args.outfile, 'w')

#find the index of the last line starting with ##INFO:
with open(args.old_header, 'r') as oldheader:
	for (i,line) in enumerate(oldheader):
		if line.startswith('##INFO'):
			p=i
			
#write all the lines before line p
with open(args.old_header, 'r') as oldheader:
	for (i,line) in enumerate(oldheader):
		if i <= p:
			newheader.write(line)
		else:
			break

#write the contig file
with open(args.contigs, 'r') as contigs:
	for line in contigs:
		newheader.write(line)

#write the rest of the first file
with open(args.old_header, 'r') as oldheader:
	for (i,line) in enumerate(oldheader):
		if i > p:
			newheader.write(line)

newheader.close()
