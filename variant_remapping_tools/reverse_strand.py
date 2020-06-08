#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam

description = ("Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints"
                "everything to a new file\n")

parser = argparse.ArgumentParser(description = description, formatter_class = RawTextHelpFormatter)
parser.add_argument("-i", "--bam", help = "bam file containing remapped variants ")
parser.add_argument("-o", "--outfile", help = "name of new file")
parser.add_argument("-p", "--old_ref_alleles", help = "name of output old ref alleles")
parser.add_argument("-f", "--flankingseqlength", type = int, help = "length of each of the flanking sequences")
parser.add_argument("-s", "--scoreperc", type = float, help = "the alignment score cut off percentage of flanking seq" "length (keeps values strictly above)")
parser.add_argument("-d", "--difference_AS_XS", type = float, help = "difference threshold % between AS and XS")

args = parser.parse_args()


bamfile = pysam.AlignmentFile(args.bam, 'rb')
outfile = open(args.outfile, 'w')
old_ref_alleles = open(args.old_ref_alleles, 'w')
flank_length = args.flankingseqlength
score_perc = args.scoreperc
diff_AS_XS_perc = args.difference_AS_XS
# Calculate the score cutoff based on flanking seq length: 
score_cutoff = -(flank_length*score_perc)
diff_cutoff = flank_length*diff_AS_XS_perc
unmapped_count = 0
primary_poor_count = 0
gap_small_count = 0
context_bad_count = 0
total_count = 0
remapped_count = 0

# Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints everything 
# into correct columns

for read in bamfile:
    name = read.query_name
    info = name.split("|")
    nucl = info[3]
    if read.is_unmapped:  # Can be decoded with bitwise flag with & 4 (4 means unmapped)
        unmapped_count += 1
        total_count += 1
        continue  # `not(read.is_unmapped)` doesn't work, the unmapped reads still get through, so a continue is needed
    if read.is_secondary:  # We only want to deal with primary reads
        continue
    total_count += 1
    # Mapped onto reverse strand:
    if read.is_reverse:  # Can be decoded with bitwise flag with & 16
        nucl = Seq(nucl, generic_dna).complement()
    # Getting AS and XS:
    AS = read.get_tag("AS")
    try:
        XS = read.get_tag("XS")  # Some variants don't have secondary alignments, which throws an error
    except KeyError:
        XS = -flank_length  # Set an arbitrary low value for the "artificial" secondary alignment 

    # Calculating correct position and filtering based on alignment around the variant position
    start = read.pos
    perf_counter = 0  # This counter will count the number of matches + mismatches around the variant position
    varpos = start  # Store the alignment start position
    cig_list = []  # This list will contain the CIGAR values: 10M will become [0,0,0,0,0,0,0,0,0,0]
    begin = 0
    local_region_size = 5
    # Create the list:
    for (cig_type, cig_length) in read.cigartuples:
        for i in range(begin, begin + cig_length):
            cig_list.append(cig_type)
        begin += cig_length
    # Deletions need to be counted in the new variant position:
    read_var_pos = start + flank_length + 1 + cig_list[0:flank_length+1].count(pysam.CDEL)
    for operator in range(0, flank_length + 2 + local_region_size + cig_list[0:flank_length+1].count(pysam.CDEL)):
        # Stop incrementing varpos once we've reached read_var_pos:
        if start + operator < read_var_pos:
            # Mismatch/matches and deletions increase the position counter:
            if cig_list[operator] == pysam.CMATCH or cig_list[operator] == pysam.CDEL: 
                varpos += 1
        # If we are in the local region around the variant, but not at the position of the variant:
        if ((start + operator >= (read_var_pos - local_region_size)) 
            and (start + operator <= (read_var_pos + local_region_size)) 
            and (start + operator != (read_var_pos))):
            if cig_list[operator] == 0:  # Match or mismatch
                perf_counter += 1  # Increase the counter for perfect local region
    if (AS <= int(score_cutoff)):
        primary_poor_count += 1
    elif (AS - XS < diff_cutoff):
        gap_small_count += 1
    elif (perf_counter != 2*local_region_size):
        context_bad_count += 1
    # Filter out AS's that are too low and XS's that are too close to AS, and those with a non-perfect local alignment:
    if (AS > int(score_cutoff)) and (AS - XS >= diff_cutoff) and (perf_counter == 2*local_region_size):
        remapped_count += 1
        # Write it all to the file:
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (read.reference_name,  varpos, info[4], nucl, info[5], info[6],\
                                                        info[7]))
        # Store old reference allele:
        old_ref_alleles.write("%s\n" % info[2])
outfile.close()
old_ref_alleles.close()
print(total_count)
print(unmapped_count)
print(primary_poor_count)
print(gap_small_count)
print(context_bad_count)
print(remapped_count)
print("% of variants rejected for:")
print("Unmapped: {:.2%}".format(unmapped_count/total_count))
print("Primary alignment too poor: {:.2%}".format(primary_poor_count/total_count))
print("Primary and secondary alignments too close: {:.2%}".format(gap_small_count/total_count))
print("Local region around variant too poorly aligned: {:.2%}".format(context_bad_count/total_count))