#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam


def calculate_variant_position(read, flank_length):
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
    read_var_pos = start + flank_length + 1 + cig_list[0:flank_length + 1].count(pysam.CDEL)
    for operator in range(0, flank_length + 2 + local_region_size + cig_list[0:flank_length+1].count(pysam.CDEL)):
        # Stop incrementing varpos once we've reached read_var_pos:
        if start + operator < read_var_pos:
            # Mismatch/matches and deletions increase the position counter:
            if cig_list[operator] == pysam.CMATCH or cig_list[operator] == pysam.CDEL:
                varpos += 1
        # If we are in the local region around the variant, but not at the position of the variant:
        if ((start + operator >= (read_var_pos - local_region_size))
                and (start + operator <= (read_var_pos + local_region_size))
                and (start + operator != read_var_pos)):
            if cig_list[operator] == 0:  # Match or mismatch
                perf_counter += 1  # Increase the counter for perfect local region
    return varpos, perf_counter, local_region_size


def is_read_valid(read, counter, flank_length, score_cutoff, diff_cutoff):
    """
    Filter out reads that are:
        - Unmapped
        - Non primary (secondary): AS's that are too low and XS's that are too close to AS
    """
    if read.is_secondary:  # We only want to deal with primary reads
        return False
    if read.is_unmapped:  # Can be decoded with bitwise flag with & 4 (4 means unmapped)
        counter['unmapped'] += 1
        return False
    # Getting AS and XS:
    AS = read.get_tag('AS')
    try:
        XS = read.get_tag('XS')  # Some variants don't have secondary alignments, which throws an error
    except KeyError:
        XS = -flank_length  # Set an arbitrary low value for the 'artificial' secondary alignment
    if AS <= int(score_cutoff):
        counter['primary_poor'] += 1
        return False
    if AS - XS < diff_cutoff:
        counter['gap_small'] += 1
        return False
    return True


def is_local_alignment_perfect(counter, perf_counter, local_region_size):
    """
    Filter out non-perfect local alignment around the variant
    """
    if perf_counter != 2 * local_region_size:
        counter['context_bad'] += 1
        return False
    # Filter out AS's that are too low and XS's that are too close to AS, and those with a non-perfect local alignment:
    counter['remapped'] += 1
    return True


def fetch_bases(fasta, contig, length, start):
    """
    Returns a subsection from a specified FASTA contig. The start coordinate is 1-based.
    """
    zero_base_start = start - 1
    end = zero_base_start + length
    new_ref = fasta.fetch(reference=contig, start=zero_base_start, end=end)
    return new_ref


def calculate_new_alleles(old_ref, new_ref, old_alt, is_reverse_strand):
    """Given allele and strand information, calculate the new (ref, alt) pair."""
    # If reverse strand -> calculate complement
    old_ref_processed = old_ref
    old_alt_processed = old_alt
    if is_reverse_strand:
        old_ref_processed = Seq(old_ref_processed, generic_dna).reverse_complement()
        old_alt_processed = Seq(old_alt_processed, generic_dna).reverse_complement()

    # No changes
    if old_ref_processed == new_ref:
        return new_ref, old_alt_processed
    # Ref and Alt are the same -> Alt changed to old Ref
    if new_ref == old_alt_processed:
        return new_ref, old_ref_processed
    else:
        return new_ref, old_alt_processed


def group_reads(bam_file_path):
    with pysam.AlignmentFile(bam_file_path, 'rb') as inbam:
        current_read_name = None
        current_reads = []
        for read in inbam:
            if read.query_name == current_read_name:
                current_reads.append(read)
            else:
                if current_read_name:
                    yield current_reads
                current_reads = [read]
            current_read_name = read.query_name
        yield current_reads


def process_bam_file(bam_file_path, output_file, flank_length, score_perc, diff_AS_XS_perc, max_nb_alignment,
                     new_genome):

    # Calculate the score cutoff based on flanking seq length
    score_cutoff = -(flank_length * score_perc)
    diff_cutoff = flank_length * diff_AS_XS_perc
    counter = Counter()
    fasta = pysam.FastaFile(new_genome)
    with open(output_file, 'w') as outfile:
        for read_group in group_reads(bam_file_path):
            if len(read_group) < max_nb_alignment:
                counter['total'] += 1
                for read in read_group:
                    if is_read_valid(read, counter, flank_length, score_cutoff, diff_cutoff):
                        varpos, perf_counter, local_region_size = calculate_variant_position(read, flank_length)
                        if is_local_alignment_perfect(counter, perf_counter, local_region_size):
                            name = read.query_name
                            info = name.split('|')
                            old_ref = info[2]
                            old_alt = info[3]

                            new_ref = fetch_bases(fasta, read.reference_name, len(old_ref), varpos)
                            new_ref, new_alt = calculate_new_alleles(old_ref, new_ref, old_alt, read.is_reverse)

                            outfile.write(
                                '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (read.reference_name, varpos, info[4], new_ref, new_alt,
                                                                  info[5], info[6], info[7]))
            else:
                counter['Too many alignments'] += 1

    print(counter['total'])
    print(counter['Too many alignments'])
    print(counter['unmapped'])
    print(counter['primary_poor'])
    print(counter['gap_small'])
    print(counter['context_bad'])
    print(counter['remapped'])
    print('% of variants rejected for:')
    print('Unmapped: {:.2%}'.format(counter['unmapped'] / counter['total']))
    print('Primary alignment too poor: {:.2%}'.format(counter['primary_poor'] / counter['total']))
    print('Primary and secondary alignments too close: {:.2%}'.format(counter['gap_small'] / counter['total']))
    print('Local region around variant too poorly aligned: {:.2%}'.format(counter['context_bad'] / counter['total']))


def main():
    description = (
        'Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints'
        'everything to a new file\n')

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--bam', type=str, required=True,
                        help='bam file containing remapped variants ')
    parser.add_argument('-o', '--outfile', type=str, required=True,
                        help='name of new file')
    parser.add_argument('-f', '--flankingseqlength', type=int, required=True,
                        help='length of each of the flanking sequences')
    parser.add_argument('-s', '--scoreperc', type=float, required=True,
                        help='the alignment score cut off percentage of flanking seq' 'length (keeps values strictly above)')
    parser.add_argument('-d', '--difference_AS_XS', type=float, required=True,
                        help='difference threshold % between AS and XS')
    parser.add_argument('-x', '--max_alignment', type=int, required=True,
                        help='the maximum number of alignment authorised for the reads to be kept')
    parser.add_argument('--newgenome', required=True,
                        help='FASTA file of the target genome')
    args = parser.parse_args()

    process_bam_file(
        bam_file_path=args.bam,
        output_file=args.outfile,
        flank_length=args.flankingseqlength,
        score_perc=args.scoreperc,
        diff_AS_XS_perc=args.difference_AS_XS,
        max_nb_alignment=args.max_alignment,
        new_genome=args.newgenome
    )


if __name__ == '__main__':
    main()
