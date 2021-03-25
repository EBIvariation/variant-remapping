#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam


def reverse_complement(sequence):
    return str(Seq(sequence, generic_dna).reverse_complement())


def calculate_new_variant_definition(left_read, right_read, ref_fasta):
    """
    Resolve the variant definition from the flanking region alignment and old variant definition
    TODO: Link to algorithm description once public
    """
    # Grab information from the read name
    name = left_read.query_name
    info = name.split('|')
    old_ref = info[2]
    old_alts = info[3].split(',')

    # Define new ref and new pos
    new_ref = fetch_bases(ref_fasta, left_read.reference_name, left_read.reference_end + 1,
                          right_read.pos - left_read.reference_end)
    new_pos = left_read.reference_end + 1

    # TODO: All th operation bellow should be recorded so they can be used for subsequent changes of the genotypes.
    # 1. Handle reference strand change
    if not left_read.is_reverse and not right_read.is_reverse:
        # Forward strand alignment
        old_ref_conv = old_ref
        old_alt_conv = old_alts
    elif left_read.is_reverse and right_read.is_reverse:
        # Reverse strand alignment
        old_ref_conv = reverse_complement(old_ref)
        old_alt_conv = [reverse_complement(alt) for alt in old_alts]
    else:
        error_msg = (f'Impossible read configuration: '
                     f'read1 is_reverse: {left_read.is_reverse}, '
                     f'read2 is_reverse: {right_read.is_reverse}, '
                     f'read1 position: {left_read.pos}, '
                     f'read2 position: {right_read.pos}')
        print(error_msg)
        raise ValueError(error_msg)

    # 2. Assign new allele sequences
    if new_ref == old_ref_conv:
        new_alts = old_alt_conv
    elif new_ref in old_alt_conv:
        old_alt_conv.remove(new_ref)
        new_alts = old_alt_conv
        new_alts.append(old_ref_conv)
    else:
        new_alts = old_alt_conv
        new_alts.append(old_ref_conv)

    # 3. Correct zero-length reference sequence
    if len(new_ref) == 0:
        new_pos -= 1
        new_ref = fetch_bases(ref_fasta, left_read.reference_name, new_pos, 1)
        new_alts = [new_ref + alt for alt in new_alts]

    return new_pos, new_ref, new_alts


def fetch_bases(fasta, contig, start, length):
    """
    Returns a subsection from a specified FASTA contig. The start coordinate is 1-based.
    """
    zero_base_start = start - 1
    end = zero_base_start + length
    new_ref = fasta.fetch(reference=contig, start=zero_base_start, end=end)
    return new_ref


def group_reads(bam_file_path):
    """
    This function assumes that the reads are sorted by query name.
    It will group reads by query name and create two subgroups of primary and secondary aligned reads.
    It returns an iterators where each element is a tuple of two lists
    :param bam_file_path: the name sorted bam file
    :return: iterator of tuples containing two lists
    """
    with pysam.AlignmentFile(bam_file_path, 'rb') as inbam:
        current_read_name = None
        primary_group = None
        secondary_group = None
        for read in inbam:
            if read.query_name == current_read_name:
                pass
            else:
                if current_read_name:
                    yield primary_group, secondary_group
                primary_group = []
                secondary_group = []
            if read.is_secondary:
                secondary_group.append(read)
            else:
                primary_group.append(read)
            current_read_name = read.query_name
        yield primary_group, secondary_group


def _order_reads(group):
    """Order read and return the most 5' (smallest coordinates) first."""
    read1, read2 = group
    if read1.pos <= read2.pos:
        return read1, read2
    else:
        return read2, read1


def pass_basic_filtering(primary_group, secondary_group, counter, filter_align_with_secondary, alignment_score_threshold):
    """Test if the alignment pass basic filtering such as presence of secondary alignments, any primary unmapped,
    primary mapped on different chromosome, or primary mapped poorly."""
    if filter_align_with_secondary and len(secondary_group):
        counter['Too many alignments'] += 1
    elif len(primary_group) < 2 or any(read.is_unmapped for read in primary_group):
        counter['Flank unmapped'] += 1
    elif len(set(read.reference_name for read in primary_group)) != 1:
        counter['Different chromosomes'] += 1
    elif any(read.get_tag('AS') <= int(alignment_score_threshold) for read in primary_group):
        counter['Poor alignment'] += 1
    else:
        return True
    return False


def pass_aligned_filtering(left_read, right_read, counter):
    """
    Test if the two reads pass the additional filters such as check for soft-clipped end next to the variant region,
    or overlapping region between the two reads.
    :param left_read: the left (or 5') most read
    :param right_read: the right (or 3') most read
    :param counter: Counter to report the number of reads filtered.
    :return: True or False
    """
    if left_read.cigartuples[-1][1] == 'S' or right_read.cigartuples[0][1] == 'S':
        counter['Soft-clipped alignments'] += 1
    elif left_read.reference_end > right_read.pos:
        counter['Overlapping alignment'] += 1
    else:
        return True
    return False


def output_failed_alignment(primary_group, outfile):
    """
    Output the original VCF entry when alignment have failed to pass all thresholds
    """
    info = primary_group[0].query_name.split('|')
    print('\t'.join(info), file=outfile)


def process_bam_file(bam_file_path, output_file, out_failed_file, new_genome, filter_align_with_secondary, alignment_score_threshold):

    # Calculate the score cutoff based on flanking seq length
    counter = Counter()
    fasta = pysam.FastaFile(new_genome)
    with open(output_file, 'w') as outfile, open(out_failed_file, 'w') as out_failed:
        for primary_group, secondary_group in group_reads(bam_file_path):
            counter['total'] += 1
            if pass_basic_filtering(primary_group, secondary_group, counter, filter_align_with_secondary,
                                    alignment_score_threshold):
                left_read, right_read = _order_reads(primary_group)
                if pass_aligned_filtering(left_read, right_read, counter):
                    counter['Remapped'] += 1
                    varpos, new_ref, new_alts = calculate_new_variant_definition(left_read, right_read, fasta)
                    info = left_read.query_name.split('|')
                    outfile.write(
                        '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (left_read.reference_name, varpos, info[4], new_ref,
                                                              ','.join(new_alts), info[5], info[6], info[7]))
                else:
                    output_failed_alignment(primary_group, out_failed)
            else:
                output_failed_alignment(primary_group, out_failed)
    for metric in ['Too many alignments', 'Flank unmapped', 'Different chromosomes',
                   'Poor alignment', 'Overlapping alignment', 'Soft-clipped alignments', 'Remapped']:
        print(f'{counter.get(metric, 0)} ({counter.get(metric, 0)/counter.get("total"):.2%}) variants rejected for {metric}')


def main():
    description = ('Process alignment results in bam format to determine the location of the variant in the new genome.'
                   ' Each variant will be either output in the new genome VCF or the old VCF will be output in a '
                   'separate file.')

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--bam', type=str, required=True,
                        help='bam file containing remapped variants ')
    parser.add_argument('-o', '--outfile', type=str, required=True,
                        help='name of new  VCF file')
    parser.add_argument('--out_failed_file', type=str, required=True,
                        help='name of the file containing reads that did not align correctly')
    parser.add_argument('-a', '--alignment_score_threshold', type=int, default=-1,
                        help='')
    parser.add_argument('-f', '--filter_align_with_secondary', action='store_true', default=False,
                        help='Filter out alignment that have one or several secondary alignments.')
    parser.add_argument('-n', '--newgenome', required=True, help='Fasta file of the target genome')
    args = parser.parse_args()

    process_bam_file(
        bam_file_path=args.bam,
        output_file=args.outfile,
        out_failed_file=args.out_failed_file,
        new_genome=args.newgenome,
        filter_align_with_secondary=args.filter_align_with_secondary,
        alignment_score_threshold=args.alignment_score_threshold
    )


if __name__ == '__main__':
    main()
