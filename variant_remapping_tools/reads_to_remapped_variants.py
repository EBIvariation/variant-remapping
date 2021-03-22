#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam


def reverse_complement(sequence):
    return Seq(sequence, generic_dna).reverse_complement()


def calculate_new_variant_definition(read1, read2, ref_fasta):
    """
    Resolve the variant definition from the flanking region alignment and old variant definition
    """
    # grab information from the read name
    name = read1.query_name
    info = name.split('|')
    old_ref = info[2]
    old_alts = info[3].split(',')

    # 1. Handle reference strand change
    if not read1.is_reverse and not read2.is_reverse and read1.pos < read2.pos:
        # Forward strand alignment
        new_ref = fetch_bases(ref_fasta, read1.reference_name, read1.reference_end + 1, read2.pos - read1.reference_end)
        new_pos = read1.reference_end + 1
        old_ref_conv = old_ref
        old_alt_conv = old_alts

    elif read1.is_reverse and read2.is_reverse and read1.pos > read2.pos:
        # Reverse strand alignment
        new_ref = fetch_bases(ref_fasta, read1.reference_name, read2.reference_end + 1, read1.pos - read2.reference_end)
        new_pos = read2.reference_end + 1
        old_ref_conv = reverse_complement(old_ref)
        old_alt_conv = [reverse_complement(alt) for alt in old_alts]
    else:
        error_msg = (f'Impossible read configuration: '
                     f'read1 is_reverse: {read1.is_reverse}, '
                     f'read2 is_reverse: {read2.is_reverse}, '
                     f'read1 position: {read1.pos}, '
                     f'read2 position: {read2.pos}')
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
        new_ref = fetch_bases(ref_fasta, new_pos, 1)
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
    """Order read in a group of two"""
    read1, read2 = group
    if read1.is_read1 and read2.is_read2:
        return read1, read2
    else:
        return read2, read1


def pass_filtering(primary_group, secondary_group, counter, filter_align_with_secondary, alignment_score_threshold):
    if filter_align_with_secondary and len(secondary_group):
        counter['Too many alignments'] += 1
    elif len(primary_group) < 2 or all(read.is_unmapped for read in primary_group):
        counter['Flank unmapped'] += 1
    elif len(set(read.reference_name for read in primary_group)) != 1:
        counter['Different chromosomes'] += 1
    elif any(read.get_tag('AS') <= int(alignment_score_threshold) for read in primary_group):
        counter['Poor alignment'] += 1
    else:
        return True
    return False


def process_bam_file(bam_file_path, output_file, new_genome, filter_align_with_secondary, alignment_score_threshold):

    # Calculate the score cutoff based on flanking seq length
    counter = Counter()
    fasta = pysam.FastaFile(new_genome)
    with open(output_file, 'w') as outfile:
        for primary_group, secondary_group in group_reads(bam_file_path):
            counter['total'] += 1
            if pass_filtering(primary_group, secondary_group, counter, filter_align_with_secondary, alignment_score_threshold):
                read1, read2 = _order_reads(primary_group)
                varpos, new_ref, new_alts = calculate_new_variant_definition(read1, read2, fasta)
                info = read1.query_name.split('|')
                outfile.write(
                    '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (read1.reference_name, varpos, info[4], new_ref,
                                                          ','.join(new_alts), info[5], info[6], info[7]))
    for metric in ['Too many alignments', 'Flank unmapped', 'Different chromosomes', 'Poor alignment']:
        print(f'{counter.get(metric, 0)} ({counter.get(metric, 0)/counter.get("total"):.2%}) variants rejected for {metric}')


def main():
    description = (
        'Reverses the allele for variants that got mapped onto the reverse strand of the new genome, and prints'
        'everything to a new file\n')

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--bam', type=str, required=True,
                        help='bam file containing remapped variants ')
    parser.add_argument('-o', '--outfile', type=str, required=True,
                        help='name of new file')
    parser.add_argument('-a', '--alignment_score_threshold', type=int, default=0,
                        help='')
    parser.add_argument('--filter_align_with_secondary', action='store_true', default=False,
                        help='Filter out alignment that have one or several secondary alignments.')
    parser.add_argument('--newgenome', required=True,
                        help='FASTA file of the target genome')
    args = parser.parse_args()

    process_bam_file(
        bam_file_path=args.bam,
        output_file=args.outfile,
        new_genome=args.newgenome,
        filter_align_with_secondary=args.filter_align_with_secondary,
        alignment_score_threshold=args.alignment_score_threshold
    )


if __name__ == '__main__':
    main()
