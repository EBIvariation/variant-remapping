#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from collections import Counter, defaultdict

import yaml
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pysam


def reverse_complement(sequence):
    return str(Seq(sequence, generic_dna).reverse_complement())


def calculate_new_variant_definition(left_read, right_read, ref_fasta, original_vcf_rec):
    """
    Resolve the variant definition from the flanking region alignment and old variant definition
    TODO: Link to algorithm description once public
    """
    # Flag to highlight low confidence in an event detected
    failure_reason = None

    # Grab information from the read name
    name = left_read.query_name

    old_ref = original_vcf_rec[3]
    old_alts = original_vcf_rec[4].split(',')

    operations = []
    # Define new ref and new pos
    new_ref = fetch_bases(ref_fasta, left_read.reference_name, left_read.reference_end + 1,
                          right_read.reference_start - left_read.reference_end)
    new_pos = left_read.reference_end + 1

    # TODO: apply the changes to the genotypes as well if they are present.
    # 1. Handle reference strand change
    if not left_read.is_reverse and not right_read.is_reverse:
        # Forward strand alignment
        old_ref_conv = old_ref
        old_alt_conv = old_alts
        operations.append('st=+')
    elif left_read.is_reverse and right_read.is_reverse:
        # Reverse strand alignment
        old_ref_conv = reverse_complement(old_ref)
        old_alt_conv = [reverse_complement(alt) for alt in old_alts]
        operations.append('st=-')
    else:
        # This case should be handled by the filtering but raise just in case...
        error_msg = (f'Impossible read configuration: '
                     f'read1 is_reverse: {left_read.is_reverse}, '
                     f'read2 is_reverse: {right_read.is_reverse}, '
                     f'read1 position: {left_read.pos}, '
                     f'read2 position: {right_read.pos}')
        raise ValueError(error_msg)

    # 2. Assign new allele sequences
    if new_ref == old_ref_conv:
        new_alts = old_alt_conv
    elif new_ref in old_alt_conv:
        old_alt_conv.remove(new_ref)
        new_alts = old_alt_conv
        new_alts.append(old_ref_conv)
        operations.append('rac=' + old_ref_conv + '-' + new_ref)
        if len(old_ref_conv) != len(new_ref):
            failure_reason = 'Reference Allele length change'
    else:
        new_alts = old_alt_conv
        new_alts.append(old_ref_conv)
        operations.append('rac=' + old_ref_conv + '-' + new_ref)
        operations.append('nra')
        failure_reason = 'Novel Reference Allele'

    # 3. Correct zero-length reference sequence
    if len(new_ref) == 0:
        new_pos -= 1
        new_ref = fetch_bases(ref_fasta, left_read.reference_name, new_pos, 1)
        new_alts = [new_ref + alt for alt in new_alts]
        operations.append('zlr')

    return new_pos, new_ref, new_alts, operations, failure_reason


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
    It will group reads by query name and create three subgroups of primary, supplementary and secondary aligned reads.
    It returns an iterators where each element is a tuple of the three lists
    :param bam_file_path: the name sorted bam file
    :return: iterator of tuples containing three lists
    """
    with pysam.AlignmentFile(bam_file_path, 'rb') as inbam:
        current_read_name = None
        primary_group = None
        secondary_group = None
        supplementary_group = None
        for read in inbam:
            if read.query_name == current_read_name:
                pass
            else:
                if current_read_name:
                    yield primary_group, supplementary_group, secondary_group
                primary_group = []
                secondary_group = []
                supplementary_group = []
            if read.is_secondary:
                secondary_group.append(read)
            elif read.is_supplementary:
                supplementary_group.append(read)
            else:
                primary_group.append(read)
            current_read_name = read.query_name
        if primary_group:
            yield primary_group, supplementary_group, secondary_group


def order_reads(primary_group, primary_to_supplementary):
    """
    Order read and return the most 5' (smallest coordinates) first.
    if a supplementary read exists and is closer to the other read then it is used in place of the primary
    """
    read1, read2 = primary_group
    suppl_read1 = suppl_read2 = None
    if read1 in primary_to_supplementary:
        suppl_read1 = primary_to_supplementary.get(read1)[0]
    if read2 in primary_to_supplementary:
        suppl_read2 = primary_to_supplementary.get(read2)[0]
    if read1.reference_start <= read2.reference_start:
        if suppl_read1 and suppl_read1.reference_start > read1.reference_start:
            read1 = suppl_read1
        if suppl_read2 and suppl_read2.reference_start < read2.reference_start:
            read2 = suppl_read2
        return read1, read2
    else:
        if suppl_read1 and suppl_read1.reference_start < read1.reference_start:
            read1 = suppl_read1
        if suppl_read2 and suppl_read2.reference_start > read2.reference_start:
            read2 = suppl_read2
        return read2, read1


def pass_basic_filtering(primary_group, secondary_group, primary_to_supplementary, counter, filter_align_with_secondary):
    """Test if the alignment pass basic filtering such as presence of secondary alignments, any primary unmapped,
    primary mapped on different chromosome, or primary mapped poorly."""
    if filter_align_with_secondary and len(secondary_group):
        counter['Too many alignments'] += 1
    elif len(primary_group) < 2 or any(read.is_unmapped for read in primary_group):
        counter['Flank unmapped'] += 1
    elif len(set(read.reference_name for read in primary_group)) != 1:
        counter['Different chromosomes'] += 1
    elif any(len(suppl) > 1 for suppl in primary_to_supplementary.values()):
        counter['Too many supplementary'] += 1
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
    # in CIGAR tuples the operation is coded as an integer
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    if left_read.cigartuples[-1][0] == pysam.CSOFT_CLIP or right_read.cigartuples[0][0] == pysam.CSOFT_CLIP:
        counter['Soft-clipped alignments'] += 1
    elif left_read.reference_end > right_read.reference_start:
        counter['Overlapping alignment'] += 1
    elif left_read.is_reverse != right_read.is_reverse:
        counter['Unexpected orientation'] += 1
    else:
        return True
    return False


def output_failed_alignment(original_vcf_rec, outfile):
    """
    Output the original VCF entry when alignment have failed to pass all thresholds
    """
    print('\t'.join(original_vcf_rec), file=outfile)


def link_supplementary(primary_group, supplementary_group):
    """Link supplementary alignments to their primary."""
    if not supplementary_group:
        # No supplementary so no linking required
        return {}
    supplementary_dict = {}
    primary_to_supplementary = defaultdict(list)
    for supplementary_read in supplementary_group:
        supplementary_dict[supplementary_read.reference_name + str(supplementary_read.reference_start + 1)] = supplementary_read
    for primary in primary_group:
        # chr2,808117,+,1211M790S,60,1;
        if primary.has_tag('SA'):
            for other_alignment in primary.get_tag('SA').split(';'):
                if other_alignment:
                    rname, pos = other_alignment.split(',')[:2]
                    primary_to_supplementary[primary].append(
                        supplementary_dict[rname + pos]
                    )
    return dict(primary_to_supplementary)


def get_vcf_entry(name, invcf_file_path):
    """Retrieve the full VCF record either from the read name or by fetching it in the original file"""
    # Remove the trailing underscore that have been added in convertVCFToBed
    name_info = name.rstrip('_').split('|')
    if len(name_info) < 4:
        # Fetch the line from the file
        with open(invcf_file_path) as open_file:
            open_file.seek(int(name_info[2]))
            return open_file.readline().strip().split('\t')
    else:
        return name_info[0:2] + name_info[3:]


def process_bam_file(invcf_file_path, bam_file_path, output_file, out_failed_file, new_genome,
                     filter_align_with_secondary, flank_length, summary_file):
    counter = Counter()
    fasta = pysam.FastaFile(new_genome)

    with open(output_file, 'w') as outfile, open(out_failed_file, 'w') as out_failed:
        for primary_group, supplementary_group, secondary_group in group_reads(bam_file_path):
            counter['total'] += 1
            primary_to_supplementary = link_supplementary(primary_group, supplementary_group)
            original_vcf_rec = get_vcf_entry(primary_group[0].query_name, invcf_file_path)
            if pass_basic_filtering(primary_group, secondary_group, primary_to_supplementary, counter, filter_align_with_secondary):
                left_read, right_read = order_reads(primary_group, primary_to_supplementary)
                if pass_aligned_filtering(left_read, right_read, counter):
                    varpos, new_ref, new_alts, ops, failure_reason = \
                        calculate_new_variant_definition(left_read, right_read, fasta, original_vcf_rec)
                    if not failure_reason:
                        counter['Remapped'] += 1
                        if original_vcf_rec[7] != '.':
                            original_vcf_rec[7] = ';'.join(original_vcf_rec[7].strip(';').split(';') + ops)
                        else:
                            original_vcf_rec[7] = ';'.join(ops)
                        outfile.write(
                            '%s\t%s\t%s\t%s\t%s\t%s\n' % (left_read.reference_name, varpos, original_vcf_rec[2],
                                                          new_ref, ','.join(new_alts), '\t'.join(original_vcf_rec[5:]))
                        )
                    else:
                        # Currently the alignment is not precise enough to ensure that the allele change for INDEL and
                        # novel reference allele are correct. So we skip them.
                        # TODO: add realignment confirmation see #14 and EVA-2417
                        counter[failure_reason] += 1
                        output_failed_alignment(original_vcf_rec, out_failed)
                else:
                    output_failed_alignment(original_vcf_rec, out_failed)
            else:
                output_failed_alignment(original_vcf_rec, out_failed)
    with open(summary_file, 'w') as open_summary:
        yaml.safe_dump({f'Flank_{flank_length}': dict(counter)}, open_summary)


def main():
    description = ('Process alignment results in bam format to determine the location of the variant in the new genome.'
                   ' Each variant will be either output in the new genome VCF or the old VCF will be output in a '
                   'separate file.')

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-v', '--vcf', type=str, required=True,
                        help='original vcf file containing the source data')
    parser.add_argument('-i', '--bam', type=str, required=True,
                        help='Input BAM file with remapped flanking regions')
    parser.add_argument('-o', '--outfile', type=str, required=True,
                        help='Output VCF file with remapped variants')
    parser.add_argument('--out_failed_file', type=str, required=True,
                        help='Name of the file containing reads that did not align correctly')
    parser.add_argument('--flank_length', type=int, required=True,
                        help='Length of the flanking region used.')
    parser.add_argument('--summary', type=str, required=True,
                        help='YAML files containing the summary metrics')

    parser.add_argument('-f', '--filter_align_with_secondary', action='store_true', default=False,
                        help='Filter out alignments that have one or several secondary alignments.')
    parser.add_argument('-n', '--newgenome', required=True, help='FASTA file of the target genome')
    args = parser.parse_args()

    process_bam_file(
        invcf_file_path=args.vcf,
        bam_file_path=args.bam,
        output_file=args.outfile,
        out_failed_file=args.out_failed_file,
        new_genome=args.newgenome,
        filter_align_with_secondary=args.filter_align_with_secondary,
        flank_length=args.flank_length,
        summary_file=args.summary
    )


if __name__ == '__main__':
    main()
