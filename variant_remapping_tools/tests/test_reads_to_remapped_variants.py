import os
from collections import Counter
from unittest.mock import Mock, patch

import pysam
from unittest import TestCase
from variant_remapping_tools.reads_to_remapped_variants import fetch_bases, process_bam_file, \
    calculate_new_variant_definition, _order_reads, link_supplementary, pass_aligned_filtering, group_reads


class TestProcess(TestCase):

    def mk_read(self, tags=None, **kwargs):
        if not tags:
            tags = {}
        read = Mock(**kwargs)
        if tags:
            read.get_tag = lambda x: tags[x]
        read.has_tag = lambda x: x in tags
        return read

    def test_fetch_old_bases(self):
        fasta_path = self.get_test_resource('genome.fa')
        fasta = pysam.FastaFile(fasta_path)
        assert fetch_bases(fasta=fasta, contig='chr1', start=1, length=2) == 'CA'
        assert fetch_bases(fasta=fasta, contig='chr1', start=48, length=1) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr1', start=98, length=1) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr1', start=353, length=3) == 'GGA'
        assert fetch_bases(fasta=fasta, contig='chr1', start=1078, length=1) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr1', start=1404, length=1) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr1', start=1818, length=3) == 'AAC'
        assert fetch_bases(fasta=fasta, contig='chr1', start=2030, length=1) == 'A'

    def test_fetch_bases(self):
        fasta_path = self.get_test_resource('new_genome.fa')
        fasta = pysam.FastaFile(fasta_path)
        assert fetch_bases(fasta=fasta, contig='chr2', start=1, length=2) == 'CA'
        assert fetch_bases(fasta=fasta, contig='chr2', start=98, length=1) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr2', start=353, length=1) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr2', start=1078, length=1) == 'A'
        assert fetch_bases(fasta=fasta, contig='chr2', start=1404, length=1) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr2', start=1818, length=3) == 'AAC'
        assert fetch_bases(fasta=fasta, contig='chr2', start=2030, length=1) == 'A'

    def test_process_bam_file(self):
        """
        Given:
        chr1	48	.	C	A	50	PASS	.	GT	1/1
        chr1	98	.	C	CG	50	PASS	.	GT	1/1
        chr1	1078	.	G	A	50	PASS	.	GT	1/1
        chr1	1818	.	AAC	T	50	PASS	.	GT	1/1
        chr1	2030	.	A	TCC	50	PASS	.	GT	1/1
        """

        bamfile = self.get_test_resource("reads_aligned_test.bam")
        fasta_path = self.get_test_resource('new_genome.fa')
        output_file = '/tmp/remapped.vcf'
        summary_file = '/tmp/summary.yml'
        out_failed_file = '/tmp/unmapped.vcf'
        process_bam_file(bamfile, output_file, out_failed_file, fasta_path, True, 50, summary_file)

        expected = [
            'chr2	48	.	C	A	50	PASS	st=+\n',
            'chr2	48	.	C	T	50	PASS	st=+\n',
            'chr2	98	.	C	CG	50	PASS	st=+\n',
            'chr2	1078	.	A	G	50	PASS	st=+;rac=G-A\n',
            'chr2	1818	.	AAC	A	50	PASS	st=+\n',
            'chr2	2030	.	A	TCC	50	PASS	st=+\n'
        ]
        with open(output_file, 'r') as vcf:
            for(i, line) in enumerate(vcf):
                assert line == expected[i]
            # Expected and Generated VCF have the same number of lines
            assert i+1 == len(expected)

    def test_pass_aligned_filtering(self):
        left_read = self.mk_read(reference_start=0, reference_end=47, cigartuples=[(0, 47)])
        right_read = self.mk_read(reference_start=51, reference_end=95,  cigartuples=[(4, 3), (0, 44)])
        counter = Counter()
        assert not pass_aligned_filtering(left_read, right_read, counter)
        assert counter['Soft-clipped alignments'] == 1

    def test_link_supplementary(self):
        primary_group = [
            self.mk_read(tags={'SA': 'chr2,100,+;'}),
            self.mk_read()
        ]
        supplementary_group = [
            self.mk_read(reference_name='chr2', reference_start=99)
        ]
        assert link_supplementary(primary_group, supplementary_group) == {primary_group[0]: [supplementary_group[0]]}

    def test_order_reads(self):
        primary_group = [
            self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False),
            self.mk_read(reference_name='chr2', reference_start=48, reference_end=58, is_reverse=False)
        ]
        primary_to_supplementary = {}
        left_read, right_read = _order_reads(primary_group, primary_to_supplementary)
        assert left_read == primary_group[0]
        assert right_read == primary_group[1]

        primary_group = [
            self.mk_read(reference_name='chr2', reference_start=48, reference_end=58, is_reverse=True),
            self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        ]
        left_read, right_read = _order_reads(primary_group, primary_to_supplementary)
        assert left_read == primary_group[1]
        assert right_read == primary_group[0]

        primary_group = [
            self.mk_read(reference_name='chr2', reference_start=1, reference_end=10, is_reverse=False),
            self.mk_read(reference_name='chr2', reference_start=48, reference_end=58, is_reverse=False)
        ]
        supplementary_read = self.mk_read(reference_name='chr2', reference_start=20, reference_end=47, is_reverse=False)
        primary_to_supplementary = {primary_group[0]: [supplementary_read]}
        left_read, right_read = _order_reads(primary_group, primary_to_supplementary)
        assert left_read == supplementary_read
        assert right_read == primary_group[1]

    def test_calculate_new_variant_definition(self):
        fasta = 'fasta_path'

        # Forward strand alignment for SNP
        left_read = self.mk_read(query_name='chr1|48|C|A', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(query_name='chr1|48|C|A', reference_name='chr2', reference_start=48, reference_end=108, is_reverse=False)
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='C'):
            assert calculate_new_variant_definition(left_read, right_read, fasta) == (48, 'C', ['A'], ['st=+'], False)

        # Reverse strand alignment for SNP
        left_read = self.mk_read(query_name='chr1|48|C|A,T', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        right_read = self.mk_read(query_name='chr1|48|C|A,T', reference_name='chr2', reference_start=48, reference_end=108, is_reverse=True)
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='G'):
            assert calculate_new_variant_definition(left_read, right_read, fasta) == (48, 'G', ['T', 'A'], ['st=-'], False)

        # Forward strand alignment for SNP with novel allele
        left_read = self.mk_read(query_name='chr1|48|T|A', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(query_name='chr1|48|T|A', reference_name='chr2', reference_start=48, reference_end=108, is_reverse=False)
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='C'):
            assert calculate_new_variant_definition(left_read, right_read, fasta) == (48, 'C', ['A', 'T'], ['st=+', 'rac=T-C', 'nra'], True)

        # Forward strand alignment for Deletion
        left_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=50, reference_end=110, is_reverse=False)
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='CAA'):
            assert calculate_new_variant_definition(left_read, right_read, fasta) == (48, 'CAA', ['C'], ['st=+'], False)

        # Reverse strand alignment for a deletion
        # REF  AAAAAAAAAAAAAAAAATTGCCCCCCCCCCCCCCCCC
        #            read 2             read 1
        #      <----------------TTG<----------------
        #                       G
        left_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        right_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=50, reference_end=110, is_reverse=True)
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='TTG'):
            assert calculate_new_variant_definition(left_read, right_read, fasta) == (48, 'TTG', ['G'], ['st=-'], False)


    @staticmethod
    def get_test_resource(resource_name):
        """
        Gets full path to the test resource located in the same directory as the test module.
        """
        # Full path to this module.
        this_module = os.path.abspath(__file__)
        # Full path to the directory where it is contained.
        module_directory = os.path.dirname(this_module)
        # Full path to the directory where the test resources are (variant-remapping/tests/resources)
        resource_directory = os.path.dirname(os.path.dirname(module_directory)) + '/tests'
        # Full path to the requested resource.
        return os.path.join(resource_directory, 'resources', resource_name)
