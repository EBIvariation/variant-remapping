import os
from collections import Counter
from unittest.mock import Mock, patch

import pysam
from unittest import TestCase


from variant_remapping_tools.reads_to_remapped_variants import fetch_bases, process_bam_file, \
    calculate_new_variants_definition, order_reads, link_supplementary, pass_aligned_filtering, \
    update_vcf_record


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
        process_bam_file([bamfile], output_file, out_failed_file, fasta_path, True, 50, summary_file)

        expected = [
            'chr2	98	.	C	CG	50	PASS	st=+	GT:GQ	1/1:0\n',
            'chr2	1078	.	A	G	50	PASS	st=+;rac=chr2|1078|G-A	GT	0/0\n',
            'chr2	1818	.	AAC	A	50	PASS	st=+	GT:GQ	1/1:0\n',
            'chr2	2030	.	A	TCC	50	PASS	st=+	GT:GQ	1/1:0\n'
        ]
        with open(output_file, 'r') as vcf:
            for(i, line) in enumerate(vcf):
                assert line == expected[i]
            # Expected and Generated VCF have the same number of lines
            assert i+1 == len(expected)

    def test_pass_aligned_filtering(self):
        left_read = self.mk_read(reference_start=0, reference_end=47, cigartuples=[(pysam.CMATCH, 47)])
        right_read = self.mk_read(
            reference_start=51, reference_end=95,  cigartuples=[(pysam.CSOFT_CLIP, 3), (pysam.CMATCH, 44)]
        )
        counter = Counter()
        # Soft clipped alignment on the side of the variant are filtered out
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
        left_read, right_read = order_reads(primary_group, primary_to_supplementary)
        assert left_read == primary_group[0]
        assert right_read == primary_group[1]

        primary_group = [
            self.mk_read(reference_name='chr2', reference_start=48, reference_end=58, is_reverse=True),
            self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        ]
        left_read, right_read = order_reads(primary_group, primary_to_supplementary)
        assert left_read == primary_group[1]
        assert right_read == primary_group[0]

        primary_group = [
            self.mk_read(reference_name='chr2', reference_start=1, reference_end=10, is_reverse=False),
            self.mk_read(reference_name='chr2', reference_start=48, reference_end=58, is_reverse=False)
        ]
        supplementary_read = self.mk_read(reference_name='chr2', reference_start=20, reference_end=47, is_reverse=False)
        primary_to_supplementary = {primary_group[0]: [supplementary_read]}
        left_read, right_read = order_reads(primary_group, primary_to_supplementary)
        assert left_read == supplementary_read
        assert right_read == primary_group[1]

    def test_calculate_new_variants_definition(self):
        fasta = 'fasta_path'

        # Forward strand alignment for SNP
        left_read = self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(reference_name='chr2', reference_start=48, reference_end=108, is_reverse=False)
        vcf_rec = ['chr1', '48', '.', 'C', 'A']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='C'):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (48, 'C', ['A'], {'st': '+'}, None)

        # Reverse strand alignment for SNP
        left_read = self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        right_read = self.mk_read(reference_name='chr2', reference_start=48, reference_end=108, is_reverse=True)
        vcf_rec = ['chr1', '48', '.', 'C', 'A,T']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='G'):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (48, 'G', ['T', 'A'], {'st': '-'}, None)

        # Forward strand alignment for SNP with novel allele
        left_read = self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(reference_name='chr2', reference_start=48, reference_end=108, is_reverse=False)
        vcf_rec = ['chr1', '48', '.', 'T', 'A']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='C'):
            var_generator = calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)
            assert next(var_generator) == \
                   (48, 'C', ['A'], {'st': '+', 'rac': 'chr2|48|T-C'}, None)
            assert next(var_generator) == \
                   (48, 'C', ['T'], {'st': '+', 'rac': 'chr2|48|T-C', 'nra': 'T'}, None)

        # Forward strand alignment for Deletion
        left_read = self.mk_read(reference_name='chr2', reference_start=1, reference_end=47, is_reverse=False)
        right_read = self.mk_read(reference_name='chr2', reference_start=50, reference_end=110, is_reverse=False)
        vcf_rec = ['chr1', '48', '.', 'CAA', 'C']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='CAA'):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (48, 'CAA', ['C'], {'st': '+'}, None)

        # Reverse strand alignment for a deletion
        # REF  AAAAAAAAAAAAAAAAATTGCCCCCCCCCCCCCCCCC
        #            read 2             read 1
        #      <----------------TTG<----------------
        #                       G
        left_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=1, reference_end=47, is_reverse=True)
        right_read = self.mk_read(query_name='chr1|48|CAA|C', reference_name='chr2', reference_start=50, reference_end=110, is_reverse=True)
        vcf_rec = ['chr1', '48', '.', 'CAA', 'C']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', side_effect=['TTG', 'ATT']):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (47, 'ATT', ['A'], {'st': '-'}, None)

    def test_calculate_new_variant_definition2(self):
        fasta = 'fasta_path'
        #  reverse strand alignment for Deletion with novel allele
        # variant:  CTGTG -> C
        # OLD REF  UUUUUUUUUUUUUUUUUCTGTGDDDDDDDDDDDDDDD
        # OLD ALT  UUUUUUUUUUUUUUUUUC----DDDDDDDDDDDDDDD

        #                read 2              read 1
        #          <----------------CACAG<--------------
        # NEW REF  DDDDDDDDDDDDDDDDACATAGUUUUUUUUUUUUUUU
        # NEW ALT  DDDDDDDDDDDDDDDDA----GUUUUUUUUUUUUUUU
        # variant:  ACATA -> A
        left_read = self.mk_read(query_name='chr1|48|CTGTG|C', reference_name='chr2', reference_start=1, reference_end=46, is_reverse=True)
        right_read = self.mk_read(query_name='chr1|48|CTGTG|C', reference_name='chr2', reference_start=50, reference_end=110, is_reverse=True)
        vcf_rec = ['chr1', '48', '.', 'CTGTG', 'C']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', side_effect=['CATAG', 'ACATA']):
            var_generator = calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)

            assert next(var_generator) == \
                   (46, 'ACATA', ['A'], {'st': '-', 'rac': 'chr2|46|ACACA-ACATA'}, None)
            assert next(var_generator) == \
                   (46, 'ACATA', ['ACACA'], {'st': '-', 'rac': 'chr2|46|ACACA-ACATA', 'nra': 'ACACA'}, None)

        #  reverse strand alignment for Insertion
        # For insertion there cannot be a reference allele change on the negative strand because the reference allele
        # is only the contex base
        # variant:  C -> CTGTG
        # OLD REF  UUUUUUUUUUUUUUUUUC----DDDDDDDDDDDDDDD
        # OLD ALT  UUUUUUUUUUUUUUUUUCTGTGDDDDDDDDDDDDDDD

        #                read 2              read 1
        #          <----------------****G<--------------
        # NEW REF  DDDDDDDDDDDDDDDDA----GUUUUUUUUUUUUUUU
        # NEW ALT  DDDDDDDDDDDDDDDDACACAGUUUUUUUUUUUUUUU
        # variant:  A -> ACACA
        left_read = self.mk_read(query_name='chr1|48|C|CTGTG', reference_name='chr2', reference_start=1,
                                 reference_end=46, is_reverse=True)
        right_read = self.mk_read(query_name='chr1|48|C|CTGTG', reference_name='chr2', reference_start=50,
                                  reference_end=110, is_reverse=True)
        vcf_rec = ['chr1', '48', '.', 'C', 'CTGTG']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases',
                   side_effect=['G', 'A']):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (46, 'A', ['ACACA'], {'st': '-'}, None)

        # Forward strand alignment for Insertion
        # Reference context base changes should not create reference allele change because it's not part of the variant.
        # variant:  C -> CTGTG
        # OLD REF  UUUUUUUUUUUUUUUUUC----DDDDDDDDDDDDDDD
        # OLD ALT  UUUUUUUUUUUUUUUUUCTGTGDDDDDDDDDDDDDDD

        #                read 1              read 2
        #          ---------------->****<--------------
        # NEW REF  UUUUUUUUUUUUUUUUUT----DDDDDDDDDDDDDDD
        # NEW ALT  UUUUUUUUUUUUUUUUUTTGTGDDDDDDDDDDDDDDD
        # variant:  T -> TTGTG
        left_read = self.mk_read(query_name='chr1|48|C|CTGTG', reference_name='chr2', reference_start=1,
                                 reference_end=46, is_reverse=False)
        right_read = self.mk_read(query_name='chr1|48|C|CTGTG', reference_name='chr2', reference_start=50,
                                  reference_end=110, is_reverse=False)
        vcf_rec = ['chr1', '48', '.', 'C', 'CTGTG']
        with patch('variant_remapping_tools.reads_to_remapped_variants.fetch_bases', return_value='T'):
            assert next(calculate_new_variants_definition(left_read, right_read, fasta, vcf_rec)) == \
                   (47, 'T', ['TTGTG'], {'st': '+'}, None)

    def test_update_vcf_record(self):
        # Allele swap no genotype change
        original_vcf_rec = ['chr1', '10', '.', 'A', 'T', '50', '.', '.', 'GT', '0/1']
        update_vcf_record('1', 11, 'T', 'A', {'rac': 'A-T'}, original_vcf_rec)
        assert original_vcf_rec == ['1', '11', '.', 'T', 'A', '50', '.', 'rac=A-T', 'GT', '0/1']

        # Allele swap with genotype change
        original_vcf_rec = ['chr1', '10', '.', 'A', 'T', '50', '.', '.', 'GT', '1/1']
        update_vcf_record('1', 11, 'T', 'A', {'rac': 'A-T'}, original_vcf_rec)
        assert original_vcf_rec == ['1', '11', '.', 'T', 'A', '50', '.', 'rac=A-T', 'GT', '0/0']

        # Novel reference allele with genotype change
        original_vcf_rec = ['chr1', '10', '.', 'A', 'T', '50', '.', '.', 'GT', '0/1']
        update_vcf_record('1', 11, 'C', ['A', 'T'], {'rac': 'A-C', 'nra': 'A'}, original_vcf_rec)
        assert original_vcf_rec == ['1', '11', '.', 'C', 'A,T', '50', '.', 'rac=A-C;nra=A', 'GT', '1/2']


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
