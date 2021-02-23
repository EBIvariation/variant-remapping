import os
from unittest import TestCase

from variant_remapping_tools.replace_refs import process_variant


class TestReplaceRef(TestCase):

    def test_process_variant(self):
        # Reference allele from new genome was inserted in the vcf file and as a result the ref and alt are the same
        line = '\t'.join(['chr1', '1', '.', 'G', 'G', '.', '.', '.'])
        oldrefalleles_file = os.path.join(os.path.dirname(__file__), 'oldrefalleles.txt')
        line_num = 0
        # The old ref allele was reinserted
        assert process_variant(line, line_num, oldrefalleles_file) == 'chr1	1	.	C	G	.	.	.'
        # Note that I think this is wrong: the old ref allele should have been inserted in the ALT not in REF column

    def test_process_variant_deletion(self):
        # Reference allele from new genome was inserted in the file which make the deletion disappear
        line = '\t'.join(['chr2', '1', '.', 'T', 'C', '.', '.', '.'])
        oldrefalleles_file = os.path.join(os.path.dirname(__file__), 'oldrefalleles.txt')
        line_num = 1
        # The old ref allele was reinserted because the second line in oldref file is TAG
        assert process_variant(line, line_num, oldrefalleles_file) == 'chr2	1	.	TAG	C	.	.	.'

    def test_process_variant_no_change(self):
        # The bases are different so no changes is applied
        line = '\t'.join(['chr1', '1', '.', 'T', 'G', '.', '.', '.'])
        oldrefalleles_file = os.path.join(os.path.dirname(__file__), 'oldrefalleles.txt')
        line_num = 0
        assert process_variant(line, line_num, oldrefalleles_file) == 'chr1	1	.	T	G	.	.	.'
