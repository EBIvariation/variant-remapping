import os
import pysam
from unittest import TestCase
from variant_remapping_tools.reads_to_remapped_variants import fetch_bases, calculate_new_alleles, process_bam_file


class TestProcess(TestCase):

    def test_fetch_old_bases(self):
        fasta_path = self.get_test_resource('genome.fa')
        fasta = pysam.FastaFile(fasta_path)
        assert fetch_bases(fasta=fasta, contig='chr1', length=2, start=1) == 'CA'
        assert fetch_bases(fasta=fasta, contig='chr1', length=1, start=48) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr1', length=1, start=98) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr1', length=3, start=353) == 'GGA'
        assert fetch_bases(fasta=fasta, contig='chr1', length=1, start=1078) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr1', length=1, start=1404) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr1', length=3, start=1818) == 'AAC'
        assert fetch_bases(fasta=fasta, contig='chr1', length=1, start=2030) == 'A'

    def test_fetch_bases(self):
        fasta_path = self.get_test_resource('new_genome.fa')
        fasta = pysam.FastaFile(fasta_path)
        assert fetch_bases(fasta=fasta, contig='chr2', length=2, start=1) == 'CA'
        assert fetch_bases(fasta=fasta, contig='chr2', length=1, start=98) == 'C'
        assert fetch_bases(fasta=fasta, contig='chr2', length=1, start=353) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr2', length=1, start=1078) == 'A'
        assert fetch_bases(fasta=fasta, contig='chr2', length=1, start=1404) == 'G'
        assert fetch_bases(fasta=fasta, contig='chr2', length=3, start=1818) == 'AAC'
        assert fetch_bases(fasta=fasta, contig='chr2', length=1, start=2030) == 'A'

    def test_calculate_new_alleles(self):
        # No changes
        assert calculate_new_alleles(old_ref='A', new_ref='A', old_alt='C', is_reverse_strand=False) == ('A', 'C')
        # Alignment reverse strand
        assert calculate_new_alleles(old_ref='A', new_ref='T', old_alt='C', is_reverse_strand=True) == ('T', 'G')
        # Allele changes
        assert calculate_new_alleles(old_ref='A', new_ref='C', old_alt='C', is_reverse_strand=False) == ('C', 'A')
        # Alignment reverse strand and allele change
        assert calculate_new_alleles(old_ref='A', new_ref='G', old_alt='C', is_reverse_strand=True) == ('G', 'T')
        # Insertion with no allele changes
        assert calculate_new_alleles(old_ref='A', new_ref='A', old_alt='TCC', is_reverse_strand=False) == ('A', 'TCC')
        # Deletion with no allele changes
        assert calculate_new_alleles(old_ref='AAC', new_ref='AAC', old_alt='T', is_reverse_strand=False) == ('AAC', 'T')

    def test_process_bam_file(self):
        """
        Given:
        chr1	48	.	C	A	50	PASS	.	GT	1/1
        chr1	98	.	C	CG	50	PASS	.	GT	1/1
        chr1	1078	.	G	A	50	PASS	.	GT	1/1
        chr1	1818	.	AAC	T	50	PASS	.	GT	1/1
        chr1	2030	.	A	TCC	50	PASS	.	GT	1/1
        """

        bamfile = self.get_test_resource("reads_aligned_test.sorted.bam")
        fasta_path = self.get_test_resource('new_genome.fa')
        output_file = '/tmp/remapped.vcf'
        process_bam_file(bamfile, output_file, 50, 0.6, 0.04, 5, fasta_path)

        expected = ["chr2	98	.	C	CG	50	PASS	.\n",
                    "chr2	1078	.	A	G	50	PASS	.\n",
                    "chr2	1818	.	AAC	T	50	PASS	.\n",
                    "chr2	2030	.	A	TCC	50	PASS	.\n"]

        with open(output_file, 'r') as vcf:
            for(i, line) in enumerate(vcf):
                assert line == expected[i]
            # Expected and Generated VCF have the same number of lines
            assert i+1 == len(expected)

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
