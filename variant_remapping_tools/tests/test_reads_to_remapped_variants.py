from unittest import TestCase
from variant_remapping_tools.reads_to_remapped_variants import fetch_bases, calculate_new_alleles, process_bam_file


class TestProcess(TestCase):

    def test_fetch_old_bases(self):
        fasta_path = '../../tests/resources/genome.fa'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=2, start=1) == 'CA'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=1, start=48) == 'C'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=1, start=98) == 'C'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=3, start=353) == 'GGA'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=1, start=1078) == 'G'
        assert fetch_bases(fasta=fasta_path, contig='chr1', length=1, start=1404) == 'G'

    def test_fetch_bases(self):
        fasta_path = '../../tests/resources/new_genome.fa'
        assert fetch_bases(fasta=fasta_path, contig='chr2', length=2, start=1) == 'CA'
        assert fetch_bases(fasta=fasta_path, contig='chr2', length=1, start=98) == 'C'
        assert fetch_bases(fasta=fasta_path, contig='chr2', length=1, start=353) == 'G'
        assert fetch_bases(fasta=fasta_path, contig='chr2', length=1, start=1078) == 'A'
        assert fetch_bases(fasta=fasta_path, contig='chr2', length=1, start=1404) == 'G'

    def test_compare_old_vs_new(self):
        # No changes
        assert calculate_new_alleles(old_ref='A', new_ref='A', old_alt='C', is_reverse_strand=False) == ('A', 'C')
        # Alignment reverse strand
        assert calculate_new_alleles(old_ref='A', new_ref='T', old_alt='C', is_reverse_strand=True) == ('T', 'G')
        # Allele changes
        assert calculate_new_alleles(old_ref='A', new_ref='C', old_alt='C', is_reverse_strand=False) == ('C', 'A')
        # Alignment reverse strand and allele change
        assert calculate_new_alleles(old_ref='A', new_ref='G', old_alt='C', is_reverse_strand=True) == ('G', 'T')

    def test_process_bam_file(self):
        self.skipTest("INDELS are not working yet")

        """
        Given:
        chr1	48	.	C	A	50	PASS	.	GT	1/1
        chr1	98	.	C	CG	50	PASS	.	GT	1/1
        chr1	353	.	GGA	G	50	PASS	.	GT	1/1
        chr1	1078	.	G	A	50	PASS	.	GT	1/1
        chr1	1404	.	G	GCG	50	PASS	.	GT	1/1
        """

        bamfile = "../../tests/resources/reads_aligned.sorted.bam"
        fasta_path = '../../tests/resources/new_genome.fa'
        output_file = '/tmp/remapped.vcf'
        process_bam_file(bamfile, output_file, 50, 0.6, 0.04, fasta_path)

        """
        Expected:
        chr2	98	.	C	CG	50	PASS	.
        chr2	353	.	G	GGA	50	PASS	.
        chr2	1078	.	A	G	50	PASS	.
        chr2	1404	.	G	GCG	50	PASS	.        
        """

        expected = ["chr2	98	.	C	CG	50	PASS	.\n",
                    "chr2	353	.	G	GGA	50	PASS	.\n",
                    "chr2	1078	.	A	G	50	PASS	.\n",
                    "chr2	1404	.	G	GCG	50	PASS	.\n"]

        with open(output_file, 'r') as vcf:
            for(i, line) in enumerate(vcf):
                assert line == expected[i]
