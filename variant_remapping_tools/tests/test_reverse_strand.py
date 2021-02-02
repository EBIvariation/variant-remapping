from collections import Counter
from unittest import TestCase
from unittest.mock import Mock

import pysam
from pysam.libcalignedsegment import AlignedSegment

from variant_remapping_tools.reverse_strand import calculate_variant_position, is_read_valid


class TestReverseStrand(TestCase):

    def _mock_read(self, pos=1000, is_unmapped=False, is_secondary=False, cigartuples=((pysam.CMATCH, 101),), tags=None):
        if tags is None:
            tags = {'AS': -6}
        read = Mock(pos=pos, is_unmapped=is_unmapped, is_secondary=is_secondary, cigartuples=cigartuples)
        read.get_tag = lambda x: tags[x]
        return read

    def setUp(self) -> None:
        self.counter = Counter()
        self.score_cutoff = -20
        self.diff_cutoff = .06

    def test_calculate_variant_position(self):
        # 101M: Perfect alignment
        read = self._mock_read()
        assert calculate_variant_position(read, 50) == (1051, 10, 5)

    def test_calculate_variant_position_with_close_del(self):
        # 51M2D50M: Deletion next to the variant position
        read = self._mock_read(cigartuples=((pysam.CMATCH, 51), (pysam.CDEL, 2), (pysam.CMATCH, 50)))
        # The script considers this to be enough to filter at the moment
        assert calculate_variant_position(read, 50) == (1051, 9, 5)

    def test_calculate_variant_position_with_far_del(self):
        # 11M2D90M: Deletion far from the variant position
        read = self._mock_read(cigartuples=((pysam.CMATCH, 11), (pysam.CDEL, 2), (pysam.CMATCH, 90)))
        assert calculate_variant_position(read, 50) == (1053, 10, 5)

    def test_calculate_variant_position_with_close_ins(self):
        # 51M2I50M: Insertion next to the variant position
        read = self._mock_read(cigartuples=((pysam.CMATCH, 51), (pysam.CINS, 2), (pysam.CMATCH, 50)))
        # The script considers this to be enough to filter at the moment
        assert calculate_variant_position(read, 50) == (1051, 9, 5)

    def test_calculate_variant_position_with_soft_clip(self):
        read = self._mock_read(cigartuples=((pysam.CSOFT_CLIP, 10), (pysam.CMATCH, 91)))
        assert calculate_variant_position(read, 50) == (1041, 10, 5)

    def test_is_read_valid(self):
        read = self._mock_read()
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is True
        assert self.counter == Counter({'total': 1, 'remapped': 1})

    def test_is_read_invalid_unmapped(self):
        read = self._mock_read(is_unmapped=True)
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is False
        assert self.counter == Counter({'total': 1, 'unmapped': 1})

    def test_is_read_invalid_secondary(self):
        read = self._mock_read(is_secondary=True)
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is False
        assert self.counter == Counter()

    def test_is_read_invalid_AS(self):
        read = self._mock_read(tags={'AS': -21})
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is False
        assert self.counter == Counter({'total': 1, 'primary_poor': 1})

    def test_is_read_invalid_diff_AS_XS(self):
        read = self._mock_read(tags={'AS': -6, 'XS': -6})
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is False
        assert self.counter == Counter({'total': 1, 'gap_small': 1})

    def test_is_read_invalid_alignment_context(self):
        read = self._mock_read(cigartuples=((pysam.CMATCH, 51), (pysam.CDEL, 2), (pysam.CMATCH, 50)))
        assert is_read_valid(read, self.counter, 50, self.score_cutoff, self.diff_cutoff) is False
        assert self.counter == Counter({'total': 1, 'context_bad': 1})


