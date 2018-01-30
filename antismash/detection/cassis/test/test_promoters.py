# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Tests for CASSIS promoters. """

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from minimock import mock, restore

from antismash.common import secmet
from antismash.common.test import helpers
from antismash.detection.cassis import promoters as cassis


class TestGetPromoters(unittest.TestCase):
    """ This class tests each section of get_promoters in turn. Each section has an
        explicit test to ensure its coverage (hence the mocking of enumerate)."""
    def setUp(self):
        self.record = helpers.DummyRecord(seq=Seq("ACGT"*25))

    def tearDown(self):
        restore()

    def add_gene(self, name, start, end, strand=1):
        gene = secmet.feature.Gene(FeatureLocation(start, end, strand), locus_tag=name)
        self.record.add_gene(gene)
        return gene

    def get_promoters(self, upstream, downstream):
        return cassis.get_promoters(self.record, self.record.get_genes(), upstream, downstream)

    def check_single_promoter(self, promoter, name, start, end):
        print(promoter)
        assert isinstance(promoter, cassis.Promoter)
        assert not isinstance(promoter, cassis.CombinedPromoter)
        assert name == promoter.gene_name
        assert start == promoter.start
        assert end == promoter.end

    def check_combined_promoter(self, promoter, first_gene, second_gene, start, end):
        print(promoter)
        assert isinstance(promoter, cassis.CombinedPromoter)
        assert promoter.get_gene_names() == [first_gene, second_gene]
        assert promoter.get_id() == "+".join(promoter.get_gene_names())
        assert promoter.start == start
        assert promoter.end == end

    def test_single_gene_forward(self):
        self.add_gene("A", 10, 20, 1)
        assert len(self.record.get_genes()) == 1

        promoters = self.get_promoters(15, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 0, 15)

        promoters = self.get_promoters(15, 25)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 0, 20)

        promoters = self.get_promoters(5, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 5, 15)

        promoters = self.get_promoters(5, 25)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 5, 20)

    def test_single_gene_reverse(self):
        self.add_gene("B", 10, 80, -1)
        assert len(self.record.get_genes()) == 1

        promoters = self.get_promoters(5, 75)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "B", 10, 85)

        promoters = self.get_promoters(25, 75)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "B", 10, 99)

        promoters = self.get_promoters(5, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "B", 75, 85)

        promoters = self.get_promoters(25, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "B", 75, 99)

    def test_first_gene_forward(self):
        # ensure coverage only considers this gene of interest
        gene_of_interest = self.add_gene("A", 10, 20, 1)
        mock("enumerate", returns=[(0, gene_of_interest)])
        other_gene = self.add_gene("B", 30, 40, 1)

        for strand in [1, -1]:
            start, end = sorted([30, 40])[::strand]
            other_gene.location = FeatureLocation(start, end, strand)
            print(other_gene.location)

            promoters = self.get_promoters(5, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "A", 5, 20)

            promoters = self.get_promoters(25, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "A", 0, 20)

            promoters = self.get_promoters(5, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "A", 5, 15)

            promoters = self.get_promoters(25, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "A", 0, 15)

    def test_first_gene_reverse(self):
        gene_of_interest = self.add_gene("A", 20, 10, -1)
        mock("enumerate", returns=[(0, gene_of_interest)])
        other_gene = self.add_gene("B", 30, 40, 1)

        # counts as a shared promoter if strands are opposite and ranges are large
        # so a small test of the opposite strand
        promoters = self.get_promoters(2, 4)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 16, 22)

        # then reverse the gene and do more thorough testing
        other_gene.location = FeatureLocation(40, 30, -1)
        promoters = self.get_promoters(5, 75)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 10, 25)

        promoters = self.get_promoters(5, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 15, 25)

        promoters = self.get_promoters(25, 5)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 15, 29)

        promoters = self.get_promoters(25, 75)
        assert len(promoters) == 1
        self.check_single_promoter(promoters[0], "A", 10, 29)

    def test_last_gene_forward(self):
        other_gene = self.add_gene("A", 10, 20, 1)
        # ensure coverage only considers this gene of interest
        gene_of_interest = self.add_gene("B", 30, 40, 1)
        mock("enumerate", returns=[(1, gene_of_interest)])

        for strand in [1, -1]:
            start, end = sorted([10, 20])[::strand]
            other_gene.location = FeatureLocation(start, end, strand)
            print(other_gene.location)

            promoters = self.get_promoters(5, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 25, 40)

            promoters = self.get_promoters(25, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 21, 40)

            promoters = self.get_promoters(5, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 25, 35)

            promoters = self.get_promoters(25, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 21, 35)

    def test_last_gene_reverse(self):
        other_gene = self.add_gene("A", 10, 20, 1)
        # ensure coverage only considers this gene of interest
        gene_of_interest = self.add_gene("B", 40, 30, -1)
        mock("enumerate", returns=[(1, gene_of_interest)])

        for strand in [1, -1]:
            start, end = sorted([10, 20])[::strand]
            other_gene.location = FeatureLocation(start, end, strand)
            print(other_gene.location)

            promoters = self.get_promoters(5, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 30, 45)

            promoters = self.get_promoters(65, 75)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 30, 99)

            promoters = self.get_promoters(5, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 35, 45)

            promoters = self.get_promoters(65, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 35, 99)

    def test_special_case(self):
        first = self.add_gene("A", 20, 10, -1)
        second = self.add_gene("B", 30, 50, 1)

        promoters = self.get_promoters(5, 5)
        assert len(promoters) == 1
        assert second.location.end > second.location.start + 5
        assert first.location.start > first.location.end + 5
        self.check_combined_promoter(promoters[0], "A", "B", 15, 35)

        promoters = self.get_promoters(5, 15)
        assert len(promoters) == 1
        assert not first.location.start > first.location.end + 15
        assert second.location.end > second.location.start + 15
        self.check_combined_promoter(promoters[0], "A", "B", 10, 45)

        first.location = FeatureLocation(35, 10, -1)
        second.location = FeatureLocation(45, 60, 1)

        promoters = self.get_promoters(25, 20)
        assert len(promoters) == 1
        assert first.location.start > first.location.end + 20
        assert not second.location.end > second.location.start + 20
        self.check_combined_promoter(promoters[0], "A", "B", 15, 60)

        promoters = self.get_promoters(25, 30)
        assert not second.location.end > second.location.start + 30
        assert not first.location.start > first.location.end + 30
        assert len(promoters) == 1
        self.check_combined_promoter(promoters[0], "A", "B", 10, 60)

    def test_normal_case_forward(self):
        other = self.add_gene("A", 10, 20, 1)
        gene_of_interest = self.add_gene("B", 40, 60, 1)
        self.add_gene("C", 70, 80, 1)
        mock("enumerate", returns=[(1, gene_of_interest)])

        for strand in [-1, 1]:
            start, end = sorted([other.location.start, other.location.end])[::strand]
            other.location = FeatureLocation(start, end, strand)

            promoters = self.get_promoters(5, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 35, 45)

            promoters = self.get_promoters(5, 25)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 35, 60)

            promoters = self.get_promoters(25, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 21, 45)

            promoters = self.get_promoters(25, 25)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 21, 60)

    def test_normal_case_reverse(self):
        self.add_gene("A", 10, 20, 1)
        gene_of_interest = self.add_gene("B", 40, 60, -1)
        other = self.add_gene("C", 80, 70, -1)
        mock("enumerate", returns=[(1, gene_of_interest)])

        for strand in [-1]:
            start, end = sorted([other.location.start, other.location.end])[::strand]
            other.location = FeatureLocation(start, end, strand)

            promoters = self.get_promoters(5, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 55, 65)

            promoters = self.get_promoters(5, 25)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 40, 65)

            promoters = self.get_promoters(25, 5)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 55, 69)

            promoters = self.get_promoters(25, 25)
            assert len(promoters) == 1
            self.check_single_promoter(promoters[0], "B", 40, 69)

    def test_bad_strand(self):
        self.add_gene("A", 10, 20, 0)
        with self.assertRaisesRegex(ValueError, "Gene A has unknown strand: 0"):
            self.get_promoters(5, 5)
