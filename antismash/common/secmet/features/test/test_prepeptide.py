# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from antismash.common.secmet import Prepeptide
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation


class TestConversion(unittest.TestCase):
    def test_basic_conversion(self):
        old = Prepeptide(
            FeatureLocation(5, 95),
            peptide_class="test_class",
            core="coreseq...",
            locus_tag="loc",
            tool="test tool",
            peptide_subclass="test_subclass",
            score=20.4,
            monoisotopic_mass=6.7,
            molecular_weight=0.5,
            alternative_weights=[5.2, 6.7, 20.5],
            leader="leaderseq.",
            tail="tailseq..."
        )
        leader, core, tail = old.to_biopython()

        assert leader.location.start == 5
        assert leader.location.end == 35
        assert leader.qualifiers["prepeptide"] == ["leader"]

        assert core.location.start == 35
        assert core.location.end == 65
        assert core.qualifiers["prepeptide"] == ["core"]

        assert tail.location.start == 65
        assert tail.location.end == 95
        assert tail.qualifiers["prepeptide"] == ["tail"]

        with self.assertRaisesRegex(ValueError, "can only be reconstructed from core feature"):
            Prepeptide.from_biopython(leader)
        with self.assertRaisesRegex(ValueError, "can only be reconstructed from core feature"):
            Prepeptide.from_biopython(tail)

        new = Prepeptide.from_biopython(core)
        assert isinstance(new, Prepeptide)
        assert str(new.location) == str(old.location)
        assert new.peptide_class == old.peptide_class
        assert new.core == old.core
        assert new.locus_tag == old.locus_tag
        assert new.peptide_subclass == old.peptide_subclass
        assert new.score == old.score
        assert new.monoisotopic_mass == old.monoisotopic_mass
        assert new.molecular_weight == old.molecular_weight
        assert new.alternative_weights == old.alternative_weights
        assert new.leader == old.leader
        assert new.tail == old.tail

    def test_compound_location(self):
        old = Prepeptide(
            CompoundLocation([FeatureLocation(10, 50, 1),
                              FeatureLocation(130, 180, 1)],
                             operator="join"),
            peptide_class="test_class",
            core="coreseq...",
            locus_tag="loc",
            tool="test tool",
            leader="10chleader",
            tail="10chartail"
        )

        leader, core, tail = old.to_biopython()
        assert leader.location.start == 10
        assert leader.location.end == 40
        assert isinstance(core.location, CompoundLocation)
        assert core.location.start == 40
        assert core.location.end == 150
        assert tail.location.start == 150
        assert tail.location.end == 180

        new = Prepeptide.from_biopython(core)
        assert str(new.location) == str(old.location)
