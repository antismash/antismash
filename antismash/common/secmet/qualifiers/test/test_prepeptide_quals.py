# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from ..prepeptide_qualifiers import (
    RiPPQualifier,
    LanthiQualifier,
    rebuild_qualifier,
)


class TestRebuilder(unittest.TestCase):
    def test_empty(self):
        assert rebuild_qualifier(None, "test") is None
        assert rebuild_qualifier({}, "test") is None

    def test_unknown(self):
        quals = {"RODEO_score": ["5"]}
        assert RiPPQualifier.from_biopython_qualifiers(quals.copy())
        with self.assertRaisesRegex(ValueError, "no known qualifier builder"):
            rebuild_qualifier(quals, "boop")


class TestRiPP(unittest.TestCase):
    def test_no_rodeo(self):
        qual = RiPPQualifier()
        assert qual.rodeo_score == 0
        assert qual.to_biopython_qualifiers() == {}

    def test_rodeo(self):
        qual = RiPPQualifier(rodeo_score=5)
        assert qual.rodeo_score == 5
        assert qual.to_biopython_qualifiers() == {"RODEO_score": ["5"]}

    def test_regeneration(self):
        old = RiPPQualifier(rodeo_score=3)
        new = RiPPQualifier.from_biopython_qualifiers(old.to_biopython_qualifiers())
        assert old.rodeo_score == new.rodeo_score


class TestLanthi(unittest.TestCase):
    def test_mods(self):
        def check_mods(mods, ami, chl, oxy, lac):
            assert ("AviCys" in mods) == ami
            assert ("Cl" in mods) == chl
            assert ("OH" in mods) == oxy
            assert ("Lac" in mods) == lac
        for amino_group in [True, False]:
            for chlorinated in [True, False]:
                for oxygenated in [True, False]:
                    for lactonated in [True, False]:
                        qual = LanthiQualifier(2, 5, amino_group, chlorinated,
                                               oxygenated, lactonated)
                        assert qual.aminovinyl_group == amino_group
                        assert qual.chlorinated == chlorinated
                        assert qual.oxygenated == oxygenated
                        assert qual.lactonated == lactonated

                        check_mods(qual.get_modifications(), amino_group, chlorinated,
                                   oxygenated, lactonated)

                        new = LanthiQualifier.from_biopython_qualifiers(qual.to_biopython_qualifiers())
                        check_mods(new.get_modifications(), amino_group, chlorinated,
                                   oxygenated, lactonated)

    def test_numerics(self):
        qual = LanthiQualifier(2, 5, False, False, False, False)
        assert qual.lan_bridges == 2
        assert qual.rodeo_score == 5

        new = LanthiQualifier.from_biopython_qualifiers(qual.to_biopython_qualifiers())
        assert new.lan_bridges == 2
        assert new.rodeo_score == 5

    def test_rebuilder(self):
        qual = LanthiQualifier(2, 5, False, True, False, True)
        bio = qual.to_biopython_qualifiers()
        assert bio
        bio["unrelated"] = ["qualifiers"]
        new = rebuild_qualifier(bio, "lanthipeptide")
        assert list(bio) == ["unrelated"]  # related ones should have been consumed
        assert isinstance(new, LanthiQualifier)
        assert new.rodeo_score == 5
        assert new.lan_bridges == 2
