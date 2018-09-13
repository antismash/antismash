# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from ..prepeptide_qualifiers import (
    RiPPQualifier,
    LanthiQualifier,
    LassoQualifier,
    ThioQualifier,
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
    def test_rodeo(self):
        for score in [-2, 0, 5]:
            qual = RiPPQualifier(rodeo_score=score)
            assert qual.rodeo_score == score
            bio = qual.to_biopython_qualifiers()
            assert bio == {"RODEO_score": [str(score)]}
            assert RiPPQualifier.from_biopython_qualifiers(bio).rodeo_score == score
            assert not bio

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


class TestThio(unittest.TestCase):
    def test_rebuilder(self):
        qual = ThioQualifier(5, True, "test macro", "core feats", [2.5, 5.7])
        bio = qual.to_biopython_qualifiers()
        assert bio
        bio["unrelated"] = ["qualifiers"]
        new = rebuild_qualifier(bio, "thiopeptide")
        assert isinstance(new, ThioQualifier)
        assert list(bio) == ["unrelated"]  # related ones should have been consumed
        assert new.rodeo_score == qual.rodeo_score
        assert new.macrocycle == qual.macrocycle
        assert new.core_features == qual.core_features
        assert new.mature_weights == qual.mature_weights
        assert new.amidation == qual.amidation

    def test_tail(self):
        qual = ThioQualifier(5, False, "test macro", "core feats", [2.5, 5.7])
        assert qual.tail_reaction == ""
        assert "tail_reaction" not in qual.to_biopython_qualifiers()
        qual.amidation = True
        assert qual.amidation
        assert qual.tail_reaction == "dealkylation of C-Terminal residue; amidation"
        bio = qual.to_biopython_qualifiers()
        assert bio["tail_reaction"] == [qual.tail_reaction]


class TestLasso(unittest.TestCase):
    def test_rebuilder(self):
        def check_values(qual):
            assert isinstance(qual, LassoQualifier)
            assert qual.rodeo_score == 7
            assert qual.num_bridges == 2
            assert qual.macrolactam == "test macro"
            assert qual.cut_mass == 5.6
            assert qual.cut_weight == 8.3

        old = LassoQualifier(7, 2, "test macro", 5.6, 8.3)
        check_values(old)
        bio = old.to_biopython_qualifiers()
        bio["unrelated"] = ["qualifiers"]
        new = LassoQualifier.from_biopython_qualifiers(bio)
        assert list(bio) == ["unrelated"]  # related ones should have been consumed
        check_values(new)

        bio = old.to_biopython_qualifiers()
        bio["unrelated"] = ["qualifiers"]
        new = rebuild_qualifier(bio, "lassopeptide")
        assert list(bio) == ["unrelated"]
        check_values(new)
