# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.detection.genefunctions.tools.halogenases.data_structures import (
    HalogenaseHit,
)

from antismash.detection.genefunctions.tools.halogenases.flavin_dependent import (
    substrate_analysis,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent.substrate_analysis import (
    SUBGROUPS,
    categorise_on_substrate_level,
)

pyrrolic = SUBGROUPS["pyrrolic"]
PYRROLE = pyrrolic.profiles["pyrrole"]


class TestPyrrolic(unittest.TestCase):
    def setUp(self):
        self.profile = PYRROLE
        self.strong_pyrrole_hit = HalogenaseHit(
            query_id="query",
            bitscore=1000,
            evalue=1e-200,
            reference_id=PYRROLE.profile_name,
            profile=PYRROLE.filename
        )

        self.bad_pyrrole_hit = HalogenaseHit(
            query_id="query",
            bitscore=1,
            evalue=0.5,
            reference_id=PYRROLE.profile_name,
            profile=PYRROLE.filename
        )

    def test_score_below_cutoff(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={}) as patched:
            result = categorise_on_substrate_level("MAGIC", [self.bad_pyrrole_hit])
            patched.assert_called_once_with("MAGIC", self.bad_pyrrole_hit, self.profile.motifs)
        assert not result

    def test_no_motif(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues", return_value={}) as patched:
            result = categorise_on_substrate_level("MAGIC", [self.strong_pyrrole_hit])
            patched.assert_called_once_with("MAGIC", self.strong_pyrrole_hit, self.profile.motifs)
        assert not result

    def test_conventional_mono_di(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"mono_di": "DRSVFW"}) as patched:
            matches = categorise_on_substrate_level("MAGIC", [self.strong_pyrrole_hit])
            patched.assert_called_once_with("MAGIC", self.strong_pyrrole_hit, self.profile.motifs)
        match = matches[0]
        assert match.confidence == 1
        assert match.number_of_decorations == "mono_di"
        assert match.profile == PYRROLE.profile_name
        assert match.substrate == PYRROLE.name

    def test_unconventional_mono_di(self):
        name = "unconv_mono_di"
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={name: "YRRNFN"}) as patched:
            matches = categorise_on_substrate_level("MAGIC", [self.strong_pyrrole_hit])
            patched.assert_called_once_with("MAGIC", self.strong_pyrrole_hit, self.profile.motifs)
        match = matches[0]
        assert match.confidence == 1
        assert match.number_of_decorations == name
        assert match.profile == PYRROLE.profile_name
        assert match.substrate == PYRROLE.name

    def test_tetra(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues",
                          return_value={"tetra_mono_di": "RRYFFA"}) as patched:
            matches = categorise_on_substrate_level("MAGIC", [self.strong_pyrrole_hit])
            match = matches[0]
            assert match.confidence == 1
            assert match.number_of_decorations == "tetra"
            assert match.profile == PYRROLE.profile_name
            assert match.substrate == PYRROLE.name
            patched.assert_called_once_with("MAGIC", self.strong_pyrrole_hit, self.profile.motifs)

    def test_residue_mismatch(self):
        residues = {PYRROLE.motifs[0].name: "ABCDEF"}
        print(self.strong_pyrrole_hit)
        assert not PYRROLE.get_matches_from_hit(residues, self.strong_pyrrole_hit)
