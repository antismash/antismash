# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import tempfile
import unittest
from unittest.mock import patch

from antismash.detection.subclusters import compounds
from antismash.detection.subclusters.compounds import CompoundInfo, get_compound


class TestReadCompounds(unittest.TestCase):
    def read(self, text):
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".txt") as handle:
            handle.write(text)
            handle.flush()
            return compounds._read_compounds(handle.name)

    def test_simple(self):
        result = self.read("SCG0001\tname\tCCO\tamino acid\n")
        assert result == {"SCG0001": CompoundInfo(name="name", smiles="CCO",
                                                  classification=["amino acid"])}

    def test_comments_and_blank_lines(self):
        result = self.read("# rule\tname\tsmiles\tclasses\n\n   \nSCG0001\tname\tCCO\tsugar\n")
        assert list(result) == ["SCG0001"]

    def test_empty_smiles(self):
        result = self.read("SCG0001\tname\t\tsugar\n")
        assert result["SCG0001"].smiles is None

    def test_classification(self):
        result = self.read("SCG0001\tname\tCCO\t sugar ; amino acid;;\n")
        assert result["SCG0001"].classification == ["sugar", "amino acid"]

    def test_too_few_columns(self):
        with self.assertRaisesRegex(ValueError, "Invalid line"):
            self.read("SCG0001\tname\tCCO\n")


class TestGetCompound(unittest.TestCase):
    def setUp(self):
        cache_patch = patch.dict(compounds._COMPOUND_CACHE, clear=True)
        cache_patch.start()
        self.addCleanup(cache_patch.stop)
        self.compound = CompoundInfo(name="name", smiles=None, classification=[])

    def test_loads_once(self):
        with patch.object(compounds, "_read_compounds",
                          return_value={"SCG0001": self.compound}) as patched:
            assert get_compound("SCG0001") is self.compound
            assert get_compound("SCG0001") is self.compound
        patched.assert_called_once()

    def test_unknown_rule(self):
        with patch.object(compounds, "_read_compounds", return_value={"SCG0001": self.compound}):
            with self.assertRaises(KeyError):
                get_compound("SCG9999")
