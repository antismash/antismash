# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import tempfile
import unittest
from unittest.mock import patch

from antismash.detection.subclusters import signatures
from antismash.detection.subclusters.signatures import SubclusterHmmSignature, get_signatures


def build_signature(name):
    return SubclusterHmmSignature(name, "desc", 10, f"{name}.hmm", accession="PF00001.1")


class TestReadSignatures(unittest.TestCase):
    def read_with_dir(self, text):
        """ Returns the signatures read from the given text, along with the
            directory the details file was in
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, "hmmdetails.txt")
            with open(filename, "w", encoding="utf-8") as handle:
                handle.write(text)
            return signatures._read_signatures(filename), temp_dir

    def read(self, text):
        return self.read_with_dir(text)[0]

    def test_simple(self):
        result = self.read("ECH_1\tEnoyl-CoA hydratase\t27\trule_hmms/ECH_1.hmm\tPF00378.26\n")
        assert len(result) == 1
        signature = result[0]
        assert signature.name == "ECH_1"
        assert signature.description == "Enoyl-CoA hydratase"
        assert signature.cutoff == 27
        assert signature.accession == "PF00378.26"

    def test_path_relative_to_details_file(self):
        line = "ECH_1\tdesc\t27\trule_hmms/ECH_1.hmm\tPF00378.26\n"
        result, details_dir = self.read_with_dir(line)
        assert result[0].hmm_file == os.path.join(details_dir, "rule_hmms", "ECH_1.hmm")

    def test_comments_and_blank_lines(self):
        result = self.read("# name\tdesc\tcutoff\tpath\taccession\n\n  \n"
                           "ECH_1\tdesc\t27\tECH_1.hmm\tPF00378.26\n")
        assert [signature.name for signature in result] == ["ECH_1"]

    def test_without_accession(self):
        for line in ["ECH_1\tdesc\t27\tECH_1.hmm", "ECH_1\tdesc\t27\tECH_1.hmm\t"]:
            signature = self.read(line + "\n")[0]
            assert signature.name == "ECH_1"
            assert signature.accession is None

    def test_wrong_column_count(self):
        for line in ["ECH_1\tdesc\t27",
                     "ECH_1\tdesc\t27\tECH_1.hmm\tPF00378.26\textra",
                     "ECH_1 desc 27 ECH_1.hmm PF00378.26"]:  # spaces instead of tabs
            with self.assertRaisesRegex(ValueError, "Invalid lines in HMM detail file"):
                self.read(line + "\n")

    def test_bad_lines_reported(self):
        good = "ECH_1\tdesc\t27\tECH_1.hmm\n"
        with self.assertRaisesRegex(ValueError, "first bad(.|\n)*second bad") as context:
            self.read(f"first bad\n{good}second bad\n")
        assert "ECH_1" not in str(context.exception)


class TestGetSignatures(unittest.TestCase):
    def setUp(self):
        cache_patch = patch.dict(signatures._SIGNATURE_CACHE, clear=True)
        cache_patch.start()
        self.addCleanup(cache_patch.stop)

    def test_loads_once(self):
        sigs = [build_signature("A"), build_signature("B")]
        with patch.object(signatures, "_read_signatures", return_value=sigs) as patched:
            first = get_signatures()
            assert get_signatures() is first
        patched.assert_called_once_with(signatures.DETAILS_FILE)
        assert first == {"A": sigs[0], "B": sigs[1]}

    def test_duplicate_names(self):
        sigs = [build_signature("A"), build_signature("A")]
        with patch.object(signatures, "_read_signatures", return_value=sigs):
            with self.assertRaisesRegex(ValueError, "Duplicate signature names"):
                get_signatures()
