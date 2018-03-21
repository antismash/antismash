# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from Bio import SearchIO
import unittest

from antismash.common import path

class TestBiopythonParseBug(unittest.TestCase):
    def test_searchio(self):
        hmm_results = open(path.get_full_path(__file__, 'data', 'biopython_bug.output'))
        try:
            results = list(SearchIO.parse(hmm_results, 'hmmer2-text'))
            self.fail("This isn't really a failure as much as it is that the biopython bug "
                      "keeping a try/except in subprocessing.run_hmmpfam() is no longer "
                      "around and the block can be removed")
        except ValueError as err:
            assert "Sequence lengths do not match" in str(err)
