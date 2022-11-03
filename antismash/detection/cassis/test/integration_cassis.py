# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

import os
from tempfile import TemporaryDirectory
import unittest

from antismash.common import path
from antismash.config import destroy_config, build_config
from antismash.detection import cassis

from .test_cassis import create_fake_record


class TestCassisMainMethod(unittest.TestCase):
    def setUp(self):
        # because the fake record has a very repetitive sequence, bump the max percentage to 100
        self.old_max_perc = cassis.MAX_PERCENTAGE
        cassis.MAX_PERCENTAGE = 100.
        self.tempdir = TemporaryDirectory(prefix="as_cassis")
        self.options = build_config(["--cpus", "2", "--output-dir", self.tempdir.name],
                                    isolated=True, modules=[cassis])

    def tearDown(self):
        cassis.MAX_PERCENTAGE = self.old_max_perc
        destroy_config()
        self.tempdir.cleanup()
        assert not os.path.exists(path.get_full_path(__file__, "data", "meme"))
        assert not os.path.exists(path.get_full_path(__file__, "data", "fimo"))
        assert not os.path.isfile(path.get_full_path(__file__, "data", "test_promoter_positions.csv"))
        assert not os.path.isfile(path.get_full_path(__file__, "data", "test_promoter_sequences.fasta"))

    def test_detect(self):
        assert cassis.MAX_PERCENTAGE == 100.0, "if lower, detect stops early"
        seq_record = create_fake_record()
        cassis.detect(seq_record, self.options)
        # TODO it would be cool to capture the output (logging messages) of cassis.detect
        # one could than check if, e.g. "Best prediction (most abundant): 'gene1' -- 'gene9'" is in the output
        # but, if using nosetests, it cannot be captured, because nosetests already collects all output
