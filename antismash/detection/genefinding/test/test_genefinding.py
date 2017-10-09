# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest
from argparse import Namespace

from antismash.detection.genefinding import check_options, is_enabled


class TestCore(unittest.TestCase):
    def test_check_options(self):
        options = Namespace()
        options.taxon = 'bacteria'
        options.genefinding_tool = "none"
        with self.assertRaises(AttributeError):
            check_options(options)
        options.genefinding_gff3 = '/nonexistant/path/to.gff'
        assert len(check_options(options)) == 1
        options.genefinding_gff3 = '/dev/null'
        assert not check_options(options)

    def test_is_enabled(self):
        options = Namespace()
        options.taxon = 'bacteria'
        options.genefinding_tool = 'none'
        options.genefinding_gff3 = False
        assert not is_enabled(options)

        options.genefinding_tool = 'prodigal'
        assert is_enabled(options)

        options.genefinding_gff3 = '/some/path'
        assert is_enabled(options)

        options.genefinding_tool = 'none'
        assert not is_enabled(options)
