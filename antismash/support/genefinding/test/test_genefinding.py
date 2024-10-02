# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from argparse import Namespace
import unittest
from unittest.mock import patch

from antismash.common.errors import AntismashInputError
from antismash.support import genefinding
from antismash.support.genefinding import check_options, is_enabled, run_on_record


class TestCore(unittest.TestCase):
    def test_check_options(self):
        options = Namespace()
        options.taxon = 'bacteria'
        options.genefinding_tool = "none"
        assert not check_options(options)

        options.taxon = 'fungi'
        assert not check_options(options)

        options.genefinding_tool = "prodigal"
        assert check_options(options)

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

    def test_run_on_record_error_cases(self):
        class FakeRecord:
            id = "fake"
        record = FakeRecord()
        options = Namespace()
        options.genefinding_tool = "none"

        for taxon in ("fungi", "bacteria"):
            options.taxon = taxon
            assert run_on_record(record, options) is None

        options.genefinding_tool = "error"
        for taxon in ("fungi", "bacteria"):
            options.taxon = taxon
            with self.assertRaisesRegex(AntismashInputError, "contains no genes and"):
                run_on_record(record, options)

        with patch.object(genefinding, "run_prodigal") as patched:
            options.taxon = "bacteria"
            options.genefinding_tool = "prodigal"
            run_on_record(record, options)
            patched.assert_called_once()

        with self.assertRaisesRegex(ValueError, "Unknown genefinding tool"):
            options.taxon = "bacteria"
            options.genefinding_tool = "bob"
            run_on_record(record, options)
