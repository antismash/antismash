# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from unittest import TestCase

from antismash.config import get_config, update_config
from antismash.common.record_processing import parse_input_sequence, pre_process_sequences
from antismash.common.test.helpers import get_simple_options, get_path_to_nisin_fasta
from antismash.detection import genefinding


class TestProdigal(TestCase):
    def setUp(self):
        self.options = update_config(get_simple_options(genefinding,
                                                        ['--genefinding-tool', 'prodigal',
                                                         '--cpus', '1']))

    def tearDown(self):
        get_config().__dict__.clear()

    def test_nisin(self):
        record = parse_input_sequence(get_path_to_nisin_fasta())[0]
        assert record.get_feature_count() == 0
        record = pre_process_sequences([record], self.options, genefinding)[0]
        assert record.get_feature_count() == 12
        # and make sure they're all CDS features
        assert len(record.get_cds_features()) == 12
