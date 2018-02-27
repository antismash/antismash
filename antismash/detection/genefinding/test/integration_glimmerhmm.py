# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
from unittest import TestCase

from antismash.config import get_config, update_config
from antismash.common.record_processing import parse_input_sequence, pre_process_sequences
from antismash.common.path import get_full_path
from antismash.common.test.helpers import get_simple_options
from antismash.detection import genefinding


class TestGlimmerHMM(TestCase):
    def setUp(self):
        self.options = update_config(get_simple_options(genefinding,
                                     ['--taxon', 'fungi',
                                      '--genefinding-tool', 'glimmerhmm',
                                      '--cpus', '1']))
        self.data_location = get_full_path(__file__, "data")

    def tearDown(self):
        get_config().__dict__.clear()

    def data_file(self, filename):
        return os.path.join(self.data_location, filename)

    def test_fumigatus_cluster(self):
        record = parse_input_sequence(self.data_file('fumigatus.cluster1.fna'), taxon="fungi")[0]
        assert record.get_feature_count() == 0
        record = pre_process_sequences([record], self.options, genefinding)[0]
        assert record.get_feature_count() == 11
        # and make sure they're all CDS features
        assert len(record.get_cds_features()) == 11
