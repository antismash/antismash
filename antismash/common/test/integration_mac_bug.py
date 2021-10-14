# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from unittest import mock, TestCase

from antismash import config
from antismash.common import record_processing
from antismash.common.test.helpers import get_path_to_nisin_fasta
from antismash.support import genefinding


class TestMacBug(TestCase):
    """ For python 3.8+ on mac, parallel behaviour of Pool() changed
        from using fork() to spawn(), which meant that config objects
        failed to serialise in the same way.

        Since the parallel section involved is in record_processing, the
        entry point is tested there.
    """

    def tearDown(self):
        config.destroy_config()

    def test_mac_bad_parallel(self):
        # single cpu is important, the mock itself fails to pickle under py3.7
        options = config.build_config(["--cpus", "1", "--genefinding-tool", "prodigal"],
                                      isolated=True, modules=[genefinding])
        config.update_config({
            "genefinding_gff": "",
            "genefinding_tool": "prodigal",
            "taxon": "bacteria",
        })

        original = record_processing.ensure_cds_info

        def wrapper(*args, **kwargs):
            config.destroy_config()  # fake the bug's effect
            return original(*args, **kwargs)

        records = record_processing.parse_input_sequence(get_path_to_nisin_fasta())
        assert len(records) == 1
        assert not records[0].get_cds_features()

        with mock.patch.object(record_processing, "ensure_cds_info", wraps=wrapper):
            record_processing.pre_process_sequences(records, options, genefinding)
        assert records[0].get_cds_features()
