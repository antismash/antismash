# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

import antismash
from antismash.common import record_processing
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import pfam2go

class PfamToGoTest(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--clusterhmmer", "--pfam2go", "--minimal"], isolated=True,
                                    modules=antismash.get_all_modules())
        # database needed since working with clusterhmmer?

    def tearDown(self):
        destroy_config()

    def check_add_to_record(self, input_file, results):
        record = record_processing.parse_input_sequence(input_file)[0]
        assert record.get_pfam_domains()
        for domain in record.get_pfam_domains():
            assert not domain.gene_ontologies
        results.add_to_record(record)
        # test it's been added to the record correctly

    def test_reuse(self):
        nisin = helpers.get_path_to_nisin_genbank()
        record = record_processing.parse_input_sequence(nisin)[0]

        results = helpers.run_and_regenerate_results_for_module(nisin, pfam2go, self.options)
        #  are the expected go ids for pfams found/no wrong ids for pfams? Only some samples first
        expected_pfams_and_gos_with_descs = {"PF00005": {"GO:0005524": "ATP binding", "GO:0016887": "ATPase activity"},
                                             "PF00072": {"GO:0000160": "phosphorelay signal transduction system"},
                                             "PF00486": {"GO:0003677": "DNA binding",
                                                         "GO:0000160": "phosphorelay signal transduction system",
                                                         "GO:0006355": "regulation of transcription, DNA-templated"}}
        expected_pfams_without_gos = ["PF05147", "PF04738"]
        pfams_found_with = []
        for pfam, all_ontologies in results.pfam_domains_with_gos.items():
            for pfam_id in pfam.db_xref:
                pfam_id_without_version = pfam_id.partition('.')[0]
                pfams_found_with.append(pfam_id_without_version)
                if pfam_id_without_version in expected_pfams_and_gos_with_descs:
                    for ontologies in all_ontologies:
                        if ontologies.pfam == pfam_id_without_version:
                            go_ids = [str(go_entry) for go_entry in ontologies.go_entries]
                            for go_id in go_ids:
                                assert go_id in expected_pfams_and_gos_with_descs[pfam_id_without_version]
        for expected_pfam, expected_gos in expected_pfams_and_gos_with_descs.items():
            assert expected_pfam in pfams_found_with
        for pfam_without_gos in expected_pfams_without_gos:
            assert pfam_without_gos not in pfams_found_with
