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
        results.add_to_record(record)
        for domain in record.get_pfam_domains():
            if domain.gene_ontologies:
                assert sorted(domain.gene_ontologies.ids) == sorted([go_entry.id
                                                               for ontologies in results.pfam_domains_with_gos[domain]
                                                               for go_entry in ontologies.go_entries])

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
        expected_pfams_found = set()
        for pfam, all_ontologies in results.pfam_domains_with_gos.items():
            pfam_ids_without_versions = [pfam_id.partition(".")[0] for pfam_id in pfam.db_xref]
            # make sure the Pfams without gos aren't in the results
            assert "PF05147" not in pfam_ids_without_versions and "PF04738" not in pfam_ids_without_versions
            for ontologies in all_ontologies:
                # make sure GeneOntologies' pfam id actually is one found in the domain's ids
                assert ontologies.pfam in pfam_ids_without_versions
                # did it find the right amount of GO IDs for the sample Pfams, and did it find the right ones?
                if ontologies.pfam in expected_pfams_and_gos_with_descs:
                    expected_pfams_found.add(ontologies.pfam)
                    go_ids = [str(go_entry) for go_entry in ontologies.go_entries]
                    assert len(go_ids) == len(expected_pfams_and_gos_with_descs[ontologies.pfam])
                    for go_id in go_ids:
                        assert go_id in expected_pfams_and_gos_with_descs[ontologies.pfam]
        # make sure all expected pfams have been found
        assert len(expected_pfams_found) == len(expected_pfams_and_gos_with_descs)

        self.check_add_to_record(nisin, results)
