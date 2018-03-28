# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from typing import Dict, List
from antismash.common import path
from antismash.common.secmet.feature import PFAMDomain
from antismash.common.secmet.record import Record
from antismash.common.test.helpers import DummyRecord
from antismash.modules.pfam2go import pfam2go


class PfamToGoTest(unittest.TestCase):
    known_connections = {'PF00015': ['GO:0004871', 'GO:0007165', 'GO:0016020'],
                         'PF00351': ['GO:0016714', 'GO:0055114'],
                         'PF02364': ['GO:0003843', 'GO:0006075', 'GO:0000148', 'GO:0016020']}
    working_descs = {'GO:0004871': 'signal transducer activity', 'GO:0007165': 'signal transduction',
                     'GO:0016020': 'membrane',
                     'GO:0016714': 'oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced pteridine as one donor, and incorporation of one atom of oxygen',
                     'GO:0055114': 'oxidation-reduction process', 'GO:0003843': '1,3-beta-D-glucan synthase activity',
                     'GO:0006075': '(1->3)-beta-D-glucan biosynthetic process',
                     'GO:0000148': '1,3-beta-D-glucan synthase complex'}

    def set_dummy_with_pfams(self, pfam_ids: Dict[str, FeatureLocation]) -> DummyRecord:  # Dict id: FeatureLocation
        pfam_domains = []
        for pfam_id in pfam_ids:
            pfam_domain = PFAMDomain(location=pfam_ids[pfam_id], description='FAKE')
            pfam_domain.db_xref = [pfam_id]
            pfam_domains.append(pfam_domain)
        fake_record = DummyRecord(features=pfam_domains)
        return fake_record


    def test_gene_ontologies(self):
        # does it use arguments given? How is bad input handled?
        sample_ontology = pfam2go.GeneOntology('GO:0004871', 'signal transducer activity')
        sample_pfam = 'PF00015'
        sample_ontologies = pfam2go.GeneOntologies(sample_pfam, [sample_ontology])
        assert sample_ontologies.pfam == sample_pfam
        assert sample_ontologies.go_entries == [sample_ontology]
        all_entries = [str(go_entry) for go_entry in sample_ontologies.go_entries]
        assert sample_ontology.id in all_entries

    def test_gene_ontologies_fail(self):
        fail_ontology = {'GO:0004871': 'signal transducer activity'}
        fail_pfam = 15
        sample_ontology = pfam2go.GeneOntology('GO:0004871', 'signal transducer activity')
        sample_pfam = 'PF00015'
        with self.assertRaises(AssertionError):  # could likely do that more prettily?
            pfam2go.GeneOntologies(sample_pfam, fail_ontology)
        with self.assertRaises(AssertionError):
            pfam2go.GeneOntologies(fail_pfam, [sample_ontology])

    def test_gene_ontology(self):
        sample_go_id = 'GO:0004871'
        sample_description = 'signal transducer activity'
        sample_ontology = pfam2go.GeneOntology(sample_go_id, sample_description)
        assert sample_ontology.id == sample_go_id
        assert sample_ontology.description == sample_description
        assert str(sample_ontology) == sample_go_id

    def test_gene_ontology_fail(self):
        fail_id = '0004871'
        working_id = 'GO:0004871'
        working_description = 'signal transducer activity'
        fail_description = ['signal transducer activity']
        with self.assertRaisesRegex(ValueError, "Invalid Gene Ontology ID: 0004871"):
            pfam2go.GeneOntology(fail_id, working_description)
        with self.assertRaises(AssertionError):
            pfam2go.GeneOntology(working_id, fail_description)

    def test_build_as_i_go(self):
        data = path.get_full_path(os.path.dirname(__file__), 'data/pfam2go-march-2018.txt')
        ontologies_per_pfam = pfam2go.build_as_i_go(data)
        for ontology in ontologies_per_pfam.values():
            self.assertIsInstance(ontology, pfam2go.GeneOntologies)
        for pfam, go_ids in self.known_connections.items():
            self.assertEqual(str(go_ids), str(ontologies_per_pfam[pfam]))

    def test_get_gos(self):
        pfams = {'PF00015':FeatureLocation(0, 3)}
        fake_record = self.set_dummy_with_pfams(pfams)
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        for pfam, all_ontologies in gos_for_fake_pfam.items():
            for ontologies in all_ontologies:
                go_ids = [str(go_entry) for go_entry in ontologies.go_entries]
                for go_id in go_ids:
                    assert go_id in self.known_connections[ontologies.pfam]

    def test_blank_records(self):
        blank_no_pfams = DummyRecord()
        blank_no_ids = Record(Seq("ATGTTATGAGGGTCATAACAT", generic_dna))
        fake_pfam_location = FeatureLocation(0, 12)
        fake_pfam = PFAMDomain(location=fake_pfam_location, description='MCPsignal')
        blank_no_ids.add_pfam_domain(fake_pfam)
        with self.assertLogs(level='INFO') as log_cm:
            gos_for_no_pfams = pfam2go.get_gos_for_pfams(blank_no_pfams)
            assert 'No Pfam domains found' in str(log_cm.output)
            gos_for_no_ids = pfam2go.get_gos_for_pfams(blank_no_ids)
            assert 'No Pfam ids found' in str(log_cm.output)

    def test_get_gos_id_handling(self):
        pfams = {'PF00015.42': FeatureLocation(0, 3), 'PF00015_42': FeatureLocation(0, 3)}
        fake_record = self.set_dummy_with_pfams(pfams)
        # are wrong PFAM ids logged?
        with self.assertLogs() as log_cm:
            gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
            assert 'Pfam id PF00015_42 is not a valid Pfam id, skipping' in str(log_cm.output)
        for pfam_domain, all_ontologies in gos_for_fake_pfam.items():
            for ontologies in all_ontologies:
                # catches both stripping version number and broken ID not being included
                assert ontologies.pfam.isalnum()
                assert ontologies.pfam.startswith('PF')
                for ontology in ontologies.go_entries:
                    assert ontology.id in self.known_connections[ontologies.pfam]

    def test_results(self):
        pfams = {'PF00015': FeatureLocation(0, 3)}
        fake_record = self.set_dummy_with_pfams(pfams)
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        fake_results = pfam2go.Pfam2GoResults(fake_record.id, gos_for_fake_pfam)
        assert gos_for_fake_pfam == fake_results.pfam_domains_with_gos
        assert fake_record.id == fake_results.record_id
        for pfam, all_ontologies in fake_results.pfam_domains_with_gos.items():
            pfam_ids_without_versions = [pfam_id.partition('.')[0] for pfam_id in pfam.db_xref]
            for ontologies in all_ontologies:
                assert ontologies.pfam in pfam_ids_without_versions

    def test_add_results_to_record(self):
        pfams = {'PF00015.2': FeatureLocation(0, 3), 'PF00351.1': FeatureLocation(0, 3),
                 'PF00015.27': FeatureLocation(3, 6)}
        fake_record = self.set_dummy_with_pfams(pfams)
        fake_duplicate_pfam = PFAMDomain(location=FeatureLocation(6, 9), description='DUPLICATE')
        fake_duplicate_pfam.db_xref = ['PF00015.2']
        fake_record.add_pfam_domain(fake_duplicate_pfam)
        assert fake_duplicate_pfam in fake_record.get_pfam_domains()
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        fake_results = pfam2go.Pfam2GoResults(fake_record.id, gos_for_fake_pfam)
        fake_results.add_to_record(fake_record)
        assert fake_duplicate_pfam.db_xref == ['PF00015.2']
        for pfam in fake_record.get_pfam_domains():
            assert pfam.gene_ontologies['pfam2go'] == fake_results.pfam_domains_with_gos[pfam]
            # make sure identical pfams (with different version numbers) all have the same gene ontologies
            for pfam_id in pfam.db_xref:
                if pfam_id.startswith('PF00015'):
                    assert pfam.gene_ontologies['pfam2go'] == fake_results.pfam_domains_with_gos[fake_duplicate_pfam]

    def test_adding_to_wrong_record(self):
        pfams = {'PF00015': FeatureLocation(0, 3)}
        fake_record = self.set_dummy_with_pfams(pfams)
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        fake_results = pfam2go.Pfam2GoResults('FAKEID', gos_for_fake_pfam)
        with self.assertRaisesRegex(ValueError, "Record to store in and record analysed don't match"):
            fake_results.add_to_record(fake_record)

    def test_to_json(self):
        fake_pfam_location = FeatureLocation(0, 12)
        pfams = {'PF00015': fake_pfam_location, 'PF00351': fake_pfam_location}
        fake_record = self.set_dummy_with_pfams(pfams)
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        fake_results = pfam2go.Pfam2GoResults(fake_record.id, gos_for_fake_pfam)
        result_json = fake_results.to_json()
        expected_result = {"pfams": {"PF00015": [("GO:0004871", "signal transducer activity"),
                                                 ("GO:0007165", "signal transduction"),
                                                 ("GO:0016020", "membrane")],
                                     "PF00351": [("GO:0016714", "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced pteridine as one donor, and incorporation of one atom of oxygen"),
                                                 ("GO:0055114", "oxidation-reduction process")]},
                           "record_id": fake_record.id}
        assert result_json["record_id"] == expected_result["record_id"]
        for pfam in expected_result["pfams"]:
            assert expected_result["pfams"][pfam] == result_json["pfams"][pfam]

    def test_from_json(self):
        fake_pfam_location = FeatureLocation(0, 12)
        pfams = {'PF00015': fake_pfam_location, 'PF00351': fake_pfam_location, 'PF05147': fake_pfam_location}
        fake_record = self.set_dummy_with_pfams(pfams)
        gos_for_fake_pfam = pfam2go.get_gos_for_pfams(fake_record)
        fake_results = pfam2go.Pfam2GoResults(fake_record.id, gos_for_fake_pfam)
        result_json = fake_results.to_json()
        results_from_json = pfam2go.Pfam2GoResults.from_json(result_json, fake_record)
        assert 'PF05147' not in result_json["pfams"]
        for pfam in results_from_json.pfam_domains_with_gos:
            for pfam_id in pfam.db_xref:
                assert pfam_id in result_json["pfams"]
        from_json_to_json = results_from_json.to_json()
        assert result_json == from_json_to_json  # JSONception


if __name__ == '__main__':
    unittest.main()
