# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from antismash.common import path
from antismash.modules import pfam2go


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

    def test_parse_mappings(self):
        data = path.get_full_path(os.path.dirname(__file__), 'pfam2go-march-2018.txt')
        print(data)
        go_ids_by_pfam, go_desc_by_id = pfam2go.parse_all_mappings(data)
        # make sure all go ids actually have associated descriptions
        for go_ids_found in go_ids_by_pfam.values():
            for go_id in go_ids_found:
                assert go_id in go_desc_by_id
        # see if GO ids work
        for pfam, go_ids in self.known_connections.items():
            self.assertEqual(go_ids, go_ids_by_pfam[pfam])

    def test_build_at_the_end(self):
        data = path.get_full_path(os.path.dirname(__file__), 'pfam2go-march-2018.txt')
        go_ids_by_pfam, go_desc_by_id = pfam2go.parse_all_mappings(data)
        # for now: sanity check, is this even working off good data?
        for go_ids_found in go_ids_by_pfam.values():
            for go_id in go_ids_found:
                self.assertIn(go_id, go_desc_by_id)
        ontologies_per_pfam = pfam2go.build_at_the_end(go_ids_by_pfam, go_desc_by_id)
        for ontologies in ontologies_per_pfam.values():
            self.assertIsInstance(ontologies, pfam2go.GeneOntologies)

    def test_construct_ontologies(self):
        data = path.get_full_path(os.path.dirname(__file__), 'pfam2go-march-2018.txt')
        go_ids_by_pfam, go_desc_by_id = pfam2go.parse_all_mappings(data)
        for pfam in go_ids_by_pfam:
            ontology_tested = pfam2go.construct_gene_ontologies(pfam, go_ids_by_pfam[pfam], go_desc_by_id)
            self.assertIsInstance(ontology_tested, pfam2go.GeneOntologies)

    def test_gene_ontologies(self):
        # does it use arguments given? How is bad input handled?
        sample_ontology = pfam2go.GeneOntology('GO:0004871', 'signal transducer activity')
        sample_pfam = 'PF00015'
        sample_ontologies = pfam2go.GeneOntologies(sample_pfam, [sample_ontology])
        assert sample_ontologies.pfam == sample_pfam
        assert sample_ontologies.go_entries == [sample_ontology]
        # check the string stuff too?

    def test_gene_ontologies_fail(self):
        fail_ontology = {'GO:0004871': 'signal transducer activity'}
        fail_pfam = 15
        with self.assertRaises(AssertionError):  # could likely do that more prettily?
            pfam2go.GeneOntologies(fail_pfam, fail_ontology)

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
        with self.assertRaisesRegex(ValueError, "Invalid Gene Ontology ID"):
            pfam2go.GeneOntology(fail_id, working_description)
        with self.assertRaises(AssertionError):
            pfam2go.GeneOntology(working_id, fail_description)

    def test_ontology_fails_when_bad_input(self):
        bad_ids = {'PF00015': ['GO:004871', 'GO:0007165', 'GO:0016020']}
        data = path.get_full_path(os.path.dirname(__file__), 'pfam2go-march-2018.txt')
        go_desc_by_id = pfam2go.parse_all_mappings(data)[1]
        with self.assertRaisesRegex(ValueError, "Gene Ontology ID has no associated description: GO:004871"):
            pfam2go.construct_gene_ontologies('PF00015', bad_ids['PF00015'], go_desc_by_id)

    def test_build_as_i_go(self):
        data = path.get_full_path(os.path.dirname(__file__), 'pfam2go-march-2018.txt')
        ontologies_per_pfam = pfam2go.build_as_i_go(data)
        for ontology in ontologies_per_pfam.values():
            self.assertIsInstance(ontology, pfam2go.GeneOntologies)
        for pfam, go_ids in self.known_connections.items():
            self.assertEqual(str(go_ids), str(ontologies_per_pfam[pfam]))



if __name__ == '__main__':
    unittest.main()
