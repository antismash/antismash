# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"Test suite for gff_parser module"

import os
from unittest import TestCase
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation

import antismash
from antismash.common import deprecated, gff_parser, path
from antismash.common.test.helpers import FakeRecord, FakeFeature

class GffParserTest(TestCase):
    def setUp(self):
        self.config = antismash.config.args.simple_options(None, [])
        self.config.gff3 = path.get_full_path(__file__, "data/test_gff.gff")
        self.config.single_entries = False
        contig1 = FakeRecord(seq="A"*2000, real_seq=True)
        contig1.id = "CONTIG_1"
        contig2 = FakeRecord(seq="A"*2000, real_seq=True)
        contig2.id = "CONTIG_2"
        self.sequences = [contig1, contig2]

    def test_run(self):
        for sequence in self.sequences:
            gff_parser.run(sequence, self.config)
        len_cds_1 = len(deprecated.get_cds_features(self.sequences[0]))
        len_cds_2 = len(deprecated.get_cds_features(self.sequences[1]))
        detected_result = (len_cds_1, len_cds_2)
        expected_result = (1, 0)
        self.assertEqual(detected_result, expected_result,
                         msg="\nResult : %s\nExpected : %s" % (detected_result, expected_result))

    def test_top_level_cds(self):
        self.config.gff3 = path.get_full_path(__file__, "data/single_cds.gff")
        gff_parser.run(self.sequences[0], self.config)
        assert len(deprecated.get_cds_features(self.sequences[0])) == 1

    def test_features_from_file(self):
        filename = path.get_full_path(__file__, 'data/fumigatus.cluster1.gff')
        record = FakeRecord(real_seq=True)
        features = gff_parser.get_features_from_file(record, open(filename))
        assert len(features) == 11
        for feature in features:
            assert feature.type == 'CDS'
            assert isinstance(feature.location, CompoundLocation)

    def test_suitability(self):
        self.config.all_record_ids = []
        with self.assertRaises(ValueError) as err:
            gff_parser.check_gff_suitability(self.config, self.sequences)
        assert "GFF3 record IDs don't match sequence file record IDs" in str(err.exception)

        # doesn't test very much
        gff_parser.run(self.sequences[0], self.config) # insert the features
        self.config.all_record_ids = ['CRO_000001']
        gff_parser.check_gff_suitability(self.config, self.sequences)
