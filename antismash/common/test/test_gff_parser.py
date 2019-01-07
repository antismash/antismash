# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from unittest import TestCase
from Bio.SeqFeature import CompoundLocation

from antismash.common import errors, gff_parser, path
from antismash.common.test.helpers import get_simple_options, DummyRecord


class GffParserTest(TestCase):
    def setUp(self):
        self.config = get_simple_options(None, [])
        self.config.genefinding_gff3 = path.get_full_path(__file__, "data", "test_gff.gff")
        self.single_entry = False
        contig1 = DummyRecord(seq="A"*2000)
        contig1.id = "CONTIG_1"
        contig2 = DummyRecord(seq="A"*2000)
        contig2.id = "CONTIG_2"
        self.sequences = [contig1, contig2]

    def test_run(self):
        for sequence in self.sequences:
            gff_parser.run(sequence, self.single_entry, self.config)
        len_cds_1 = len(self.sequences[0].get_cds_features())
        len_cds_2 = len(self.sequences[1].get_cds_features())
        detected_result = (len_cds_1, len_cds_2)
        expected_result = (1, 0)
        self.assertEqual(detected_result, expected_result,
                         msg="\nResult : %s\nExpected : %s" % (detected_result, expected_result))

    def test_top_level_cds(self):
        self.config.genefinding_gff3 = path.get_full_path(__file__, "data", "single_cds.gff")
        gff_parser.run(self.sequences[0], self.single_entry, self.config)
        assert len(self.sequences[0].get_cds_features()) == 1

    def test_features_from_file(self):
        filename = path.get_full_path(__file__, 'data', 'fumigatus.cluster1.gff')
        record = DummyRecord()
        features = gff_parser.get_features_from_file(record, open(filename))
        assert len(features) == 11
        for feature in features:
            assert feature.type == 'CDS'
            assert isinstance(feature.location, CompoundLocation)

    def test_suitability(self):
        self.sequences[0].id = "NOT_CONTIG_1"
        with self.assertRaisesRegex(errors.AntismashInputError,
                                    "GFF3 record IDs don't match sequence file record IDs"):
            gff_parser.check_gff_suitability(self.config, self.sequences)

        # doesn't test very much
        self.sequences[0].id = "CONTIG_1"
        gff_parser.run(self.sequences[0], self.single_entry, self.config)  # insert the features
        assert not gff_parser.check_gff_suitability(self.config, self.sequences)

        # test force correlation
        self.sequences = self.sequences[1:]  # CONTIG_2
        assert gff_parser.check_gff_suitability(self.config, self.sequences)
