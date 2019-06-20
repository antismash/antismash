# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from unittest import TestCase
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from antismash.common import errors, gff_parser, path


class GffParserTest(TestCase):
    def setUp(self):
        self.gff_file = path.get_full_path(__file__, "data", "test_gff.gff")
        self.single_entry = False
        contig1 = SeqRecord(seq="A"*2000)
        contig1.id = "CONTIG_1"
        contig2 = SeqRecord(seq="A"*2000)
        contig2.id = "CONTIG_2"
        self.sequences = [contig1, contig2]

    def test_run(self):
        results = gff_parser.run(self.gff_file)
        first = results["CONTIG_1"]
        assert len(first) == 1
        assert isinstance(first[0], SeqFeature)

        assert "CONTIG_2" not in results

    def test_top_level_cds(self):
        self.gff_file = path.get_full_path(__file__, "data", "single_cds.gff")
        cds_features = gff_parser.run(self.gff_file)["CONTIG_1"]
        assert len(cds_features) == 1

    def test_features_from_file(self):
        filename = path.get_full_path(__file__, 'data', 'fumigatus.cluster1.gff')
        features = gff_parser.get_features_from_file(open(filename))["cluster0"]
        assert len(features) == 11
        for feature in features:
            assert feature.type == 'CDS'
            assert isinstance(feature.location, CompoundLocation)

    def test_suitability(self):
        self.sequences[0].id = "NOT_CONTIG_1"
        with self.assertRaisesRegex(errors.AntismashInputError,
                                    "GFF3 record IDs don't match sequence file record IDs"):
            gff_parser.check_gff_suitability(self.gff_file, self.sequences)

        self.sequences[0].id = "CONTIG_1"
        gff_parser.check_gff_suitability(self.gff_file, self.sequences)

        # test force correlation
        self.sequences = self.sequences[1:]  # CONTIG_2
        gff_parser.check_gff_suitability(self.gff_file, self.sequences)
