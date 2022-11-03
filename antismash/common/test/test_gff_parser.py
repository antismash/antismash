# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

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
        assert {feat.type for feat in first} == {
            "CDS",
            "exon",
            "mRNA",
            "five_prime_UTR",
            "gene",
        }
        assert len(first) == 7  # 10 features in the file, but the CDS features combine
        assert isinstance(first[0], SeqFeature)
        assert first[0].type == "gene"
        assert isinstance(first[6], SeqFeature)
        assert first[6].type == "CDS"
        # ensure the CDS components of the gene are properly combined
        assert str(first[6].location) == "join{[1124:1291](+), [1396:1450](+), [1753:1871](+)}"
        assert len(list(filter(lambda x: x.type == "CDS", first))) == 1
        # and check parent linkages are correct
        for feature in first[1:]:  # don't look at the gene if we're looking at gene references
            assert feature.qualifiers["gene"] == first[0].qualifiers["Name"]

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
