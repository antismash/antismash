# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from io import StringIO
from unittest import TestCase
from unittest.mock import patch

from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from antismash.common import errors, gff_parser, path
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation


class GffParserTest(TestCase):
    def setUp(self):
        self.gff_file = path.get_full_path(__file__, "data", "test_gff.gff")
        self.single_entry = False
        contig1 = SeqRecord(seq=Seq("A"*2000))
        contig1.id = "CONTIG_1"
        contig2 = SeqRecord(seq=Seq("A"*2000))
        contig2.id = "CONTIG_2"
        self.sequences = [contig1, contig2]

    def test_run(self):
        gff_parser.update_records(self.gff_file, self.sequences)
        first = self.sequences[0].features
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

        assert not self.sequences[1].features

    def test_top_level_cds(self):
        self.gff_file = path.get_full_path(__file__, "data", "single_cds.gff")
        gff_parser.update_records(self.gff_file, self.sequences)
        assert len(self.sequences[0].features) == 1
        assert self.sequences[0].features[0].type == "CDS"

    def test_features_from_file(self):
        filename = path.get_full_path(__file__, 'data', 'fumigatus.cluster1.gff')
        features = gff_parser.get_features_from_file(open(filename, encoding="utf-8"))["cluster0"]
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

    def test_circular_suitability(self):
        # the suitability check should not reject coordinates outside the record
        # at this stage, that happens with feature parsing when circularity is known
        record = SeqRecord(id="DUMMY.1", seq="A" * 300)
        text = (
            f"{record.id}	Genbank	region	1	400	.	+	.	ID=DUMMY.1:1..400;Is_circular=true\n"
            f"{record.id}	Genbank	CDS	200	500	.	+	0	ID=cds-1"
        )

        def dummy_open(*_args, **_kwargs):
            return StringIO(text)

        with patch("builtins.open", side_effect=dummy_open) as mocked:
            gff_parser.check_gff_suitability("dummy_filename", [record])
            mocked.assert_called()

    def test_any_have_circular(self):
        explicit_circular = {"Is_circular": ["true"]}
        explicit_linear = {"Is_circular": ["false"]}
        unknown_state = {"key": ["val"]}
        dummy = SeqFeature(type="dummy")

        # without a relevant type
        # no relevant annotations
        feature = SeqFeature(FeatureLocation(1, 10, 1), type="CDS", qualifiers={})
        assert not gff_parser.any_have_circularity([feature, dummy])
        # relevant annotations
        feature.qualifiers = explicit_circular
        assert not gff_parser.any_have_circularity([dummy, feature])

        # relevant type
        # no relevant annotations
        feature.type = "source"  # at this stage, the type has changed from raw GFF
        feature.qualifiers = unknown_state
        assert not gff_parser.any_have_circularity([dummy, feature])
        # explicitly circular
        feature.qualifiers = explicit_circular
        assert gff_parser.any_have_circularity([feature, dummy])
        # explicity linear
        feature.qualifiers = explicit_linear
        assert not gff_parser.any_have_circularity([dummy, feature])

    def test_region_feature_rename(self):
        record_id = "DUMMY"
        text = StringIO(
            f"{record_id}	Genbank	region	1	400	.	+	.	ID=DUMMY.1:1..400;\n"
            f"{record_id}	Genbank	CDS	200	300	.	+	0	ID=cds-1"
        )
        features_by_record = gff_parser.get_features_from_file(text)
        assert len(features_by_record) == 1
        features = features_by_record[record_id]
        assert features[0].type == "source"
        assert features[0].location == FeatureLocation(0, 400, 1)
        assert features[1].type == "CDS"
        assert features[1].location == FeatureLocation(199, 300, 1)

    def test_multiple_region_features(self):
        record = SeqRecord(id="DUMMY.1", seq="A" * 300)
        text = (
            f"{record.id}	Genbank	region	1	400	.	+	.	ID=DUMMY.1:1..400;Is_circular=true\n"
            f"{record.id}	Genbank	region	500	700	.	+	0	ID=DUMMY.1:500..700;Is_circular=true"
        )

        def dummy_open(*_args, **_kwargs):
            return StringIO(text)

        with patch("builtins.open", side_effect=dummy_open) as mocked:
            with self.assertRaisesRegex(errors.AntismashInputError, "feature already defined"):
                gff_parser.update_records("dummy_filename", [record])
            mocked.assert_called()


class TestSplitLocation(TestCase):
    def setUp(self):
        self.feature = SeqFeature()

    def split(self, length, features=None):
        if not features:
            features = [self.feature]
        assert all(feature.location for feature in features)
        return gff_parser.split_cross_origin_locations(features, length)

    def test_all_after_point(self):
        self.feature.location = FeatureLocation(20, 30, 1)
        with self.assertRaisesRegex(ValueError, "entirely outside record"):
            self.split(10)
        self.feature.location = CompoundLocation([FeatureLocation(20, 30, 1), FeatureLocation(40, 50, 1)])
        with self.assertRaisesRegex(ValueError, "entirely outside record"):
            self.split(10)

    def test_simple_all_before_point(self):
        original = FeatureLocation(1, 5, 1)
        self.feature.location = original
        self.split(10)
        assert self.feature.location == original

    def test_compound_all_before_point(self):
        original = CompoundLocation([FeatureLocation(1, 5, 1), FeatureLocation(8, 12, 1)])
        self.feature.location = original
        self.split(20)
        assert self.feature.location == original

    def test_point_between_exons(self):
        original = CompoundLocation([FeatureLocation(11, 15, 1), FeatureLocation(18, 22, 1)])
        expected = CompoundLocation([FeatureLocation(11, 15, 1), FeatureLocation(1, 5, 1)])
        self.feature.location = original
        self.split(17)
        assert self.feature.location == expected

    def test_point_in_exons(self):
        original = CompoundLocation([
            FeatureLocation(11, 15, 1),
            FeatureLocation(18, 22, 1),
            FeatureLocation(25, 30, 1),
        ])
        expected = CompoundLocation([
            FeatureLocation(11, 15, 1),
            FeatureLocation(18, 20, 1),
            FeatureLocation(0, 2, 1),
            FeatureLocation(5, 10, 1),
        ])
        self.feature.location = original
        self.split(20)
        assert self.feature.location == expected

    def test_point_in_exons_reverse(self):
        original = CompoundLocation([
            FeatureLocation(120, 130, -1),
            FeatureLocation(80, 110, -1),
            FeatureLocation(60, 70, -1),
        ])
        expected = CompoundLocation([
            FeatureLocation(20, 30, -1),
            FeatureLocation(0, 10, -1),
            FeatureLocation(80, 100, -1),
            FeatureLocation(60, 70, -1),
        ])
        self.feature.location = original
        self.split(100)
        assert self.feature.location == expected


class TestUpdateRecords(TestCase):
    def setUp(self):
        self.record = SeqRecord(seq="A" * 100, id="dummy")
        self.circular_feature = SeqFeature(
            location=FeatureLocation(0, len(self.record), 1),
            type="source",
            qualifiers={"Is_circular": ["true"]},
        )

    def run_faked_file(self, features, records=None):
        if records is None:
            records = [self.record]
        with patch.object(gff_parser, "get_features_from_file", return_value=features) as patched:
            with patch("builtins.open"):
                gff_parser.update_records("unused_filename", records)
            patched.assert_called_once()

    def test_no_records(self):
        with self.assertRaisesRegex(ValueError, "no records provided"):
            self.run_faked_file({self.record.id: [self.circular_feature]}, records=[])

    def test_update_with_existing_topology(self):
        # existing linear topology conflicts with the GFF feature's circularity
        self.record.annotations["topology"] = "linear"
        with self.assertRaisesRegex(errors.AntismashInputError, "existing and incompatible topology"):
            self.run_faked_file({self.record.id: [self.circular_feature]})
        # existing circular topology is fine
        self.record.annotations["topology"] = "circular"
        self.run_faked_file({self.record.id: [self.circular_feature]})

    def test_update_sets_circular(self):
        assert self.record.annotations.get("topology") != "circular"
        self.run_faked_file({self.record.id: [self.circular_feature]})
        assert self.record.annotations["topology"] == "circular"
