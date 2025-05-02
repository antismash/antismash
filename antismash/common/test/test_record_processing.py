# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from tempfile import NamedTemporaryFile
import os
import unittest
from unittest import mock

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from antismash import config
from antismash.common import record_processing, path
from antismash.common.errors import AntismashInputError
from antismash.common.secmet import Record
from antismash.common.secmet.test.helpers import DummyCDSMotif, DummyFeature, DummyRecord
from antismash.common.test import helpers


class TestParseRecords(unittest.TestCase):
    def test_nisin(self):
        nisin_path = helpers.get_path_to_nisin_genbank()
        records = record_processing.parse_input_sequence(nisin_path)
        assert len(records) == 1
        assert isinstance(records[0], Record)
        assert len(records[0].get_cds_features()) == 11
        assert len(records[0].seq) == 15016

    def test_trim(self):
        nisin_path = helpers.get_path_to_nisin_genbank()
        records = record_processing.parse_input_sequence(nisin_path, start=10, end=5000)
        assert len(records) == 1
        assert isinstance(records[0], Record)
        assert len(records[0].get_cds_features()) == 2
        assert len(records[0].seq) == 4990

    def test_cross_origin_trim(self):
        record = DummyRecord(length=20, circular=True)
        bio_record = record.to_biopython()
        with self.assertRaisesRegex(AntismashInputError, "cannot be used for a cross-origin"):
            with mock.patch.object(record_processing, "_strict_parse", return_value=[bio_record]):
                record_processing.parse_input_sequence("dummy_filename", start=15, end=5)

    def test_minimum_length(self):
        nisin_path = helpers.get_path_to_nisin_genbank()
        records = record_processing.parse_input_sequence(nisin_path,
                                                         minimum_length=-16)
        assert len(records) == 1

        records = record_processing.parse_input_sequence(nisin_path,
                                                         minimum_length=15016)
        assert len(records) == 1

        with self.assertRaisesRegex(AntismashInputError, "smaller than minimum length"):
            record_processing.parse_input_sequence(nisin_path, minimum_length=15017)

        for bad_len in [5.6, None, "5"]:
            with self.assertRaisesRegex(TypeError, "minimum_length must be an int"):
                record_processing.parse_input_sequence(path, minimum_length=bad_len)

    def test_nonexistant(self):
        with self.assertRaisesRegex(AntismashInputError, "No such file or directory"):
            record_processing.parse_input_sequence("does_not_exist.gbk")

    def test_empty(self):
        with NamedTemporaryFile(suffix=".gbk") as temp:
            with self.assertRaisesRegex(AntismashInputError, "no valid records found"):
                record_processing.parse_input_sequence(temp.name)

    def test_malformed_qualifiers(self):
        filepath = path.get_full_path(__file__, "data", "incomplete_qualifier.gbk")
        warning = "double-quote characters like \" should be escaped"
        with self.assertRaisesRegex(AntismashInputError, warning):
            record_processing.parse_input_sequence(filepath)

    def test_partial_invalid(self):
        # pretend there's two inputs, and that the first has an error
        # with the ignore invalid option turned on, the valid input should still come back
        dummy_error = record_processing.SecmetInvalidInputError("some reason")
        dummy_inputs = [SeqRecord(Seq("ACGT"), id="bad"), SeqRecord(Seq("TGCA"), id="good")]
        expected_outputs = [dummy_error, Record.from_biopython(dummy_inputs[1], taxon="bacteria")]
        for value in [False, True]:
            with mock.patch.object(record_processing, "_strict_parse", return_value=dummy_inputs):
                with mock.patch.object(Record, "from_biopython", side_effect=expected_outputs):
                    # with the option off, it should be an error
                    if not value:
                        with self.assertRaises(AntismashInputError):
                            record_processing.parse_input_sequence("dummy file", ignore_invalid_records=value)
                    else:
                        records = record_processing.parse_input_sequence("dummy file", ignore_invalid_records=value)
                        assert len(records) == 1
                        assert records[0].id == dummy_inputs[1].id

    def test_all_invalid_and_ignored(self):
        dummy_error = record_processing.SecmetInvalidInputError("some reason")
        dummy_inputs = [SeqRecord(Seq("ACGT"), id="bad")]
        with mock.patch.object(Record, "from_biopython", side_effect=dummy_error):
            with mock.patch.object(record_processing, "_strict_parse", return_value=dummy_inputs):
                with self.assertRaisesRegex(AntismashInputError, "no valid records"):
                    record_processing.parse_input_sequence("dummy file", ignore_invalid_records=True)

    def test_absolute_name_limit(self):
        record = DummyRecord()
        record.id = "A" * os.pathconf("/", "PC_NAME_MAX")
        with mock.patch.object(record_processing, "_strict_parse", return_value=[record]):
            with self.assertRaisesRegex(AntismashInputError, "too long"):
                record_processing.parse_input_sequence("dummy file")


class TestGapNotation(unittest.TestCase):
    def setUp(self):
        config.build_config([], isolated=True)

    def tearDown(self):
        config.destroy_config()

    def test_no_conversion_required(self):
        gaps = "acgtNacgtN"
        no_gaps = "acgtacgt"
        records = [Record(Seq(gaps)), Record(Seq(no_gaps))]
        records = [record_processing.sanitise_sequence(record) for record in records]
        assert len(records) == 2
        assert records[0].seq == gaps.upper()
        assert records[1].seq == no_gaps.upper()

    def test_conversion(self):
        for seq in ("acgtZacgtX", "acgtxacgtz", "acgtoacgtf"):
            record = Record(Seq(seq))
            record = record_processing.sanitise_sequence(record)
            assert record.seq == "ACGTNACGTN"

    def test_removal(self):
        record = Record(Seq("acg-ta--"))
        record = record_processing.sanitise_sequence(record)
        assert record.seq == "ACGTA"

        # test modified in place
        record = Record(Seq("acg-ta--"))
        record_processing.sanitise_sequence(record)
        assert record.seq == "ACGTA"

    def test_mix(self):
        record = Record(Seq("acg-ta-F-x"))
        record = record_processing.sanitise_sequence(record)
        assert record.seq == "ACGTANN"


class TestTrimSequence(unittest.TestCase):
    def setUp(self):
        self.seq = "0123456789"
        self.record = SeqRecord(Seq("0123456789"))

    def trim_seq(self, start, end):
        new = record_processing.trim_sequence(self.record, start, end)
        # since func called, alter start/end for these slice checks
        start = max(start, 0)
        if end < 0:
            end = len(self.seq)
        # check that Bio.Seq slices as expected
        assert str(self.record.seq[start:end]) == self.seq[start:end]
        # check that the record's seq was sliced
        assert str(new.seq) == self.seq[start:end]
        # return for further checking
        return new

    def test_bad_values(self):
        # start outside seq
        with self.assertRaises(ValueError):
            record_processing.trim_sequence(self.record, start=10, end=-1)
        # end outside seq
        with self.assertRaises(ValueError):
            record_processing.trim_sequence(self.record, start=-1, end=20)
        # start after end
        with self.assertRaises(ValueError):
            record_processing.trim_sequence(self.record, start=7, end=4)

    def test_good_values(self):
        self.trim_seq(-1, 5)
        self.trim_seq(4, -1)
        self.trim_seq(2, 9)

    def test_features_removed(self):
        self.record.features.extend([SeqFeature(FeatureLocation(s, e)) for s, e in [
                                                    (1, 4), (4, 7), (8, 10)]])
        assert len(self.trim_seq(-1, -1).features) == len(self.record.features)
        assert len(self.trim_seq(-1, 7).features) == 2
        assert len(self.trim_seq(-1, 5).features) == 1

        assert len(self.trim_seq(0, 5).features) == 1
        assert len(self.trim_seq(3, 7).features) == 1
        assert len(self.trim_seq(3, 9).features) == 1

        assert not self.trim_seq(9, 10).features

    def test_no_change(self):
        new = record_processing.trim_sequence(self.record, start=-1, end=-1)
        assert new.seq == self.record.seq


class TestPreprocessRecords(unittest.TestCase):
    def setUp(self):
        class DummyModule:
            def __init__(self):
                self.was_run = False

            def get_arguments(self):
                args = config.args.ModuleArgs("genefinding", "genefinding")
                args.add_option("gff3", default="", type=str, help="dummy", dest="gff3")
                args.add_option("tool", default="", type=str, help="dummy", dest="tool")
                return args

            def run_on_record(self, _record, options):
                assert options.genefinding_gff3 is not None
                assert options.genefinding_tool is not None
                assert options.taxon
                self.was_run = True

        self.genefinding = DummyModule()

        options = config.build_config(["--cpus", "1"], isolated=True, modules=[self.genefinding])
        config.update_config({"triggered_limit": False})
        self.options = options

    def read_nisin(self):
        records = record_processing.parse_input_sequence(helpers.get_path_to_nisin_genbank())
        assert len(records) == 1
        return records

    def read_double_nisin(self):
        records = self.read_nisin() + self.read_nisin()
        assert len(records) == 2
        return records

    def run_on_records(self, records):
        record_processing.pre_process_sequences(records, self.options, self.genefinding)

    def tearDown(self):
        config.destroy_config()

    def test_nisin_genbank(self):
        self.run_on_records(self.read_nisin())
        assert not self.genefinding.was_run

    def test_nisin_fasta_only(self):
        config.update_config({"genefinding_tool": "none"})
        filepath = path.get_full_path(__file__, "data", "nisin.fasta")
        records = record_processing.parse_input_sequence(filepath)
        assert len(records) == 1
        assert not records[0].get_cds_features()
        # make sure genefinding wasn't run with default options
        with self.assertRaisesRegex(AntismashInputError, "all records skipped"):
            record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert not self.genefinding.was_run
        assert not records[0].get_cds_features()

        # make sure genefinding was run when not 'none'
        records[0].skip = False
        config.update_config({"genefinding_tool": "not-none"})
        # due to no genes actually being marked, it'll raise an error
        with self.assertRaisesRegex(AntismashInputError, "all records skipped"):
            record_processing.pre_process_sequences(records, self.options, self.genefinding)
        # but genefinding was still run
        assert self.genefinding.was_run
        # still no features because we used dummy genefinding
        for record in records:
            assert not record.get_cds_features()
            assert record.skip.lower() == "no genes found"

    def test_nisin_fasta_gff(self):
        fasta = path.get_full_path(__file__, "data", "nisin.fasta")
        gff = path.get_full_path(__file__, "data", "nisin.gff3")
        config.update_config({"genefinding_gff3": gff})
        records = record_processing.parse_input_sequence(fasta, gff_file=gff)
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert not self.genefinding.was_run
        assert len(records[0].get_cds_features()) == 11

    def test_shotgun(self):
        filepath = path.get_full_path(__file__, "data", "wgs.gbk")
        with self.assertRaisesRegex(AntismashInputError, "incomplete whole genome shotgun records are not supported"):
            record_processing.parse_input_sequence(filepath)

    def test_duplicate_record_ids(self):
        records = self.read_double_nisin()
        assert records[0].id == records[1].id
        records = record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert len(records) == 2
        assert records[0].id != records[1].id

    def test_limit(self):
        records = self.read_double_nisin()
        assert all(rec.skip is None for rec in records)
        assert not self.options.triggered_limit
        config.update_config({"limit": 1})
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert records[0].skip is None
        assert records[1].skip.startswith("skipping all but largest 1")
        assert self.options.triggered_limit

    def test_limit_uses_largest(self):
        cds = helpers.DummyCDS(start=1, end=7)
        records = [helpers.DummyRecord(seq="A" * i, features=[cds]) for i in [50, 10, 70]]
        config.update_config({
            "limit": 2,
            "minlength": 0,
        })
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert records[0].skip is None
        assert records[1].skip.startswith("skipping all but largest")
        assert records[2].skip is None
        assert self.options.triggered_limit

    def test_limit_to_record_partial(self):
        records = self.read_double_nisin()
        assert all(rec.skip is None for rec in records)
        config.update_config({"limit_to_record": records[0].id})
        records[0].id += "_changed"
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert not records[1].skip
        assert records[0].skip.startswith("did not match filter")

    def test_limit_to_record_complete(self):
        records = self.read_double_nisin()
        config.update_config({"limit_to_record": "bad_id"})
        with self.assertRaisesRegex(AntismashInputError, "no sequences matched filter"):
            record_processing.pre_process_sequences(records, self.options, self.genefinding)

    def test_no_record_id(self):
        records = self.read_nisin()
        records[0].id = ""
        with self.assertRaisesRegex(AntismashInputError, "record has no name"):
            self.run_on_records(records)

    def test_long_names(self):
        record = self.read_nisin()[0]
        record.id = "A" * 16
        record.name = record.id
        self.run_on_records([record])
        assert record.id == record.name == "A" * 16

        record.id = "A" * 17
        record.name = record.id
        self.run_on_records([record])
        assert record.id == record.name == "A" * 17

        config.update_config({"allow_long_headers": False})
        record.id = "A" * 17
        record.name = record.id
        self.run_on_records([record])
        assert len(record.id) <= 16
        assert len(record.name) <= 16

    def test_mac_bad_parallel(self):
        """ For python 3.8+ on mac, parallel behaviour of Pool() changed
            from using fork() to spawn(), which meant that config objects
            failed to serialise in the same way.
            This test exists to ensure that such failures are worked around
            for at least the record_processing portions.
        """
        filepath = path.get_full_path(__file__, "data", "nisin.fasta")
        records = record_processing.parse_input_sequence(filepath)

        # to mimic the parallel section not keeping config, wrap the function
        # called in the subprocesses and destroy the config prior to the original
        # being called
        original = record_processing.ensure_cds_info

        def wrapper(*args, **kwargs):
            config.destroy_config()
            assert len(self.options) == 0
            results = original(*args, **kwargs)
            records[0].skip = False  # prevent raising an input error so the path completes fully
            return results

        # then run the full process with the newly wrapped function
        with mock.patch.object(record_processing, "ensure_cds_info", wraps=wrapper):
            self.run_on_records(records)
        # then ensure genefinding would have run normally
        assert self.genefinding.was_run

    def test_missing_sequence(self):
        # a record missing a sequence shouldn't crash in WGS testing
        # and shouldn't report as a WGS
        with NamedTemporaryFile(suffix=".fasta") as handle:
            handle.write(">R1\nACGT\n>R2\n\n>R3\nACGT\n".encode())
            handle.flush()
            with self.assertRaisesRegex(AntismashInputError, "no sequence .*R2.*"):
                record_processing.parse_input_sequence(handle.name)


class TestUniqueID(unittest.TestCase):
    def test_bad_starts(self):
        for bad_start in ["start", None, {}, []]:
            with self.assertRaises((ValueError, TypeError)):
                record_processing.generate_unique_id("pref", [], bad_start)

    def test_bad_collections(self):
        for bad_existing in [None, object(), 2]:
            with self.assertRaises((ValueError, TypeError)):
                record_processing.generate_unique_id("pref", bad_existing, 1)

    def test_bad_max(self):
        for bad_max in ["start", None, {}, []]:
            with self.assertRaises((ValueError, TypeError)):
                record_processing.generate_unique_id("pref", {}, 1, bad_max)

    def test_generation(self):
        existing = {f"a_{i}" for i in range(15)}
        new, counter = record_processing.generate_unique_id("a", existing)
        assert len(existing) == 15 and new not in existing
        assert new == "a_15" and counter == 15

        new, counter = record_processing.generate_unique_id("a", existing, start=17)
        assert len(existing) == 15 and new not in existing
        assert new == "a_17" and counter == 17

        new, counter = record_processing.generate_unique_id("b", existing)
        assert len(existing) == 15 and new not in existing
        assert new == "b_0" and counter == 0

    def test_overlong(self):
        existing = {f"a_{i}" for i in range(150)}
        # prefix itself too long
        with self.assertRaisesRegex(RuntimeError, "Could not generate .*"):
            record_processing.generate_unique_id("aaa", existing, start=0, max_length=3)

        # the generated number is too long
        with self.assertRaisesRegex(RuntimeError, "Could not generate .*"):
            record_processing.generate_unique_id("a", existing, start=140, max_length=4)

    def test_existing_format_markers(self):
        for marker in ["%s", "{}", "{name}"]:
            name = f"{marker}_other"
            result, _ = record_processing.generate_unique_id(name, {name}, start=0)
            assert result == name + "_0"


class TestStripRecord(unittest.TestCase):
    def test_cds_motifs(self):
        record = helpers.DummyRecord()

        non_as_motif = DummyFeature(feature_type="CDS_motif")
        non_as_motif._qualifiers["stuff"] = ["thing"]
        record.add_feature(non_as_motif)

        as_motif = DummyCDSMotif(domain_id="as")
        as_motif.created_by_antismash = True
        record.add_cds_motif(as_motif)

        bio = record.to_biopython()
        assert len(bio.features) == 2
        assert set(feat.type for feat in bio.features) == {"CDS_motif"}
        assert record_processing.strip_record(bio) is bio  # strips and returns
        assert len(bio.features) == 1
        record = Record.from_biopython(bio, taxon="bacteria")

        motifs = record.get_cds_motifs()
        assert len(motifs) == 1
        assert motifs[0].domain_id.startswith("non_aS_motif")
