# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from tempfile import NamedTemporaryFile
import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from antismash import config
from antismash.common import record_processing, path
from antismash.common.secmet import Record
from antismash.common.test import helpers

class TestParseRecords(unittest.TestCase):
    def test_nisin(self):
        path = helpers.get_path_to_nisin_genbank()
        records = record_processing.parse_input_sequence(path)
        assert len(records) == 1
        assert isinstance(records[0], Record)
        assert len(records[0].get_cds_features()) == 11
        assert len(records[0].seq) == 15016

    def test_minimum_length(self):
        path = helpers.get_path_to_nisin_genbank()
        records = record_processing.parse_input_sequence(path, minimum_length=-16)
        assert len(records) == 1

        records = record_processing.parse_input_sequence(path, minimum_length=15016)
        assert len(records) == 1

        records = record_processing.parse_input_sequence(path, minimum_length=15017)
        assert not records

        for bad_len in [5.6, None, "5"]:
            with self.assertRaisesRegex(TypeError, "minimum_length must be an int"):
                record_processing.parse_input_sequence(path, minimum_length=bad_len)

    def test_nonexistant(self):
        with self.assertRaisesRegex(ValueError, "Sequence file not found: .*"):
            record_processing.parse_input_sequence("does_not_exist.gbk")

    def test_empty(self):
        with NamedTemporaryFile(suffix=".gbk") as temp:
            records = record_processing.parse_input_sequence(temp.name)
            assert not records


class TestGapNotation(unittest.TestCase):
    def test_no_conversion_required(self):
        gaps = "acgtNacgtN"
        no_gaps = "acgtacgt"
        records = [Record(Seq(gaps)), Record(Seq(no_gaps))]
        record_processing.ensure_gap_notation_consistent(records)
        assert len(records) == 2
        assert records[0].seq == gaps
        assert records[1].seq == no_gaps

    def test_conversion(self):
        for seq in ("acgtXacgtX", "acgtxacgtx", "acgt-acgt-"):
            record = Record(Seq(seq))
            record_processing.ensure_gap_notation_consistent([record])
            assert record.seq == "acgtNacgtN"

class TestTrimSequence(unittest.TestCase):
    def setUp(self):
        self.seq = "0123456789"
        self.record = SeqRecord(Seq("0123456789"))

    def trim_seq(self, start, end):
        new = record_processing.trim_sequence(self.record, start, end)
        # since func called, alter start/end for these slice checks
        if start < 0:
            start = 0
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

class TestIsNuclSeq(unittest.TestCase):
    def test_seq(self):
        # > 20%
        for seq in ["AGTC", "AGCTFC", "agtcfc", "AGTCFCT"]:
            assert record_processing.is_nucl_seq(Seq(seq))
        # edge case == 20% should be failure
        assert not record_processing.is_nucl_seq(Seq("AGFTC"))
        # and less than 20%
        assert not record_processing.is_nucl_seq(Seq("AGFTCF"))

    def test_str(self):
        for seq in ["AGTC", "AGCTFC", "agtcfc", "AGTCFCT"]:
            assert record_processing.is_nucl_seq(seq)
        assert not record_processing.is_nucl_seq("AGFTC")
        assert not record_processing.is_nucl_seq("AGFTCF")

class TestPreprocessRecords(unittest.TestCase):
    def setUp(self):
        options = helpers.get_simple_options(None, [])
        options.triggered_limit = False
        self.options = config.update_config(options)
        class DummyModule:
            def __init__(self):
                self.was_run = False
            def run_on_record(self, *_args, **_kwargs):
                self.was_run = True
        self.genefinding = DummyModule()

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
        filepath = path.get_full_path(__file__, "data", "nisin.fasta")
        records = record_processing.parse_input_sequence(filepath)
        assert len(records) == 1
        assert len(records[0].get_cds_features()) == 0
        # make sure genefinding wasn't run with default options
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert not self.genefinding.was_run
        assert not records[0].get_cds_features()

        # make sure genfinding was run when not 'none'
        records[0].skip = False
        config.update_config({"genefinding_tool" : "not-none"})
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert self.genefinding.was_run
        # still no features because we used dummy genefinding
        assert not records[0].get_cds_features()

    def test_nisin_fasta_gff(self):
        fasta = path.get_full_path(__file__, "data", "nisin.fasta")
        gff = path.get_full_path(__file__, "data", "nisin.gff3")
        config.update_config({"genefinding_gff3" : gff})
        records = record_processing.parse_input_sequence(fasta)
        record_processing.pre_process_sequences(records, self.options, self.genefinding)
        assert not self.genefinding.was_run
        assert len(records[0].get_cds_features()) == 11

    def test_shotgun(self):
        filepath = path.get_full_path(__file__, "data", "wgs.gbk")
        records = record_processing.parse_input_sequence(filepath)
        with self.assertRaisesRegex(RuntimeError,
                     "Incomplete whole genome shotgun records are not supported"):
            record_processing.pre_process_sequences(records, self.options, None)

    def test_duplicate_record_ids(self):
        records = self.read_double_nisin()
        assert records[0].id == records[1].id
        record_processing.pre_process_sequences(records, self.options, None)
        assert len(records) == 2
        assert records[0].id != records[1].id

    def test_limit(self):
        records = self.read_double_nisin()
        assert all(rec.skip == False for rec in records)
        assert not self.options.triggered_limit
        config.update_config({"limit" : 1})
        record_processing.pre_process_sequences(records, self.options, None)
        assert records[0].skip == False
        assert records[1].skip.startswith("skipping all but first 1")
        assert self.options.triggered_limit

    def test_limit_to_record_partial(self):
        records = self.read_double_nisin()
        assert all(rec.skip == False for rec in records)
        config.update_config({"limit_to_record" : records[0].id})
        records[0].id += "_changed"
        record_processing.pre_process_sequences(records, self.options, None)
        assert not records[1].skip
        assert records[0].skip.startswith("did not match filter")

    def test_limit_to_record_complete(self):
        records = self.read_double_nisin()
        config.update_config({"limit_to_record" : "bad_id"})
        with self.assertRaisesRegex(ValueError, "No sequences matched filter.*"):
            record_processing.pre_process_sequences(records, self.options, None)
