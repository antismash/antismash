# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock

from antismash.common.secmet.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config, update_config
from antismash.outputs import gbk
from antismash.outputs.gbk import (
    add_antismash_comments,
)


class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.rec = DummyRecord()
        self.options = build_config([], modules=[gbk])
        update_config({
            "start": -1,
            "end": -1,
            "version": "5.dummy",
            "branding": "antiSMASH",
        })
        self.bio = self.rec.to_biopython()

    def tearDown(self):
        destroy_config()

    def test_empty(self):
        # it shouldn't crash is the only real test here, since no other operations will happen
        add_antismash_comments([], self.options)

    def test_unbounded(self):
        add_antismash_comments([(self.rec, self.bio)], self.options)
        comment = self.bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert comment["Version"] == self.options.version
        assert "Original ID" not in comment
        assert "Starting at" not in comment
        assert "Ending at" not in comment

    def test_start_bounded(self):
        update_config({"start": 7})
        add_antismash_comments([(self.rec, self.bio)], self.options)
        comment = self.bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert "Original ID" not in comment
        assert comment["Starting at"] == "7"

    def test_end_bounded(self):
        update_config({
            "start": -1,
            "end": 1000,
        })
        add_antismash_comments([(self.rec, self.bio)], self.options)
        comment = self.bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert "Original ID" not in comment
        assert comment["Ending at"] == "1000"

    def test_changed_record_id(self):
        self.rec.original_id = "something else"
        add_antismash_comments([(self.rec, self.bio)], self.options)
        comment = self.bio.annotations["structured_comment"]["antiSMASH-Data"]
        assert comment["Original ID"] == self.rec.original_id

    def test_branding(self):
        alt_bio = self.rec.to_biopython()
        update_config({"branding": "alternate_branding"})
        add_antismash_comments([(self.rec, self.bio)], self.options)
        alt_structured = alt_bio.annotations["structured_comment"]
        assert len(alt_structured) == 1, list(alt_structured)
        comment = alt_structured[f"{self.options.branding}-Data"]
        assert {"Version", "Run date"}.issubset(set(comment))


class TestPaths(unittest.TestCase):
    def tearDown(self):
        destroy_config()

    def build_config(self, extra_args):
        destroy_config()
        args = [
            "--output-base", "dummybase",
            "--output-dir", "dummy_dir",
        ] + extra_args
        return build_config(args, isolated=True, modules=[gbk])

    def build_summary(self, extra_args):
        options = self.build_config(extra_args)
        return gbk.get_summary_file_path(options)

    def test_summary_default(self):
        path = self.build_summary([])
        assert path.endswith("dummy_dir/dummybase.gbk")

    def test_summary_compressed(self):
        path = self.build_summary(["--gbk-compress-summary"])
        assert path.endswith("dummy_dir/dummybase.gbk.gz")

    def test_summary_disabled(self):
        for disable in ["--no-gbk-summary", "--no-enable-gbk", "--exclude-modules=gbk"]:
            for args in [[], ["--gbk-compress-summary"]]:
                path = self.build_summary([disable] + args)
                assert path is None

    def build_region(self, rec_id, region_number, extra_args):
        options = self.build_config(extra_args)
        region = Mock(get_region_number=lambda: region_number)
        return gbk.get_region_file_path(DummyRecord(record_id=rec_id), region, options)

    def test_region_default(self):
        path = self.build_summary([])
        for rec_id in ["A", "B"]:
            path = self.build_region(rec_id, 2, [])
            assert path.endswith(f"{rec_id}.region002.gbk")
            path = self.build_region(rec_id, 11, [])
            assert path.endswith(f"{rec_id}.region011.gbk")

    def test_regions_disabled(self):
        for disable in ["--no-gbk-regions", "--no-enable-gbk", "--exclude-modules=gbk"]:
            path = self.build_region("dummy_id", 0, [disable])
            assert path is None
