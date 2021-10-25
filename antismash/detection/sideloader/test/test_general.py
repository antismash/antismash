# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring,too-many-public-methods

import unittest

from antismash.common.errors import AntismashInputError
from antismash.common.test.helpers import DummyCDS, DummyRecord
from antismash.detection.sideloader import general, _parse_arg

from .test_loading import GOOD_FILE


def make_record(name, cds_start=-1):
    features = []
    if cds_start > -1:
        features.append(DummyCDS(start=cds_start, end=cds_start + 3))
    return DummyRecord(features=features, record_id=name)


class TestSimple(unittest.TestCase):
    def test_good_arg(self):
        result = _parse_arg("HM219853.1:50-500")
        assert result.accession == "HM219853.1"
        assert result.start == 50
        assert result.end == 500

    def test_bad_args(self):
        for bad in ["a:a:1-5", ":1-5", "a:1-", "a:1-5-50", "a:", "a:1a-500", "a:50-1"]:
            with self.assertRaises(ValueError):
                _parse_arg(bad)

    def test_result(self):
        result = general.load_single_record_annotations([], make_record("AXC"), _parse_arg("AcC:1-50"))
        assert not result.subregions

        result = general.load_single_record_annotations([], make_record("AcC", 30), _parse_arg("AcC:1-50"))
        assert not result.protoclusters
        assert len(result.subregions) == 1
        sub = result.subregions[0]
        assert sub.tool.name == "manual"
        assert sub.start == 1
        assert sub.end == 50
        assert result.record_id == "AcC"

    def test_empty(self):
        with self.assertRaisesRegex(AntismashInputError, "area contains no complete CDS"):
            general.load_single_record_annotations([], make_record("test"), _parse_arg("test:1-50"))

    def test_cds_marker(self):
        cdses = [DummyCDS(locus_tag="needle", start=20, end=26),
                 DummyCDS(locus_tag="miss", start=90, end=96)]
        record = DummyRecord(seq="A" * 100, features=cdses)
        result = general.load_single_record_annotations([], record, None, cds_markers=["nothing"])
        assert not result.subregions

        padding = 50
        result = general.load_single_record_annotations([], record, None, cds_markers=["needle"],
                                                        cds_marker_padding=padding)
        assert len(result.subregions) == 1
        sub = result.subregions[0]
        assert sub.tool.name == "manual"
        assert sub.start == 0  # shouldn't be negative when padding applied
        assert sub.end == cdses[0].location.end + padding


class TestSingleFile(unittest.TestCase):
    def test_filtering_by_record_id(self):
        results = general.load_single_record_annotations([GOOD_FILE], make_record("HM219853.1", 250), None)
        assert len(results.subregions) == 1
        assert results.subregions[0].label == "Polyketide"
        assert len(results.protoclusters) == 1
        assert results.protoclusters[0].product == "T1PKS"

        results = general.load_single_record_annotations([GOOD_FILE], make_record("not-HM219853.1", 1250), None)
        assert len(results.subregions) == 1
        assert results.subregions[0].label == "unknown"
        assert len(results.protoclusters) == 1
        assert results.protoclusters[0].product == "NRPS"

        results = general.load_single_record_annotations([GOOD_FILE], make_record("nomatch"), None)
        assert not results.subregions
        assert not results.protoclusters

    def test_multi_file(self):
        results = general.load_single_record_annotations([GOOD_FILE, GOOD_FILE], make_record("HM219853.1", 250), None)
        assert len(results.subregions) == 2
        assert len(results.protoclusters) == 2
        for sub in results.subregions:
            assert sub.label == "Polyketide"
        for proto in results.protoclusters:
            assert proto.product == "T1PKS"
