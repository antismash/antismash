# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import OrderedDict
from io import StringIO
import unittest

import Bio.SeqIO as seqio
from Bio.SeqFeature import ExactPosition, BeforePosition, AfterPosition, \
                           UnknownPosition, FeatureLocation, CompoundLocation, \
                           SeqFeature
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import serialiser
from antismash.common.secmet import Record
from antismash.common.test.helpers import get_path_to_nisin_genbank


class TestLocationSerialiser(unittest.TestCase):
    # TODO: test invalid inputs
    def convert(self, location, expected_type=FeatureLocation):
        assert isinstance(location, expected_type)

        before_string = str(location)
        print(before_string)  # just for help when debugging a failing test
        after_string = serialiser.location_to_json(location)
        assert isinstance(after_string, str)
        assert before_string == after_string

        new_location = serialiser.location_from_json(after_string)
        assert isinstance(new_location, expected_type)

        return new_location

    def test_before_position(self):
        location = FeatureLocation(BeforePosition(1), ExactPosition(6), strand=-1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, BeforePosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, ExactPosition)
        assert new_location.end == 6

    def test_after_position(self):
        location = FeatureLocation(ExactPosition(1), AfterPosition(6), strand=1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, ExactPosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, AfterPosition)
        assert new_location.end == 6

    def test_unknown_position(self):
        location = FeatureLocation(ExactPosition(1), UnknownPosition(), strand=1)
        new_location = self.convert(location)

        assert isinstance(new_location.start, ExactPosition)
        assert new_location.start == 1

        assert isinstance(new_location.end, UnknownPosition)

    def test_compound(self):
        first = FeatureLocation(1, 6, strand=1)
        second = FeatureLocation(10, 16, strand=1)
        location = CompoundLocation([first, second], operator="join")
        assert 5 in location
        assert 7 not in location
        assert 15 in location

        new_location = self.convert(location, expected_type=CompoundLocation)
        assert location.start == 1
        assert 5 in new_location
        assert 7 not in new_location
        assert 15 in new_location
        assert location.end == 16
        assert new_location.operator == "join"

    def test_strands(self):
        for strand in [1, 0, -1, None]:
            location = FeatureLocation(1, 6, strand=strand)
            new_location = self.convert(location)
            assert new_location.strand == strand


class TestResultsJSON(unittest.TestCase):
    def create_data_stream(self, data):
        def force_source_ordering(records):
            for record in records:
                source = None
                for feature in record.features:
                    if feature.type == "source":
                        source = feature
                        break
                new = OrderedDict()
                for qual in sorted(list(source.qualifiers)):
                    new[qual] = source.qualifiers[qual]
                source.qualifiers = new

        stream = StringIO()
        biopython = [rec.to_biopython() for rec in data]

        # because otherwise there's a nondeterministic ordering mismatch
        force_source_ordering(biopython)

        seqio.write(biopython, stream, "genbank")
        stream.seek(0)
        return stream

    def test_record_to_json_and_back(self):
        filename = get_path_to_nisin_genbank()
        records = list(seqio.parse(open(filename), "genbank"))
        records = [Record.from_biopython(rec, taxon="bacteria") for rec in records]
        rec_results = [{}, {}, {}]
        results = serialiser.AntismashResults(filename, records, rec_results, "dummy")
        json_handle = StringIO()
        results.write_to_file(json_handle)
        json_handle.seek(0)
        new_results = serialiser.AntismashResults.from_file(json_handle, taxon="bacteria")
        assert results.to_json() == new_results.to_json()
        # check no records were lost
        assert len(new_results.records) == len(results.records)
        # check that the contents of the records is the same
        #  by converting to biopython and writing to genbanks
        original = self.create_data_stream(results.records)
        new = self.create_data_stream(new_results.records)
        oldvalue = original.getvalue()
        newvalue = new.getvalue()
        with TemporaryDirectory(change=True):
            open("old.json", "w").write(oldvalue)
            open("new.json", "w").write(newvalue)
            for oldline, newline in zip(oldvalue.split('\n'), newvalue.split('\n')):
                assert oldline == newline


class TestFeatureSerialiser(unittest.TestCase):
    def test_simple_feature(self):
        location = FeatureLocation(1, 6, strand=1)
        f_type = "test type"
        qualifiers = {"a": ["1", "2"], "b": ["3", "4"]}
        f_id = "dummy id"
        # skipping biopython deprecated members: ref, ref_db, strand, location_operator

        feature = SeqFeature(location=location, type=f_type,
                             qualifiers=qualifiers, id=f_id)
        print(str(feature))

        json = serialiser.feature_to_json(feature)
        print(json)  # for debugging failures
        new_feature = serialiser.feature_from_json(json)
        print(str(new_feature))
        assert new_feature.qualifiers == feature.qualifiers
        assert new_feature.id == feature.id
        assert new_feature.type == feature.type
        assert str(new_feature.location) == str(new_feature.location)
