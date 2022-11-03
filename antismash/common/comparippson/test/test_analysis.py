# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import glob
from io import StringIO
import json
import unittest
from unittest.mock import patch

from antismash.common.comparippson import analysis, databases
from antismash.common.comparippson.databases import (
    ComparippsonDatabase as Database,
)
from antismash.common.comparippson.data_structures import (
    build_consensus_string,
    Hit,
    Segment,
)
from antismash.common.html_renderer import RIPP_CLASSES


def create_dummy_db(name="some name", version="7.4", url="http://nowhere/@accession@",
                    id_format="@accession@", description_format="some @description@",
                    fields=None, cache=None):
    dir_name = "dummy"
    if fields is None:
        fields = ["accession", "description"]
    result = Database(name, version, url, id_format, description_format, fields, dir_name)
    if cache is not None:
        result._reference_cache = cache
    return result


def create_dummy_segment(sequence="TIVS", start=1, end=5, full_length=5):
    return Segment(sequence, start, end, full_length)


def create_dummy_hit(query_name="q_1", reference_fields=None, match_count=5,
                     query=None, reference=None):
    if query is None:
        query = create_dummy_segment()
    if reference is None:
        reference = create_dummy_segment(start=5, end=10, full_length=12)
    if reference_fields is None:
        reference_fields = {
            "accession": "ref_core_name",
            "description": "bits of data",
        }
    return Hit(query_name, reference_fields, match_count, query, reference)


class TestDBResults(unittest.TestCase):
    def test_json_conversion(self):
        dummy_db = create_dummy_db(name="conversion test")
        hits = {
            "first": [create_dummy_hit("first")],
            "second": [create_dummy_hit("second")],
        }
        results = analysis.DBResults(dummy_db, hits, {"third": "first", "fourth": "first"})
        # check there's enough variety of values stored for the test to make sense
        assert results.database.name == "conversion test"
        assert results.aliases["third"] == "first"
        assert results.hits["second"][0].match_count == hits["second"][0].match_count

        # convert to and from a JSON string
        dumped = json.dumps(results.to_json())
        rebuilt = analysis.DBResults.from_json(json.loads(dumped))

        # all those values should be the same
        assert rebuilt.database.name == "conversion test"
        assert rebuilt.aliases == results.aliases
        assert rebuilt.hits == results.hits


class TestMultiDBResults(unittest.TestCase):
    def setUp(self):
        aliases = {"third": "second"}
        self.results = analysis.MultiDBResults([
            analysis.DBResults(
                create_dummy_db(name="A"),
                {
                    "first": [
                        create_dummy_hit(query_name="first", match_count=2),
                        create_dummy_hit(query_name="first", match_count=4),
                    ],
                    "second": [create_dummy_hit(query_name="second", match_count=3)],
                },
                aliases,
            ),
            analysis.DBResults(
                create_dummy_db(name="B"),
                {
                    "first": [create_dummy_hit(query_name="first", match_count=7)],
                    "second": [
                        create_dummy_hit(query_name="second", match_count=6),
                        create_dummy_hit(query_name="second", match_count=5),
                    ],
                },
                aliases,
            ),
        ], aliases)

    def test_json(self):
        def check_instance(instance):
            assert instance.aliases["third"] == "second"
            assert [res.database.name for res in instance.db_results] == ["A", "B"]
            assert instance.db_results[1].hits["first"][0].match_count == 7

        check_instance(self.results)

        dumped = json.dumps(self.results.to_json())
        rebuilt = analysis.MultiDBResults.from_json(json.loads(dumped))

        check_instance(rebuilt)
        assert rebuilt == self.results

    def test_building_by_query(self):
        original = dict(self.results.by_query)
        assert original
        assert original
        # fake the post_init process when unset
        self.results.by_query = None
        self.results.__post_init__()
        assert self.results.by_query == original
        # fake the post_init process when it is set (and thus trusted)
        by_query = {"second": original["second"]}
        copy = dict(by_query)
        self.results.by_query = by_query
        self.results.__post_init__()
        # should have been kept without reconstruction
        assert self.results.by_query is by_query
        # and should not have been altered
        assert self.results.by_query == copy

    def test_aliases(self):
        for aliased, actual in self.results.aliases.items():
            assert aliased not in self.results.by_query
            aliased_html = self.results.get_html_for_query(aliased)
            actual_html = self.results.get_html_for_query(actual)
            # they should be different, because the requested name is used
            assert aliased_html != actual_html
            # but if that's replaced, they should be identical
            assert aliased_html.replace(aliased, actual) == actual_html

    def test_colour_subset(self):
        # S -> dha, T -> dhb
        base = self.results.get_html_for_query("first")
        assert 'class="dhb">T<' in base and len(base) > 100
        all_chars = "".join(RIPP_CLASSES)
        full = self.results.get_html_for_query("first", colour_subset=all_chars)
        assert base == full
        some = self.results.get_html_for_query("first", colour_subset="SC")
        assert 'class="dhb"' not in some
        none = self.results.get_html_for_query("first", colour_subset="")
        for val in RIPP_CLASSES.values():
            assert f'class="{val}"' not in none

    def test_empty(self):
        results = analysis.MultiDBResults([], {})
        assert results.by_query == {}
        results.get_html_for_query("some key")


class TestHits(unittest.TestCase):
    def test_parsing_single(self):
        # if this format changes, parsing will also change and the test has to update
        assert Hit.BLAST_FORMAT == "6 qacc sacc nident qseq qstart qend qlen sseq sstart send slen"
        parts = ["qname", "refname", 4, "QSEQ", 2, 5, 6, "RSEQ", 3, 7, 10]
        fields = {"test": "value", "something": "else"}
        hit = Hit.from_simple_blast_line(list(map(str, parts)), fields=dict(fields))
        assert hit.query_name == "qname"
        assert hit.query == Segment("QSEQ", 2, 5, 6)
        assert hit.reference == Segment("RSEQ", 3, 7, 10)
        assert hit.match_count == 4
        assert hit.reference_fields == fields

    def test_parsing_full(self):
        raw = "\n".join([
            "q1	1|a	32	WKSESLCTPGCVTGALQTCFLQTLTCNCKISK	1	32	32	WKSESLCTPGCVTGALQTCFLQTLTCNCKISK	1	32	32",
            "q1	2|b	6	WKSESL	1	6	32	WKSESL	5	11	32",
            "q2	1|a	30	WKSESLCTPGCVTGALQTCFLTQLTCNCKISK	1	32	32	WKSESLCTPGCVTGALQTCFLQTLTCNCKISK	1	32	32",
        ])
        entries = {
            "1": {
                "accession": "something",
                "description": "a fake",
            },
            "2": {
                "accession": "some_other_thing",
                "description": "more fakes",
            },
        }
        database = create_dummy_db(cache=entries)
        # all hits
        hits_by_query = analysis._parse_simple_hits(StringIO(raw), 0.0, database)
        assert len(hits_by_query) == 2
        best = hits_by_query["q1"][0]
        assert best.reference_fields["accession"] == "something"
        assert best.match_count == 32
        assert hits_by_query["q2"][0].match_count == 30
        assert hits_by_query["q1"][1].match_count == 6
        # and with a higher threshold, there should be less hits
        hits_by_query = analysis._parse_simple_hits(StringIO(raw), 0.5, database)
        assert len(hits_by_query) == 2
        assert hits_by_query["q1"] == [best]
        assert len(hits_by_query["q2"]) == 1

    def test_consensus(self):
        query = create_dummy_segment(sequence="TTWTVTTTGVWASTIS-NNC")
        ref = create_dummy_segment(sequence="TTWTVT---IFLSTISVNNC")
        hit = create_dummy_hit(query=query, reference=ref)
        consensus = hit.get_consensus()
        assert consensus == "TTWTVT   ++ STIS NNC"  # from blast itself


class TestConsensus(unittest.TestCase):
    def test_consensus_simple(self):
        for top, bottom, expected in ["TA ", "VI+", "AL ", "GL ", "NS+", "WG ", "NH+"]:
            assert build_consensus_string(top, bottom) == expected, f"{top} and {bottom} != {expected}"
            if expected == " ":
                assert build_consensus_string(top, bottom, space="!") == "!"


class TestDatabase(unittest.TestCase):
    def setUp(self):
        data = create_dummy_db().to_json()  # cheat a bit by building it with the class
        data.pop("dir_name")  # remove it since it's an independent input
        # then add the entries that would be found
        data["entries"] = {str(i): {k: f"{i}'s value" for k in data["fields"]} for i in range(5)}
        self.data = data

    def test_from_metadata(self):
        built = Database.from_metadata(self.data, "dummy_dir")
        assert built.name == self.data["name"]

        # make sure all keys are enforceds
        for key in self.data:
            incomplete = dict(self.data)
            incomplete.pop(key)
            with self.assertRaises(KeyError):
                Database.from_metadata(incomplete, "dummy_dir")

    def test_incomplete_entries(self):
        self.data["entries"]["1"].pop(self.data["fields"][0])
        with self.assertRaisesRegex(ValueError, "has missing fields"):
            Database.from_metadata(self.data, "dummy_dir")

    def test_entries_with_extras(self):
        key = "not_defined"
        assert key not in self.data["fields"]
        self.data["entries"]["1"][key] = "not good"
        with self.assertRaisesRegex(ValueError, "has unknown fields"):
            Database.from_metadata(self.data, "dummy_dir")


class TestGetDatabases(unittest.TestCase):
    def test_missing_files(self):
        class DummyOptions:
            thing = 1

            def __getattr__(self, key):
                self.thing += 1
                return str(self.thing)

        with patch.object(glob, "iglob", return_value=iter(["/not/a/real/path"])):
            with patch.object(databases, "find_latest_database_version"):
                with self.assertRaisesRegex(RuntimeError, "missing data"):
                    databases.get_databases(DummyOptions())


class TestFiltering(unittest.TestCase):
    def test_unique(self):
        data = {"name": "CORE"}
        cores, aliases = analysis.filter_duplicates(data)
        assert cores == data and cores is not data
        assert not aliases

    def test_duplicate(self):
        data = {"name": "CORE", "other": "CORE"}
        cores, aliases = analysis.filter_duplicates(data)
        assert cores == {"name": "CORE"}
        assert aliases == {"other": "name"}

    def test_multiple(self):
        data = {"name": "CORE", "other": "CORE", "another": "CORE"}
        cores, aliases = analysis.filter_duplicates(data)
        assert cores == {"name": "CORE"}
        assert aliases == {"other": "name", "another": "name"}

    def test_mixed(self):
        data = {"name": "CORE", "different": "DIFF", "another": "CORE"}
        cores, aliases = analysis.filter_duplicates(data)
        assert cores == {"name": "CORE", "different": "DIFF"}
        assert aliases == {"another": "name"}
