# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import json
import unittest

from ..secmet import _parse_format, SecMetQualifier


class TestReverseParse(unittest.TestCase):
    def test_simple(self):
        fmt = "score: {}, name: {}"
        scores = [0, 0., 204, 700.5, -250, -21.]
        names = ["thing", "stuff", "name:", "score:"]
        for score in scores:
            for name in names:
                p_score, p_name = _parse_format(fmt, fmt.format(score, name))
                assert float(p_score) == score
                assert p_name == name

    def test_specified(self):
        fmt = "evalue: {:g}"
        values = [1e-24, 0., 7e-398]
        for value in values:
            parsed = _parse_format(fmt, fmt.format(value))
            assert len(parsed) == 1
            assert float(parsed[0]) == value

    def test_literal_braces(self):
        fmt = "{{{}}}"
        formatted = fmt.format("stuff")
        assert formatted == "{stuff}"  # for sanity
        parsed = _parse_format(fmt, formatted)
        assert len(parsed) == 1
        assert parsed[0] == "stuff"

        fmt = "a {{ {} }} b"
        formatted = fmt.format("stuff")
        assert formatted == "a { stuff } b"  # for sanity
        parsed = _parse_format(fmt, formatted)
        assert len(parsed) == 1
        assert parsed[0] == "stuff"

    def test_none_found(self):
        fmt = "thing {}"
        with self.assertRaisesRegex(ValueError, "could not match format"):
            _parse_format(fmt, "hing stuff")


class TestDomain(unittest.TestCase):
    def test_construction(self):
        dom = SecMetQualifier.Domain("test name", 1e-5, 120.7, 57, "test tool")
        assert dom.query_id == "test name"
        assert dom.evalue == 1e-5
        assert dom.bitscore == 120.7
        assert dom.nseeds == 57
        assert dom.tool == "test tool"

    def test_json_conversion(self):
        old = SecMetQualifier.Domain("test name", 1e-5, 120.7, 57, "test tool")
        dump = json.dumps(old.to_json())
        new = SecMetQualifier.Domain.from_json(json.loads(dump))
        assert new.query_id == old.query_id == "test name"
        assert new.evalue == old.evalue == 1e-5
        assert new.bitscore == old.bitscore == 120.7
        assert new.nseeds == old.nseeds == 57
        assert new.tool == old.tool == "test tool"

    def test_string_conversion(self):
        old = SecMetQualifier.Domain("name test", 5e-235, 20.17, 30, "tool test")
        assert repr(old) == str(old)
        dump = str(old)
        new = SecMetQualifier.Domain.from_string(dump)
        assert new.query_id == old.query_id == "name test"
        assert new.evalue == old.evalue == 5e-235
        assert new.bitscore == old.bitscore == 20.17
        assert new.nseeds == old.nseeds == 30
        assert new.tool == old.tool == "tool test"

    def test_equality(self):
        first = SecMetQualifier.Domain("name test", 5e-235, 20.17, 30, "tool test")
        second = SecMetQualifier.Domain("name test", 5e-235, 20.17, 30, "tool test")

        assert first == second
        second.query_id = "tmp"
        assert first != second
        second.query_id = first.query_id

        assert first == second
        second.evalue = 1e-5
        assert first != second
        second.evalue = first.evalue

        assert first == second
        second.bitscore = 1.
        assert first != second
        second.bitscore = first.bitscore

        assert first == second
        second.tool = "tmp"
        assert first != second

        assert first != str(first)


class TestSecMetQualifier(unittest.TestCase):
    def setUp(self):
        self.domains = []
        self.domains.append(SecMetQualifier.Domain("test name", 1e-5, 120.7, 57, "test tool"))
        self.domains.append(SecMetQualifier.Domain("name test", 5e-235, 20.17, 30, "tool test"))

    def test_basics(self):
        qual = SecMetQualifier(self.domains)
        assert qual.domains == self.domains
        assert qual.domains is not self.domains

    def test_add_domains(self):
        qual = SecMetQualifier([])
        qual.add_domains(self.domains)
        assert qual.domains == self.domains
        qual.add_domains(self.domains)  # duplicates ignored
        assert qual.domains == self.domains
        self.domains[0].query_id = "new name"
        qual.add_domains(self.domains)  # non-duplicate now
        assert len(qual.domains) == 3
        assert qual.domain_ids == ["test name", "name test", "new name"]

    def test_biopython_suitability(self):
        # must behave as a list of strings or have conversion methods used
        qual = SecMetQualifier(self.domains)
        assert isinstance(qual, list)
        for item in qual:
            assert isinstance(item, str)
        assert len(qual) == 1
        assert qual[0] == "; ".join(map(str, self.domains))

    def test_regeneration(self):
        qual = SecMetQualifier(self.domains)
        bio = list(qual)
        new = SecMetQualifier.from_biopython(bio)
        assert list(new) == bio
        for domain in new.domains:
            assert isinstance(domain, SecMetQualifier.Domain)
        assert new.domains == qual.domains
        assert new.domain_ids == qual.domain_ids

        with self.assertRaisesRegex(ValueError, "Cannot parse qualifier"):
            SecMetQualifier.from_biopython(bio + ["something else"])
