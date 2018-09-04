# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from ..secmet import _parse_format


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
