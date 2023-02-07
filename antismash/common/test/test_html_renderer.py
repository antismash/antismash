# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import re
import unittest

from antismash.common import html_renderer as renderer
from antismash.config import build_config


def _verify_html_tags_match(html):
    regex = re.compile(r"<(/?)(.*?)>")
    matches = regex.findall(html)
    self_closing_tags = {"br", "hr", "img", "input", "link", "meta", "embed"}
    stack = []
    print(html)
    assert matches
    for closing, tag in matches:
        tag = tag.split()[0]
        if tag in self_closing_tags:
            continue
        print(closing, tag)
        if closing:
            assert stack and stack[-1] == tag, f"extra </{tag}>"
            print("closed", stack[-1])
            stack.pop()
        else:
            stack.append(tag)
            print("opened", stack[-1])
    print(stack)
    assert not stack, f'unclosed: {", ".join(stack)}'


class TestHTMLMatcher(unittest.TestCase):
    def test_simple(self):
        _verify_html_tags_match("<div>thing</div>")
        with self.assertRaises(AssertionError):
            _verify_html_tags_match("<div>thing")
        with self.assertRaises(AssertionError):
            _verify_html_tags_match("thing</div>")

    def test_self_closed(self):
        cases = ["<br>", "<hr>", "<img src='stuff'>", "<input stuff>", "<link>",
                 "<meta>", "<embed>"]
        for html in cases:
            _verify_html_tags_match(html)


class TestSafeSelector(unittest.TestCase):
    def test_no_replacement(self):
        assert renderer._safe_selector("test") == "test"

    def test_period(self):
        assert renderer._safe_selector("some_name.3") == "some_name-3"

    def test_colon(self):
        assert renderer._safe_selector("some:name-3") == "some-name-3"

    def test_combo(self):
        assert renderer._safe_selector("some:name.3") == "some-name-3"


class TestCollapser(unittest.TestCase):
    def test_simple(self):
        for level in ["all", "candidate", "protocluster", "cds", "none"]:
            html = str(renderer.collapser_start("dummy", level)) + str(renderer.collapser_end())
            _verify_html_tags_match(html)
            assert "collapser-target-dummy" in html
            assert f"collapser-level-{level}" in html

    def test_safe_classes(self):
        html = str(renderer.collapser_start("some:name.3", "all")) + str(renderer.collapser_end())
        _verify_html_tags_match(html)
        assert "some:name.3" not in html
        assert "collapser-target-some-name-3" in html

    def test_bad_level(self):
        with self.assertRaisesRegex(ValueError, "unknown collapser level"):
            renderer.collapser_start("dummy", "invalid")


class TestDocsLink(unittest.TestCase):
    def setUp(self):
        self.url = build_config([]).urls.docs_baseurl

    def test_no_subtarget(self):
        result = renderer.docs_link("label")
        assert isinstance(result, renderer.Markup)
        expected = f"<a class='external-link' href='{self.url}' target='_blank'>label</a>"
        assert str(result) == expected

    def test_subtarget(self):
        target = "modules/modulename"
        result = renderer.docs_link("label", target)
        assert isinstance(result, renderer.Markup)
        url = self.url + target
        expected = f"<a class='external-link' href='{url}' target='_blank'>label</a>"
        assert str(result) == expected


class TestSequences(unittest.TestCase):
    def test_unclassed(self):
        result = renderer.spanned_sequence("ABCD", {})
        assert result == "".join(f"<span>{char}</span>" for char in "ABCD")

    def test_classed(self):
        result = renderer.spanned_sequence("ABA", {
            "A": "aa",
            "B": "bbbb",
        })
        a_chunk = '<span class="aa">A</span>'
        b_chunk = '<span class="bbbb">B</span>'
        assert result == f"{a_chunk}{b_chunk}{a_chunk}"

    def test_mixed(self):
        result = renderer.spanned_sequence("ABCD", {
            "A": "aa",
            "D": "some-class",
        })
        a_chunk = '<span class="aa">A</span>'
        b_chunk = '<span>B</span>'
        c_chunk = '<span>C</span>'
        d_chunk = '<span class="some-class">D</span>'
        assert result == f"{a_chunk}{b_chunk}{c_chunk}{d_chunk}"

    def test_invalid(self):
        for char in "<> ?!+": # indicative, not exhaustive
            with self.assertRaisesRegex(ValueError, "invalid character in HTML class"):
                renderer.spanned_sequence("A", {"A": f"a{char}b"})

    def test_substitutions(self):
        result = renderer.spanned_sequence("ABC", {}, substitutions={
            "A": "123",
            "C": "456",
        })
        a_chunk = '<span>123</span>'
        b_chunk = '<span>B</span>'
        c_chunk = '<span>456</span>'
        assert result == f"{a_chunk}{b_chunk}{c_chunk}"

    def test_ripps_dehydration(self):
        result = renderer.coloured_ripp_sequence("TISC", dehydrate=True)
        assert result == "".join([
            '<span class="dhb">Dhb</span>',
            '<span>I</span>',
            '<span class="dha">Dha</span>',
            '<span class="cys">Cys</span>',
        ])

    def test_ripp_clean(self):
        result = renderer.coloured_ripp_sequence("TISC", dehydrate=False)
        assert result == "".join([
            '<span class="dhb">T</span>',
            '<span>I</span>',
            '<span class="dha">S</span>',
            '<span class="cys">C</span>',
        ])
        # and that dehydration is not the default
        assert result == renderer.coloured_ripp_sequence("TISC")


class TestSwitch(unittest.TestCase):
    def get_input_part(self, html):
        _verify_html_tags_match(html)
        start = html.find("<input ")
        end = html[start:].find(">")
        input_part = html[start + 1:start + end]
        assert input_part
        return input_part.strip()

    def test_basic(self):
        html = renderer.switch("Description", "some-switch-class")
        assert '<input class="some-switch-class"' in html
        assert "Description" in html

    def test_ids(self):
        for test_id in ["some-unique-id", "other"]:
            html = renderer.switch("Desc", "some-class", id_attr=test_id)
            assert f' id="{test_id}"' in self.get_input_part(html)

    def test_enabled(self):
        html = renderer.switch("Desc", "some-class", starts_on=False)
        part = self.get_input_part(html)
        assert not part.endswith(" checked")
        assert "checked" not in part

        html = renderer.switch("Desc", "some-class", starts_on=True)
        part = self.get_input_part(html)
        assert part.endswith(" checked")


class TestSelectedMarker(unittest.TestCase):
    def test_basic(self):
        for name in ["LOCUS_12345", "orfA"]:
            html = renderer.selected_cds_marker(name)
            assert f'data-locus="{name}"' in html  # for pulling names out of the element
            assert f'selected-marker-{name}"' in html  # for finding elements by name


class TestWildcards(unittest.TestCase):
    def test_trivial(self):
        assert renderer.replace_with("foo") == "@!foo!@"


class TestBlastLink(unittest.TestCase):
    def test_static(self):
        translation = "METMG"
        html = renderer.build_blastp_link("locus", "link text", translation=translation)
        assert f"QUERY={translation}&" in html

    def test_wildcarded(self):
        for name in ["LOCUS_12345", "orfA"]:
            html = renderer.build_blastp_link(name, "link text")
            # check all the information needed to ...
            assert "wildcard-container" in html  # find the element
            assert 'data-wildcard-attrs="href"' in html  # find the attribute to update
            assert f'data-locus="{name}"' in html  # with the right ORF data
            assert f"QUERY={renderer.replace_with('translation')}&" in html  # in the right attribute

    def test_link_text(self):
        for text in ["link description", "other desc."]:
            html = renderer.build_blastp_link("locus", text)
            assert html.endswith(f">{text}</a>")


class TestJS(unittest.TestCase):
    def test_version_naming(self):
        # release candidates may exist at some point, but they shouldn't be merged
        assert float(renderer.get_antismash_js_version())

    def test_url_has_version(self):
        assert renderer.get_antismash_js_version() in renderer.get_antismash_js_url()
