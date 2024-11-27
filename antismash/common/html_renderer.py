# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides a thin wrapper on Jinja2's renderer for use in module's HTML
    output.
"""

import os
from typing import Any, Dict, List, Optional, Union

import jinja2 as _jinja2
from markupsafe import Markup

from antismash.config import get_config

_ANTISMASH_JS_VERSION = "0.16"  # major and minor only, always use the latest patch
# defaults for RiPP sequence classes and substitutions
RIPP_CLASSES = {
    "S": "dha",
    "T": "dhb",
    "C": "cys",
}
RIPP_SUBSTITUTIONS = {key: val.title() for key, val in RIPP_CLASSES.items()}
WILDCARD_TEMPLATE = "@!{0}!@"


class HTMLSection:  # pylint: disable=too-few-public-methods
    """ A generic container for an HTML section. Requires a label for naming
        the section, the section content (as Markup), and a class name for the
        tab element for linking related elements in different areas.

        Sidepanel content shouldn't share class_name with other sidepanel content,
        the same is true for the body details.
    """
    def __init__(self, label: str, content: Markup, class_name: str) -> None:
        self.label = label
        self.content = content
        self.class_name = _safe_selector(class_name)


class HTMLSections:
    """ A class for containing all HTML output for a module, separated into
        zero or more detail panel sections and zero or more side panel sections.

        The name given in the constructor will be used as a prefix for HTML class
        names to allow for selecting both side and detail panels with a single click.
    """
    def __init__(self, name: str) -> None:
        self.name = name
        self.sidepanel_sections: List[HTMLSection] = []
        self.detail_sections: List[HTMLSection] = []

    def add_detail_section(self, label: str, section: Markup, class_name: str = "") -> None:
        """ Add a detail section with the given tab name, using the markup provided. """
        self.detail_sections.append(HTMLSection(label, section, class_name or self.name))

    def add_sidepanel_section(self, label: str, section: Markup, class_name: str = "") -> None:
        """ Add a sidepanel section with the given tab name, using the markup provided. """
        self.sidepanel_sections.append(HTMLSection(label, section, class_name or self.name))

    def __repr__(self) -> str:
        return f"HTMLSections({self.name})"


def _safe_selector(name: str) -> str:
    """ Returns a valid HTML-selector from the provided name

        NOTE: needs to be kept in sync with antismash-js
    """
    # . is a class separator
    # : is a state separator
    return name.replace(":", "-").replace(".", "-")


def build_blastp_link(locus: str, display_text: str, translation: str = None) -> Markup:
    """ Returns a link to a BlastP search.

        Arguments:
            locus: the name of the CDS
            translation: the translation of the CDS, if the link should be
                         static and not filled in on demand by javascript

        Returns:
            an HTML fragment as a Markup instance
    """
    url = Markup(
        "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&"
        f"PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY={replace_with('translation')}&"
        "LINK_LOC=protein&PAGE_TYPE=BlastSearch"
    )
    if translation is not None:
        static_url = url.replace(replace_with("translation"), translation)
        return Markup(f'<a href="{static_url}>{display_text}</a>')
    return Markup(
        '<a class="wildcard-container" data-wildcard-attrs="href" '
        f'href="{url}" target="_blank" '
        f'data-locus="{locus}">{display_text}</a>'
    )


def get_antismash_js_version() -> str:
    """ Returns the version (excluding patch) of 'antismash.js' to use """
    return _ANTISMASH_JS_VERSION


def get_antismash_js_url() -> str:
    """ Returns the external URL for the antiSMASH javascript """
    return f"https://dl.secondarymetabolites.org/releases/as-js/{_ANTISMASH_JS_VERSION}/antismash.js"


def collapser_start(target: str, level: str = "all") -> Markup:
    """ Builds the start of a collapser specific to the target. Must be matched
        with a collapser_end() call.

        Arguments:
            target: the name of the target (e.g. CDS name)
            level: the interaction level at which this collapser will be expanded,
                   possible levels are:
                all: always expanded, must be manually collapsed
                candidate: expands when the relevant candidate_cluster is selected
                protocluster: expands when the relevant protocluster is selected
                cds: expands when the relevant CDS is selected
                none: never expands automatically, must be manually expanded

        Returns:
            HTML fragments as a Markup instance
    """
    if level not in ["all", "candidate", "protocluster", "cds", "none"]:
        raise ValueError(f"unknown collapser level: {level}")
    classes = ["collapser", f"collapser-target-{_safe_selector(target)}"]
    classes.append(f"collapser-level-{level}")
    child = '<div class="collapser-content">'
    if level == "all":
        classes.append("expanded")
        child = child[:-1] + ' style="display:block;">'
    class_string = ' '.join(classes)
    return Markup(f'<div class="{class_string}"> </div>{child}')


def spanned_sequence(sequence: str, class_mapping: Dict[str, str],
                     substitutions: Dict[str, str] = None,
                     positional_classes: dict[int, str] = None) -> Markup:
    """ Builds an HTML fragment with spans for each character, if a character is
        present in the class mapping, the span will be given that class.
        Substitutions for characters can also be provided in the same way.
        Any provided substitutions will be automatically escaped for HTML.

        Arguments:
            sequence: the sequence string to build from
            class_mapping: a dictionary mapping sequence character to HTML class
            substitutions: a dictionary mapping sequence character to another string
            positional_classes: a dictionary mapping sequence index to any additional HTML classes to insert

        Returns:
            an HTML fragment as a Markup instance
    """
    if positional_classes is None:
        positional_classes = {}
    substitutions = substitutions or {}
    for class_label in class_mapping.values():
        if not class_label.replace("-", "").isalnum():
            raise ValueError(f"invalid character in HTML class: {class_label}")
    spans = []
    for i, char in enumerate(sequence):
        char_class = " ".join((class_mapping.get(char, ""), positional_classes.get(i, ""))).strip()
        extra = ""
        if char_class:
            extra = f' class="{char_class}"'
        spans.append(f'<span{extra}>{Markup.escape(substitutions.get(char, char))}</span>')
    return Markup("".join(spans))


def coloured_ripp_sequence(sequence: str, dehydrate: bool = False, colour_subset: str = None,
                           positional_classes: dict[int, str] = None) -> Markup:
    """ Builds an HTML fragment with spans for each character, using a predefined
        set of span classes and substitutions.

        Arguments:
            sequence: the sequence string to build from
            dehydrate: whether to substitute in dehydrations (e.g. "T" -> "Dbh")
            colour_subset: a subset of normally coloured RiPP characters to restrict colouring to
            positional_classes: a dictionary mapping sequence index to any additional HTML classes to insert

        Returns:
            an HTML fragment as a Markup instance
    """
    if colour_subset is not None:
        classes = {c: RIPP_CLASSES[c] for c in colour_subset}
    else:
        classes = RIPP_CLASSES
    if not dehydrate:
        return spanned_sequence(sequence, classes, positional_classes=positional_classes)
    return spanned_sequence(sequence, classes, substitutions=RIPP_SUBSTITUTIONS,
                            positional_classes=positional_classes)


def collapser_end() -> Markup:
    """ Builds the end of a collapser, must be matched with a collapser_start() call.
    """
    return Markup("</div>")


_TOOLTIP_COUNTER = 0


def help_tooltip(text: str, name: str) -> Markup:
    """ Constructs a help icon with tooltip, each will have a unique ID generated
        based on the given name.

        Arguments:
            text: the content of the tooltip
            name: a prefix for id generation

        Returns:
            A Markup instance with the constructed HTML
    """
    global _TOOLTIP_COUNTER  # pylint: disable=global-statement
    _TOOLTIP_COUNTER += 1
    unique_id = f"{name}-help-{_TOOLTIP_COUNTER}"
    return Markup(
        '<div class="help-container">'
        f' <div class="help-icon" data-id="{unique_id}"></div>'
        f' <span class="help-tooltip" id="{unique_id}">{text}</span>'
        '</div>'
    )


def replace_with(key: str) -> Markup:
    """ Creates a HTML construct where the given key will be replaced by CDS/ORF
        data by javascript. The key provided must still be included in the
        creation/export of JSON data.

        Arguments:
            key: the name of the attribute to fill in

        Returns:
            a Markup instance with the constructed HTML
    """
    return Markup(WILDCARD_TEMPLATE.format(key))


def selected_cds_marker(name: str) -> Markup:
    """ Constructs a small marker that will be automatically displayed when the
        CDS with the given name is selected in the gene overview.

        Arguments:
            name: the name of the CDS

        Returns:
            A Markup instance with the constructed HTML
    """
    return Markup(f'<span class="cds-selected-marker cds-selected-marker-{name}" data-locus="{name}"></span>')


def docs_link(label: str, subtarget: str = "") -> Markup:
    """ Constructs a link to the documentation website specified in the antiSMASH
        config

        Arguments:
            label: the text that will be used as the link text
            subtarget: the page within the docs that the link should be targeted at,
                       e.g. "modules/t2pks"

        Returns:
            a Markup instance with the constructed HTML
    """
    base = "<a class='external-link' href='{url}' target='_blank'>{label}</a>"
    url = get_config().urls.docs_baseurl + subtarget
    return Markup(base.format(url=url, label=label))


def switch(label: str, classname: str, id_attr: str = "", starts_on: bool = False) -> Markup:
    """ Creates a checkbox construct with the given class name on the construct
        and (along with id, if given) on the checkbox input itself

        Arguments:
            label: the text to use as a label for the switch
            classname: the class name to give to the checkbox
            id_attr: an optional id to specify, must be unique
            starts_on: whether the switch should be on by default
    """
    id_attr = f' id="{id_attr}"' if id_attr else ""
    check_attr = " checked" if starts_on else ""
    return Markup(
        f'<div class="{classname} switch-container"><span class="switch-desc">{label}</span>'
        ' <label class="switch">'
        f'  <input class="{classname}" type="checkbox"{id_attr}{check_attr}>'
        '  <span class="slider"></span>'
        ' </label>'
        '</div>'
    )


class _Template:  # pylint: disable=too-few-public-methods
    """ A jinja-based template rendered with extra builtins available.

        Non-functional on its own, requires self.template to be set to a jinja Template
    """
    def __init__(self, search_path: Optional[Union[str, list[str]]] = None) -> None:
        self.template: Optional[_jinja2.Template] = None
        if not search_path:
            loader = _jinja2.BaseLoader()
        else:
            loader = _jinja2.FileSystemLoader(search_path)
        self.env = _jinja2.Environment(loader=loader, autoescape=True,
                                       undefined=_jinja2.StrictUndefined)

    def render(self, **kwargs: Any) -> Markup:
        """ Renders the template HTML, providing any given arguments to the renderer """
        if not self.template:
            raise ValueError("attempting to render without a template")

        defaults = {
            "build_blastp_link": build_blastp_link,
            "collapser_start": collapser_start,
            "collapser_end": collapser_end,
            "coloured_ripp_sequence": coloured_ripp_sequence,
            "help_tooltip": help_tooltip,
            "replace_with": replace_with,
            "selected_cds_marker": selected_cds_marker,
            "spanned_sequence": spanned_sequence,
            "switch": switch,
        }
        defaults.update(kwargs)
        return Markup(self.template.render(**defaults))


class StringTemplate(_Template):  # pylint: disable=too-few-public-methods
    """ A template renderer for string templates """
    def __init__(self, template: str) -> None:
        super().__init__()
        self.template = self.env.from_string(template)


class FileTemplate(_Template):  # pylint: disable=too-few-public-methods
    """ A template renderer for file templates """
    def __init__(self, template_file: str, extra_paths: Optional[list[str]] = None) -> None:
        if extra_paths is None:
            extra_paths = []
        super().__init__(search_path=[os.path.dirname(template_file)] + extra_paths)
        self.template = self.env.get_template(os.path.basename(template_file))
