# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides a thin wrapper on Jinja2's renderer for use in module's HTML
    output.
"""

import os
from typing import Any, Optional
from typing import List  # comment hints, pylint: disable=unused-import

import jinja2 as _jinja2
from jinja2 import Markup


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
        self.sidepanel_sections = []  # type: List[HTMLSection]
        self.detail_sections = []  # type: List[HTMLSection]

    def add_detail_section(self, label: str, section: Markup, class_name: str = "") -> None:
        """ Add a detail section with the given tab name, using the markup provided. """
        self.detail_sections.append(HTMLSection(label, section, class_name or self.name))

    def add_sidepanel_section(self, label: str, section: Markup, class_name: str = "") -> None:
        """ Add a sidepanel section with the given tab name, using the markup provided. """
        self.sidepanel_sections.append(HTMLSection(label, section, class_name or self.name))

    def __repr__(self) -> str:
        return "HTMLSections(%s)" % (self.name)


def _safe_selector(name: str) -> str:
    """ Returns a valid HTML-selector from the provided name

        NOTE: needs to be kept in sync with antismash-js
    """
    # . is a class separator
    # : is a state separator
    return name.replace(":", "-").replace(".", "-")


def collapser_start(target: str, level: str = "all") -> Markup:
    """ Builds the start of a collapser specific to the target. Must be matched
        with a collapser_end() call.

        Arguments:
            target: the name of the target (e.g. CDS name)
            level: the interaction level at which this collapser will be expanded,
                   possible levels are:
                all: always expanded, must be manually collapsed
                supercluster: expands when the relevant supercluster is selected
                cluster: expands when the relevant cluster is selected
                cds: expands when the relevant CDS is selected
                none: never expands automatically, must be manually expanded

        Returns:
            HTML fragments as a Markup instance
    """
    if level not in ["all", "supercluster", "cluster", "cds", "none"]:
        raise ValueError("unknown collapser level: %s" % level)
    classes = ["collapser", "collapser-target-%s" % _safe_selector(target)]
    classes.append("collapser-level-%s" % level)
    child = '<div class="collapser-content">'
    if level == "all":
        classes.append("expanded")
        child = child[:-1] + ' style="display:block;">'
    return Markup('<div class="%s"> </div>%s' % (" ".join(classes), child))


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
    unique_id = "%s-help-%d" % (name, _TOOLTIP_COUNTER)
    return Markup(('<div class="help-container">'
                   ' <div class="help-icon" data-id="{0}"></div>'
                   ' <span class="help-tooltip" id="{0}">{1}</span>'
                   '</div>').format(unique_id, text))


def switch(label: str, classname: str, id_attr: str = "", starts_on: bool = False) -> Markup:
    """ Creates a checkbox construct with the given class name on the construct
        and (along with id, if given) on the checkbox input itself

        Arguments:
            label: the text to use as a label for the switch
            classname: the class name to give to the checkbox
            id_attr: an optional id to specify, must be unique
            starts_on: whether the switch should be on by default
    """
    if id_attr:
        if not id_attr.startswith("#"):
            id_attr = "#" + id_attr
        id_attr = ' id="%s"' % id_attr
    check_attr = " checked" if starts_on else ""
    return Markup(('<div class="{0} switch-container"><span class="switch-desc">{1}</span>'
                   ' <label class="switch">'
                   '  <input class="{0}" type="checkbox"{2}{3}>'
                   '  <span class="slider"></span>'
                   ' </label>'
                   '</div>').format(classname, label, id_attr, check_attr))


class _Template:  # pylint: disable=too-few-public-methods
    """ A jinja-based template rendered with extra builtins available.

        Non-functional on its own, requires self.template to be set to a jinja Template
    """
    def __init__(self, template_dir: Optional[str] = None) -> None:
        self.template = None  # type: Optional[_jinja2.Template]
        if not template_dir:
            loader = _jinja2.BaseLoader()
        else:
            loader = _jinja2.FileSystemLoader(template_dir)
        self.env = _jinja2.Environment(loader=loader, autoescape=True,
                                       undefined=_jinja2.StrictUndefined)

    def render(self, **kwargs: Any) -> Markup:
        """ Renders the template HTML, providing any given arguments to the renderer """
        if not self.template:
            raise ValueError("attempting to render without a template")

        defaults = {
            "collapser_start": collapser_start,
            "collapser_end": collapser_end,
            "help_tooltip": help_tooltip,
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
    def __init__(self, template_file: str) -> None:
        super().__init__(os.path.dirname(template_file))
        self.template = self.env.get_template(os.path.basename(template_file))
