# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Provides a thin wrapper on Jinja2's renderer for use in module's HTML
    output.
"""

import os
from typing import Any, Optional

import jinja2 as _jinja2
from jinja2 import Markup


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


class _Template:
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

    def render(self, **kwargs: Any) -> str:
        """ Renders the template HTML, providing any given arguments to the renderer """
        if not self.template:
            raise ValueError("attempting to render without a template")

        defaults = {
            "collapser_start": collapser_start,
            "collapser_end": collapser_end,
            "switch": switch,
        }
        defaults.update(kwargs)
        return self.template.render(**defaults)


class StringTemplate(_Template):
    """ A template renderer for string templates """
    def __init__(self, template: str) -> None:
        super().__init__()
        self.template = self.env.from_string(template)

class FileTemplate(_Template):
    """ A template renderer for file templates """
    def __init__(self, template_file: str) -> None:
        super().__init__(os.path.dirname(template_file))
        self.template = self.env.get_template(os.path.basename(template_file))
