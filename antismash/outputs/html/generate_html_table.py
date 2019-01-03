# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles generation of a separate HTML page for knownclusterblast's MiBIG hits
"""

import os
from typing import List

from jinja2 import FileSystemLoader, Environment, StrictUndefined

from antismash.common.path import get_full_path
from antismash.modules.clusterblast.results import MibigEntry


def generate_html_table(outfile_name: str, mibig_entries: List[MibigEntry]) -> None:
    """ Generates an HTML page containing a table for MiBIG hits for CDSes

        Arguments:
            outfile_name: the path to write the HTML page to
            mibig_entries: a list of clusterblast MibigEntry hits

        Returns:
            None
    """
    if not os.path.exists(os.path.dirname(outfile_name)):
        os.mkdir(os.path.dirname(outfile_name))

    with open(outfile_name, 'w') as handle:
        env = Environment(autoescape=True, undefined=StrictUndefined,
                          loader=FileSystemLoader(get_full_path(__file__, "templates")))
        template = env.get_template('mibig_hits_table.html')

        aux = template.render(mibig_homology_file_lines=mibig_entries)
        handle.write(aux)
