# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os

from jinja2 import FileSystemLoader, Environment, StrictUndefined

def generate_html_table(outfile_name, mibig_entries):
    if not os.path.exists(os.path.dirname(outfile_name)):
        os.mkdir(os.path.dirname(outfile_name))

    with open(outfile_name, 'w') as handle:
        env = Environment(autoescape=True, undefined=StrictUndefined,
                          loader=FileSystemLoader(os.path.dirname(__file__)))
        template = env.get_template('html_table.html')

        aux = template.render(mibig_homology_file_lines=mibig_entries)
        handle.write(aux)
