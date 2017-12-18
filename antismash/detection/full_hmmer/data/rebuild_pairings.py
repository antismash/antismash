#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
(Re)builds an importable python dictionary from a PFAM database, mapping PFAM
name to accession (e.g. CbiA -> PF01656).

By default looks for the database 'Pfam-A.hmm' in the same directory as this file,
and outputs to a python file 'name2pfamid.py' in the current working directory.

"""

import os

HEADER = '''# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Only provides a mapping of protein name to PFAM id """

NAME_TO_PFAMID = {
'''

FOOTER = """
}
"""

def rebuild_pairings(database_filename=None, output_filename="name2pfamid.py"):
    """ (Re)builds the name mapping from the given database and writes it to the
        given output file.

        Arguments:
            database_filename: the path to the database
                        (defaults to 'Pfam-A.hmm' in the same directory as this file)
            output_filename: the path for the output file

        Returns:
            the number of mappings written to file
    """
    if not database_filename:
        database_filename = os.path.join(__file__.rsplit(os.sep, 1)[0], "Pfam-A.hmm")
    names = []
    accessions = []
    with open(database_filename, "r") as database:
        for line in database:
            value = line.split()[1].strip()
            assert value, "empty value in line: %s" % line
            if line.startswith("NAME"):
                names.append(value)
            elif line.startswith("ACC"):
                # strip any trailing info from accession: e.g. PF01656.18 -> PF01656
                accessions.append(value.split(".", 1)[0])
    assert len(names) == len(accessions), "some pairings incomplete"
    lines = []
    for name, acc in zip(names, accessions):
        lines.append('    "{}" : "{}"'.format(name, acc))
    with open(output_filename, "w") as pairings:
        pairings.write(HEADER)
        pairings.write(",\n".join(lines))
        pairings.write(FOOTER)
    return len(lines)


if __name__ == "__main__":
    print("rebuilt with %d pairings" % rebuild_pairings())
