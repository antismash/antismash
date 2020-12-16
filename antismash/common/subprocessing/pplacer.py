# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running pplacer.
"""

from typing import Dict
from helperlibs.wrappers.io import TemporaryDirectory
from antismash.common.fasta import write_fasta

from .base import execute, get_config

def run_pplacer(query_name: str,
                alignment: Dict[str, str],
                reference_pkg: str,
                reference_alignment: str,
                reference_tree: str) -> str:
    """Function that uses the reference tree with the new alignment to place
    query domains onto reference tree.
    """
    with TemporaryDirectory(change=True):
        temp_aln = 'temp_aln.fasta'
        names, seqs = [], []
        for name in alignment:
            names.append(name)
            seqs.append(alignment[name])
        write_fasta(names, seqs, temp_aln) 
        temp_pplacer_jplace = 'temp_pplacer_jplace.jplace'
        pplacer_result = execute([get_config().executables.pplacer,
                                  "-t", reference_tree,
                                  "-r", reference_alignment,
                                  "-o", temp_pplacer_jplace,
                                  "-c", reference_pkg,
                                  temp_aln])
        if not pplacer_result.successful():
            raise RuntimeError("pplacer returned %d: %r while comparing query named %s" \
                               % (pplacer_result.return_code,
                                  pplacer_result.stderr.replace("\n", ""),
                                  query_name))
        temp_pplacer_tree = 'temp_pplacer_tree.tre'
        guppy_result = execute([get_config().executables.guppy,
                                "sing",
                                "-o", temp_pplacer_tree,
                                temp_pplacer_jplace])
        if not guppy_result.successful():
            raise RuntimeError("guppy (pplacer) returned %d: %r while comparing query named %s" \
                               % (guppy_result.return_code,
                                  guppy_result.stderr.replace("\n", ""),
                                  query_name))
        tfh = open(temp_pplacer_tree, "r")
        newick_tree = tfh.read()
        tfh.close()
    return newick_tree


def run_pplacer_version() -> str:
    """ Get the version of the pplacer binary """
    pplacer = get_config().executables.pplacer
    command = [
        pplacer,
        "--version",
    ]
    version_string = execute(command).stdout
    return version_string
