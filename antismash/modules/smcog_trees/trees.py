# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The phylogenetic tree generation section of the smcogs module.
    Uses the classifications generated elsewhere in the module to generate the
    trees based on other possible classifications.
"""

from collections import defaultdict
from io import StringIO
import logging
from typing import Dict, List, Union

from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError
import brawn
import matplotlib
from matplotlib import pyplot

from antismash.common import path, fasta, subprocessing
from antismash.common.secmet import CDSFeature
from antismash.config import get_config

# silence the matplotlib noisy logging (relevant when options.debug is set)
logging.getLogger('matplotlib').setLevel(logging.WARNING)


def generate_trees(genes_within_clusters: list[CDSFeature],
                   nrpspks_genes: List[CDSFeature]) -> Dict[str, str]:
    """ smCOG phylogenetic tree construction """
    pks_nrps_cds_names = set(feature.get_name() for feature in nrpspks_genes)
    logging.info("Calculating and drawing phylogenetic trees of cluster genes "
                 "with smCOG members")
    cds_features = []
    for cds in genes_within_clusters:
        cds_name = cds.get_name()
        if cds_name in pks_nrps_cds_names:
            continue
        if not cds.gene_functions.get_by_tool("smcogs"):
            continue
        cds_features.append(cds)

    args = []
    for cds in cds_features:
        smcog = cds.gene_functions.get_by_tool("smcogs")[0].description.split(":")[0]
        args.append([cds, smcog])
    return {name: tree for name, tree in subprocessing.parallel_function(smcog_tree_analysis, args)
            if tree}


def smcog_tree_analysis(cds: CDSFeature, smcog: str) -> tuple[str, str]:
    "run smCOG search on all gene cluster CDS features"
    gene_id = cds.get_name()
    alignment = align_smcogs(smcog, gene_id, cds.translation)
    trimmed = trim_alignment(alignment)
    tree = build_tree(trimmed, gene_id)
    return gene_id, tree


def align_smcogs(smcog: str, name: str, sequence: str) -> Dict[str, str]:
    """ Align to multiple sequence alignment, return as a dictionary """
    reference_path = path.get_full_path(__file__, "data", f"{smcog.lower()}_muscle.fasta")
    with open(reference_path, encoding="utf-8") as handle:
        reference = brawn.Alignment.from_file(handle)
    query = brawn.Alignment({name: sequence})
    return brawn.combine_alignments(query, reference).to_dict()


def trim_alignment(contents: dict[str, str]) -> dict[str, str]:
    """ remove all positions before the first and after the last position shared
        by at least a third of all sequences
    """

    def find_first_aa_position(conservations: List[Dict[str, int]], sequence_count: int) -> int:
        """ Finds the first position of a shared amino acid """
        for position, conservation in enumerate(conservations):
            aa = sorted(conservation.items(), key=lambda x: (x[1], x[0]), reverse=True)
            base, count = aa[0]
            # skip best hits that are gaps
            if base == "-":
                continue
            # check that the count is greater than required
            if count >= sequence_count / 3:
                return position
        return 0  # can't be earlier than the start

    # check all sequences are the same length
    sequence_length = len(list(contents.values())[0])
    for name, seq in contents.items():
        assert sequence_length == len(seq), f"{name} has different sequence length"
    # stripping ( and ) because it breaks newick tree parsing
    # and keeping only the last two fields (id and description)
    names = ["|".join(name.replace("(", "_").replace(")", "_").rsplit('|', 2)[-2:]) for name in list(contents)]
    seqs = list(contents.values())

    # store conservation of residues
    conservations: List[Dict[str, int]] = [defaultdict(lambda: 0) for i in range(sequence_length)]
    for seq in seqs:
        for position, base in enumerate(seq):
            conservations[position][base] += 1

    # Find first and last amino acids shared
    first_shared_amino = find_first_aa_position(conservations, len(seqs))

    conservations.reverse()
    last_shared_amino = sequence_length - find_first_aa_position(conservations, len(seqs))

    # Shorten sequences to detected conserved regions
    seqs = [seq[first_shared_amino:last_shared_amino] for seq in seqs]
    return dict(zip(names, seqs))


def build_tree(alignment: dict[str, str], name: str) -> str:
    """ Builds a newick tree from the given alignment

        Arguments:
            name: the identifier of the query within the alignment
            alignment: a mapping of identifier to aligned sequence

        Returns:
            a tree as a newick formatted string
    """
    command = [get_config().executables.fasttree, "-quiet", "-fastest", "-noml"]
    run_result = subprocessing.execute(command, stdin="".join(fasta.build_fasta(alignment)))
    if not run_result.successful():
        raise RuntimeError(f"Fasttree failed to run successfully: {run_result.stderr}")
    handle = StringIO(run_result.stdout)

    try:
        tree = Phylo.read(handle, 'newick')
    except NewickError:
        logging.debug('Invalid newick tree for %r', name)
        return ''

    # enforce a minimum distance between branches
    max_size = max(tree.distance(node) for node in tree.get_terminals())
    for clade in tree.get_nonterminals() + tree.get_terminals():
        if not clade.branch_length:
            clade.branch_length = max_size / 20
        else:
            clade.branch_length = abs(clade.branch_length) + max_size / 20

    handle = StringIO("")
    Phylo.NewickIO.write([tree], handle)
    return handle.getvalue()


def draw_tree(tree: Union[str, Phylo.BaseTree], file_path: str, identifier: str = "") -> None:
    """ Generates an image of the given tree and writes it to file

        Arguments:
            tree: the tree to draw
            file_path: the path of the output file
            identifier: an optional identifier to highlight within the drawn tree

        Returns:
            None
    """
    if isinstance(tree, str):
        handle = StringIO(tree)
        tree = Phylo.read(handle, "newick")
    matplotlib.use('Agg')

    label_colors = {}
    if identifier:
        label_colors[identifier] = "green"

    Phylo.draw(tree, do_show=False, label_colors=label_colors,
               label_func=lambda node: str(node).replace("|", " "))
    fig = pyplot.gcf()
    fig.set_size_inches(20, (tree.count_terminals() / 3))
    pyplot.axis('off')
    fig.savefig(file_path, bbox_inches='tight')
    pyplot.close(fig)
