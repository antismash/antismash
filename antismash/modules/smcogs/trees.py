# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


from collections import defaultdict
from io import StringIO
import glob
import logging
import shutil
import os
from typing import Dict

from Bio import Phylo
from Bio.Phylo.NewickIO import NewickError
from helperlibs.wrappers.io import TemporaryDirectory
import matplotlib
#import pylab

from antismash.common import path, fasta, subprocessing

def generate_trees(smcogs_dir, hmm_results, geneclustergenes, nrpspks_genes, options) -> Dict[str, str]:
    """ smCOG phylogenetic tree construction """
    pksnrpscoregenenames = set([feature.get_name() for feature in nrpspks_genes])
    logging.info("Calculating and drawing phylogenetic trees of cluster genes "
                 "with smCOG members")
    with TemporaryDirectory(change=True):
        cds_features = []
        for cds in geneclustergenes:
            gene_id = cds.get_name()
            if gene_id not in pksnrpscoregenenames and hmm_results.get(gene_id):
                cds_features.append(cds)
        args = []
        for index, cds in enumerate(cds_features):
            smcog = hmm_results[cds.get_name()][0].hit_id.split(":")[0]
            args.append([cds, index, smcog, smcogs_dir])
        subprocessing.parallel_function(smcog_tree_analysis, args)

    files = glob.glob("*.png")
    tree_filenames = {}
    for filename in files:
        tag = filename.rsplit(".png", 1)[0]
        tree_filenames[tag] = filename
    return tree_filenames


def smcog_tree_analysis(cds, inputnr, smcog, output_dir) -> None:
    "run smCOG search on all gene cluster CDS features"
    gene_id = cds.get_name()
    seq = cds.translation
    # create input.fasta file with single query sequence to be used as input for MSA
    fasta.write_fasta([gene_id], [seq], "input" + str(inputnr) + ".fasta")
    alignment_file = alignsmcogs(smcog, inputnr)
    # Generate trimmed alignment
    trim_alignment(inputnr, alignment_file)
    # Draw phylogenetic tree
    draw_tree(inputnr, output_dir, gene_id)


def alignsmcogs(smcog, inputnr) -> str:
    """ Align to multiple sequence alignment, output as fasta file """
    reference = path.get_full_path(__file__, "data", "%s_muscle.fasta" % str(smcog).lower())
    output_filename = "muscle%d.fasta" % inputnr
    musclecommand = ["muscle", "-quiet", "-profile", "-in1", reference,
                     "-in2", "input" + str(inputnr) + ".fasta",
                     "-out", output_filename]
    result = subprocessing.execute(musclecommand)
    if result.return_code:
        raise RuntimeError("Muscle failed to run: %s, %s", musclecommand, result.stderr[-100:])
    return output_filename


def trim_alignment(inputnr, alignment_file) -> None:
    """ remove all positions before the first and after the last position shared
        by at least a third of all sequences
    """

    def find_first_aa_position(conservations, sequence_count):
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

    contents = fasta.read_fasta(alignment_file)
    # check all sequences are the same length
    sequence_length = len(list(contents.values())[0])
    for name, seq in contents.items():
        assert sequence_length == len(seq), "%s has different sequence length" % name
    # stripping ( and ) because it breaks newick tree parsing
    # and keeping only the last two fields (id and description)
    names = ["|".join(name.replace("(", "_").replace(")", "_").rsplit('|', 2)[-2:]) for name in list(contents)]
    seqs = list(contents.values())

    # store conservation of residues
    conservations = [defaultdict(lambda: 0) for i in range(sequence_length)]
    for seq in seqs:
        for position, base in enumerate(seq):
            conservations[position][base] += 1

    # Find first and last amino acids shared
    first_shared_amino = find_first_aa_position(conservations, len(seqs))

    conservations.reverse()
    last_shared_amino = sequence_length - find_first_aa_position(conservations, len(seqs))

    # Shorten sequences to detected conserved regions
    seqs = [seq[first_shared_amino:last_shared_amino] for seq in seqs]
    seed_fasta_name = "trimmed_alignment" + str(inputnr) + ".fasta"
    fasta.write_fasta(names, seqs, seed_fasta_name)


def draw_tree(inputnr: int, output_dir: str, tag: str) -> str:
    """ Construct a PNG for display via fasttree

        Returns:
            the filename of the image generated
    """
    matplotlib.use('Agg')
    command = ["fasttree", "-quiet", "-fastest", "-noml", "trimmed_alignment%d.fasta" % inputnr]
    run_result = subprocessing.execute(command)
    if not run_result.successful():
        raise RuntimeError("Fasttree failed to run successfully:", run_result.stderr)

    handle = StringIO(run_result.stdout)
    tree_filename = os.path.join(output_dir, tag + '.png')
    try:
        tree = Phylo.read(handle, 'newick')
    except NewickError:
        logging.debug('Invalid newick tree for %r', tag)
        return ''

    # enforce a minimum distance between branches
    max_size = max(tree.distance(node) for node in tree.get_terminals())
    for clade in tree.get_nonterminals() + tree.get_terminals():
        if not clade.branch_length:
            clade.branch_length = max_size / 20
        else:
            clade.branch_length = abs(clade.branch_length) + max_size / 20
    # change the colour of the query gene
    label_colors = {tag: 'green'}

    Phylo.draw(tree, do_show=False, label_colors=label_colors,
               label_func=lambda node: str(node).replace("|", " "))
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(20, (tree.count_terminals() / 3))
    matplotlib.pyplot.axis('off')
    fig.savefig(os.path.join(output_dir, tag + '.png'), bbox_inches='tight')
    matplotlib.pyplot.close(fig)
    return tree_filename
