# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import defaultdict
import glob
import logging
import shutil
import os
from typing import Dict

from helperlibs.wrappers.io import TemporaryDirectory

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
    seq = cds.get_aa_sequence()
    # create input.fasta file with single query sequence to be used as input for MSA
    fasta.write_fasta([gene_id], [seq], "input" + str(inputnr) + ".fasta")
    alignment_file = alignsmcogs(smcog, inputnr)
    # Generate trimmed alignment
    trim_alignment(inputnr, alignment_file)
    # Draw phylogenetic tree
    draw_tree(inputnr)
    # Convert tree to draw PNG image
    try:
        convert_tree(inputnr, output_dir, gene_id)
    except RuntimeError:
        logging.critical("Error during image generation for %s", gene_id)  # TODO fix


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
    # stripping [ and ] because it breaks TreeGraph later on
    # stripping ( and ) because it breaks newick tree parsing
    names = [name.replace("[", "").replace("]", "").replace("(", "_").replace(")", "_") for name in list(contents)]
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


def draw_tree(inputnr) -> str:
    """ Construct phylogenetic tree with fasttree

        Returns:
            the filename of the newick tree generated
    """
    tree_filename = "tree%d.nwk" % inputnr
    command = ["fasttree", "-quiet", "-fastest", "-noml", "trimmed_alignment%d.fasta" % inputnr]
    run_result = subprocessing.execute(command, stdout=open(tree_filename, "w"))
    if not run_result.successful():
        raise RuntimeError("Fasttree failed to run successfully:", run_result.stderr)
    return tree_filename


def convert_tree(inputnr, smcog_dir, tag) -> None:
    """ Convert tree to XTG and draw PNG image using TreeGraph
    """
    core_command = ['java', '-Djava.awt.headless=true', '-jar',
                    path.get_full_path(__file__, 'external', 'TreeGraph.jar')]

    if not os.path.exists("tree%s.nwk" % inputnr):
        raise RuntimeError("No newick tree file exists for section %d" % inputnr)

    xtg_file = "tree%s.xtg" % inputnr
    tree_png = "%s.png" % tag.split('.')[0]

    command = core_command + ['-convert', 'tree%s.nwk' % inputnr, '-xtg', xtg_file]
    run_result = subprocessing.execute(command, timeout=60*20)
    if run_result.return_code or "exception" in run_result.stdout or "Exception" in run_result.stdout or not os.path.exists(xtg_file):
        raise RuntimeError("Tree conversion failed in external command attempting to create %s: %s\n%s\n%s:" % (tree_png, command, run_result.stdout, run_result.stderr))

    command = core_command + ['-image', xtg_file, tree_png]
    run_result = subprocessing.execute(command, timeout=60*20)
    if run_result.return_code or "exception" in run_result.stdout or "Exception" in run_result.stdout or not os.path.exists(tree_png):
        raise RuntimeError("Tree image creation failed in external command attempting to create %s: %s" % (tree_png, command))

    png = tag.split(".")[0] + '.png'
    target = os.path.join(smcog_dir, png)
    if os.path.exists(target):
        os.remove(target)
    shutil.move(png, smcog_dir)
    os.remove("tree%d.xtg" % inputnr)
    os.remove("trimmed_alignment%d.fasta" % inputnr)
