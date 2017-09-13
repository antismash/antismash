# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from collections import OrderedDict
import glob
import logging
import shutil
import os

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path, deprecated, subprocessing

def read_fasta(filename):
    # type: (str) -> Dict[str, str]
    """ reads a fasta file into a dict: id -> sequence, returns the dict """
    ids = []
    sequence_info = [] # type: List[str]
    with open(filename, "r") as fasta:
        current_seq = []
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq = []
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-", "z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        # TODO: refactor code using this func to deal with this ValueError instead
        # raise ValueError("Fasta file contains no sequences")
    return OrderedDict(zip(ids, sequence_info))

def generate_trees(smcogs_dir, hmm_results, geneclustergenes, nrpspks_genes, options):
    pksnrpscoregenenames = set([feature.get_name() for feature in nrpspks_genes])
    #smCOG phylogenetic tree construction
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

def smcog_tree_analysis(cds, inputnr, smcog, output_dir):
    "run smCOG search on all gene cluster CDS features"
    gene_id = cds.get_name()
    seq = str(deprecated.get_aa_sequence(cds))
    #create input.fasta file with single query sequence to be used as input for MSA
    deprecated.writefasta([gene_id], [seq], "input" + str(inputnr) + ".fasta")
    alignment_file = alignsmcogs(smcog, inputnr)
    #Generate trimmed alignment
    trim_alignment(inputnr, alignment_file)
    #Draw phylogenetic tree
    drawtree(inputnr)
    #Convert tree to draw PNG image
    converttree(inputnr, output_dir, gene_id)

def alignsmcogs(smcog, inputnr):
    #Align to multiple sequence alignment, output as fasta file
    reference = path.get_full_path(__file__, "data/%s_muscle.fasta" % str(smcog).lower())
    output_filename = "muscle%d.fasta" % inputnr
    musclecommand = ["muscle", "-quiet", "-profile", "-in1", reference,
                     "-in2", "input" + str(inputnr) + ".fasta",
                     "-out", output_filename]
    result = subprocessing.execute(musclecommand)
    if result.return_code:
        raise RuntimeError("Muscle failed to run: %s, %s", musclecommand, result.stderr[-100:])
    return output_filename

def trim_alignment(inputnr, alignment_file):
    #edit muscle fasta file: remove all positions before the first and after the last position shared by >33% of all sequences
    contents = read_fasta(alignment_file)
    names = list(contents)
    seqs = list(contents.values())
    #Find first and last amino acids shared conserved >33%
    #Create list system to store conservation of residues
    conservationlist = []
    lenseqs = len(seqs[0])
    nrseqs = len(seqs)
    for i in range(lenseqs):
        conservationlist.append({i : 0 for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ-"})
    for seq in seqs:
        for a, j in enumerate(seq):
            conservationlist[a][j] += 1
    firstsharedaa = 0
    lastsharedaa = lenseqs
    #Find first amino acid shared
    first = True
    nr = 0
    for i in conservationlist:
        aa = sorted(i.items(), key=lambda x: (x[1], x[0]), reverse=True)
        if aa[0][0] != "-" and i[aa[1][0]] > (nrseqs // 3) and first:
            firstsharedaa = nr
            first = False
        nr += 1
    #Find last amino acid shared
    conservationlist.reverse()
    first = True
    nr = 0
    for i in conservationlist:
        aa = sorted(i.items(), key=lambda x: (x[1], x[0]), reverse=True)
        if aa[0][0] != "-" and i[aa[1][0]] > (nrseqs // 3) and first:
            lastsharedaa = lenseqs - nr
            first = False
        nr += 1
    #Shorten sequences to detected conserved regions
    seqs = [seq[firstsharedaa:lastsharedaa] for seq in seqs]
    seq_len = len(seqs[0])
    for name, seq in zip(names, seqs):
        assert seq_len == len(seq), "%s has longer sequence" %name
    seedfastaname = "trimmed_alignment" + str(inputnr) + ".fasta"
    deprecated.writefasta(names, seqs, seedfastaname)

def drawtree(inputnr):
    #Draw phylogenetic tree with fasttree 2.1.1
    tree_filename = "tree%d.nwk" % inputnr
    command = ["fasttree", "-quiet", "-fastest", "-noml", "trimmed_alignment%d.fasta" % inputnr]
    run_result = subprocessing.execute(command, stdout=open(tree_filename, "w"))
    if not run_result.successful():
        raise RuntimeError("Fasttree failed to run successfully:", run_result.stdout)
    return tree_filename

def converttree(inputnr, smcog_dir, tag):
     #Convert tree to XTG and draw PNG image using TreeGraph
     core_command = ['java', '-Djava.awt.headless=true', '-jar', path.get_full_path(__file__, 'external/TreeGraph.jar')]

     command = core_command + ['-convert', 'tree%s.nwk'% inputnr, '-xtg', 'tree%s.xtg' % inputnr]
     run_result = subprocessing.execute(command, timeout=60*20)
     if run_result.return_code or "exception" in run_result.stderr or "Exception" in run_result.stderr:
        raise RuntimeError("Tree conversion failed in external command: %s" % command)

     command = core_command + ['-image', 'tree%s.xtg'% inputnr, "%s.png" % tag.split('.')[0]]
     run_result = subprocessing.execute(command, timeout=60*20)
     if run_result.return_code or "exception" in run_result.stderr or "Exception" in run_result.stderr:
        raise RuntimeError("Tree image creation failed in external command: %s" % command)

     shutil.copy(tag.split(".")[0] + '.png', smcog_dir)
     os.remove(tag.split(".")[0] + ".png")
     os.remove("tree%d.xtg" % inputnr)
     os.remove("trimmed_alignment%d.fasta" % inputnr)
