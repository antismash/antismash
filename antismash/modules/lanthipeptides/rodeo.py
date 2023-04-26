# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The RODEO subsection of lanthipeptide analysis. Runs a support vector machine
    using a precalculated training set and a variety of feature attributes to
    determine whether a CDS likely contains a lanthipeptide precursor.

    Requires fimo binary to be present on the system somewhere for full accuracy,
    but not required.
"""


import os
import re
from tempfile import NamedTemporaryFile
from typing import Dict, List, Set, Tuple

import joblib

from antismash.common import path, subprocessing, secmet, utils
from antismash.config import get_config as get_global_config

from .config import get_config as get_lanthi_config


def run_rodeo(record: secmet.Record, query: secmet.CDSFeature, leader: str,
              core: str, domains: List[str]) -> int:
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv = acquire_rodeo_heuristics(record, query, leader, core, domains)
    rodeo_score += heuristic_score

    fimo_scores: Dict[int, float] = {}

    if not get_global_config().without_fimo and get_lanthi_config().fimo_present:
        # Find motifs
        fimo_scores = identify_lanthi_motifs(leader, core)

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, gathered_tabs_for_csv, fimo_scores)
    svm_score = run_rodeo_svm(csv_columns)
    rodeo_score += svm_score
    return rodeo_score


def run_rodeo_svm(csv_columns: List[float]) -> int:
    """Run RODEO SVM"""
    classifier_path = path.get_full_path(__file__, "data", "lanthipeptide.classifier.pkl")
    scaler_path = path.get_full_path(__file__, "data", "lanthipeptide.scaler.pkl")
    assert os.path.exists(classifier_path) and os.path.exists(scaler_path)
    classifier = joblib.load(classifier_path)
    scaler = joblib.load(scaler_path)
    csv_cols = [[float(i) for i in csv_columns[2:]]]
    scaled = scaler.transform(csv_cols)
    if int(classifier.predict(scaled)[0]) == 1:
        return 10
    return 0


def acquire_rodeo_heuristics(record: secmet.Record, query: secmet.CDSFeature,
                             leader: str, core: str,
                             domains: List[str]) -> Tuple[int, List[float]]:
    """ Calculate heuristic scores for RODEO

        Arguments:
            record: the record instance to analyse
            query: the feature being checked
            leader: the sequence of the peptide leader
            core: the sequence of the peptide core
            domains: the domains found within CDS features of the cluster

        Returns:
            a tuple of
                the RODEO score, and
                a list of floats for use in the RODEO SVM
    """
    tabs: List[float] = []
    score = 0
    precursor = leader + core
    # Leader peptide contains FxLD motif
    if re.search('F.LD', leader):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core residue position of Sx4C motif
    match = re.search('S....C', core)
    if match:
        tabs.append(match.span()[0])
    else:
        tabs.append(0)
    # Core residue position of Tx4C motif
    match = re.search('T....C', core)
    if match:
        tabs.append(match.span()[0])
    else:
        tabs.append(0)
    # Core residue position of Sx5C motif
    match = re.search('S.....C', core)
    if match:
        tabs.append(match.span()[0])
    else:
        tabs.append(0)
    # Core residue position of Tx5C motif
    match = re.search('T.....C', core)
    if match:
        tabs.append(match.span()[0])
    else:
        tabs.append(0)
    # Precursor is within 500 nt?
    hmmer_profiles = ['LANC_like', 'Lant_dehyd_C']
    distance = utils.distance_to_pfam(record, query, hmmer_profiles)
    if distance < 500:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains LanB dehydratase domain (PF04738)
    if "Lant_dehyd_C" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains Lan C cyclase domain (PF05147)
    if "LANC_like" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster LACKS LanB dehydratase domain (PF04738)
    if "Lant_dehyd_C" not in domains:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster LACKS Lan C cyclase domain (PF05147)
    if "LANC_like" not in domains:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains LanB dehydratase elimination C-terminal domain (PF14028)
    if "PF14028" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains S8 peptidase subtilase (PF00082)
    if "Peptidase_S8" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains C39 peptidase (PF03412)
    if "Peptidase_C39" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains ABC transporter (PF00005)
    if "PF00005" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains YcaO-like protein (PF02624)
    if "YcaO" in domains:
        score -= 4
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains ThiF-like protein (PF00899)
    if "ThiF" in domains:
        score -= 4
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains PF02052 (Gallidermin)
    if set(domains).intersection({"Gallidermin", "mature_a", "mature_b", "matura_ab"}):
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains PF8130
    if "Antimicr18" in domains:
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide mass < 4000 Da
    precursor_analysis = utils.RobustProteinAnalysis(precursor,
                                                     monoisotopic=True,
                                                     ignore_invalid=True)
    if precursor_analysis.molecular_weight() < 4000:
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Core peptide mass < 2000 Da
    core_analysis = utils.RobustProteinAnalysis(core, monoisotopic=True,
                                                ignore_invalid=True)
    if core_analysis.molecular_weight() < 2000:
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide pHMMs below:
    precursor_hit = False
    # Precursor peptide hits gallidermin superfamily (cl03420) HMM
    if cds_has_domains(query, {"TIGR03731", "Gallidermin"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits lantibio_gallid (TIGR03731) HMM
    if cds_has_domains(query, {"TIGR03731"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits lanti_SCO0268 superfamily (cl22812) HMM
    if cds_has_domains(query, {"TIGR04451", "strep_PEQAXS"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits LD_lanti_pre (TIGR04363) HMM
    if cds_has_domains(query, {"LD_lanti_pre"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits Antimicrobial18 (cl06940) HMM
    if cds_has_domains(query, {"Antimicr18"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits gallidermin (PF02052) HMM
    if cds_has_domains(query, {"Gallidermin", "mature_a", "mature_ab", "mature_b"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # precursor peptide hits Antimicrobial18 (PF08130) HMM
    if cds_has_domains(query, {"Antimicr18"}):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)

    if precursor_hit:
        score += 3

    # Precursor peptide mass (unmodified)
    precursor_analysis = utils.RobustProteinAnalysis(precursor, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(precursor_analysis.molecular_weight()))

    # Unmodified leader peptide mass
    leader_analysis = utils.RobustProteinAnalysis(leader, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(leader_analysis.molecular_weight()))

    # Unmodified core peptide mass
    core_analysis = utils.RobustProteinAnalysis(core, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(core_analysis.molecular_weight()))

    # Length of leader peptide
    tabs.append(len(leader))
    # Length of core peptide
    tabs.append(len(core))
    # Length of precursor peptide
    tabs.append(len(precursor))
    # Ratio of length of leader peptide / length of core peptide
    tabs.append(float(len(leader) / float(len(core))))
    # Core peptide ≥ 35 residues
    if len(core) >= 35:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Core peptide contains CC motif (not in last 3 residues)
    if 'CC' in core[:-3]:
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader peptide has > 4 negatively charge motifs
    if sum(leader.count(aa) for aa in "DE") > 4:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader peptide has net negative charge
    charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
    if sum(charge_dict.get(aa, 0) for aa in leader) < 0:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader residue position of FxLD motif
    match = re.search('F.LD', leader)
    if match:
        tabs.append(match.span()[0])
    else:
        tabs.append(0)
    # Core peptide contains C-terminal CC (within last 3 residues)
    if 'CC' in core[-3:]:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core peptide contains DGCGxTC / SFNS / SxxLC / CTxGC / TPGC / SFNSxC motifs
    motifs = (('DGCG.TC', 2), ('SFNS', 2),
              ('S..LC', 2), ('CT.GC', 1),
              ('TPGC', 1), ('SFNS.C', 1))
    for motif, motif_score in motifs:
        if re.search(motif, core):
            score += motif_score
            tabs.append(1)
        else:
            tabs.append(0)
    # Core peptide contains < 2 or < 3 Cys
    if core.count("C") < 2:
        score -= 6
        tabs += [1, 1]
    elif core.count("C") < 3:
        score -= 3
        tabs += [1, 0]
    else:
        tabs += [0, 0]
    # No Cys/Ser/Thr in core peptide
    for amino, penalty in [("C", -10), ("S", -4), ("T", -4)]:
        if amino not in core:
            score += penalty
            tabs.append(1)
        else:
            tabs.append(0)
    # Lanthionine regex maximum ring number > 4
    numrings, profile = lanscout(core)
    if numrings > 4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Lanthionine regex maximum ring number < 3
    if numrings < 3:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Lanthionine regex 4-membered ring/5-membered ring/6-membered ring/7-membered ring/8-membered ring
    scores = [2, 2, 2, 2, 1]
    scorepos = 0
    for ringsize in profile[:2]:
        if ringsize not in [0, 1, 2]:
            score += scores[scorepos]
            tabs.append(1)
        else:
            tabs.append(0)
        scorepos += 1
    for ringsize in profile[2:]:
        if ringsize != 0:
            score += scores[scorepos]
            tabs.append(1)
        else:
            tabs.append(0)
        scorepos += 1
    return score, tabs


def generate_rodeo_svm_csv(leader: str, core: str, previously_gathered_tabs: List[float],
                           fimo_scores: Dict[int, float]) -> List[float]:
    """Generates all the items for one candidate precursor peptide"""
    precursor = leader + core
    columns: List[float] = []
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns += previously_gathered_tabs
    # Lanthionine regex maximum ring number
    numrings, profile = lanscout(core)
    columns.append(numrings)
    # Lanthionine regex (n+4)-membered ring count
    for i in range(5):
        columns.append(profile[i])
    # Ratio of number of Cys in core peptide to sum of Ser/Thr in core peptide
    if "S" in core or "T" in core:
        columns.append(core.count("C") / float(core.count("S") + core.count("T")))
    else:
        columns.append(1.0)
    # Ratio of number of Cys/Ser/Thr to length of core peptide
    columns.append(float(core.count("S") + core.count("T") + core.count("C")) / len(core))
    # log10 p-value MEME motifs
    for i in range(1, 6):
        columns.append(fimo_scores.get(i, 0))
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    amino_groups = ["FWY", "DE", "RK", "RKDE", "GAVLMI", "ST"]
    # Number in leader of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    for group in amino_groups:
        columns.append(sum(leader.count(aa) for aa in group))
    # Number in core of each amino acid
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in core of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    for group in amino_groups:
        columns.append(sum(core.count(aa) for aa in group))
    # Number in entire precursor of each amino acid
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in entire precursor of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    for group in amino_groups:
        columns.append(sum(precursor.count(aa) for aa in group))
    return columns


def lanscout(core: str) -> Tuple[int, List[int]]:
    """ define lanthionine ring with a c-terminal Cys """
    lanlower = 2
    lanupper = 6
    lan = f"(?=([T|S].{{{lanlower},{lanupper}}}C))"
    sizes = []
    totrings = re.compile(lan, re.I).findall(core)
    for ring in totrings:
        sizes.append(len(ring))
    temp = []
    for j in range(lanlower + 2, lanupper + 3):
        temp.append(sizes.count(j))

    assert len(temp) == 5

    return len(totrings), temp


def cds_has_domains(cds: secmet.CDSFeature, domains: Set[str]) -> bool:
    """ Tests whether a cds has any of the given domains

        Arguments:
            cds: the CDSFeature to check
            domains: a set of domain names

        Returns:
            True if any of the domains are present in the CDS, otherwise False
    """
    return bool(cds.sec_met and set(cds.sec_met.domain_ids).intersection(domains))


def identify_lanthi_motifs(leader: str, core: str) -> Dict[int, float]:
    """Run FIMO to identify lanthipeptide-specific motifs"""
    motifs_file = path.get_full_path(__file__, "data", "lanthi_motifs_meme.txt")
    with NamedTemporaryFile() as tempfile:
        with open(tempfile.name, "w", encoding="utf-8") as out_file:
            out_file.write(f">query\n{leader}{core}")
        fimo_output = subprocessing.run_fimo_simple(motifs_file, tempfile.name)
    fimo_scores = {}
    for motif in fimo_output:
        fimo_scores[int(motif.pattern_name)] = motif.score
    return fimo_scores
