# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" RODEO subsections of thiopeptide analysis """

import os
import re
from typing import List

import numpy as np
from sklearn.externals import joblib

from antismash.common import path, utils


def acquire_rodeo_heuristics(leader, core, domains):
    """Calculate heuristic scores for RODEO"""
    tabs = []
    score = 0
    # Contains TOMM YcaO (PF02624)
    if "YcaO" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains LanB N-terminal domain (PF04738)
    if "Lant_dehyd_N" in domains or "PF04738" in domains or "tsrC" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains LanB C-terminal domain (PF14028)
    if "Lant_dehyd_C" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains TOMM dehydrogenase (PF00881)
    if "PF00881" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains rSAM methyltransferase (PF04055)
    if "PF04055" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains P450 (PF00067)
    if "p450" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Contains ABC transporter or abhydrolase
    abc_transp_abhydrolases = ["PF00005", "PF01061", "PF12698", "PF12697", "PF00561"]
    for dom in abc_transp_abhydrolases:
        if dom in domains:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
    # CSS/CTT, SS/SSS/SSS, CC/CCC/CCCC, TT/TT/TTTT motifs
    motifs = (('[C][S]{2,}', 1), ('[C][T]{2,}', 1), ('[S]{2,}', 1), ('[S]{3,}', 1),
              ('[S]{4,}', 2), ('[C]{2,}', 1), ('[C]{3,}', 1), ('[C]{4,}', 2),
              ('[T]{2,}', 1), ('[T]{3,}', 1), ('[T]{4,}', 2))
    for motif in motifs:
        if re.search(motif[0], core):
            score += motif[1]
            tabs.append(1)
        else:
            tabs.append(0)
    # No Cys/Ser/Thr core residues
    for aa in "CST":
        if aa not in core:
            score -= 2
            tabs.append(1)
        else:
            tabs.append(0)
    # Mass of core peptide (unmodified) < 2100
    core_analysis = utils.RobustProteinAnalysis(core, monoisotopic=True)
    if core_analysis.molecular_weight() < 2100:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Sum of repeating Cys/Ser/Thr > 4
    number_of_repeating_CST, _, avg_heteroblock_length, _ = thioscout(core)
    if sum([int(nr) for nr in number_of_repeating_CST.split(", ")]) > 4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Avg heterocycle block length > 3
    if not np.isnan(avg_heteroblock_length) and avg_heteroblock_length > 3:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader net charge < 5
    charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
    leader_charge = sum([charge_dict[aa] for aa in leader if aa in charge_dict])
    if leader_charge < 5:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader net charge > 0
    if leader_charge > 0:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader contains a Cys
    if "C" in leader:
        score -= 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Peptide terminates Cys/Ser/Thr
    if core[-1] in ["C", "S", "T"]:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Core contains >= 2 positive residues
    if sum([core.count(aa) for aa in "RK"]) >= 2:
        score -= 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Number of heterocyclizable residues to core ratio > 0.4
    if float(sum([core.count(aa) for aa in "CST"])) / len(core) >= 0.4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    return score, tabs


def generate_rodeo_svm_csv(leader, core, previously_gathered_tabs):
    """Generates all the items for one candidate precursor peptide"""
    precursor = leader + core
    columns = []
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns += previously_gathered_tabs
    # Number repeating blocks of heterocyclizable residues in core
    number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks = thioscout(core)
    columns.append(number_of_repeat_blocks)
    # Number of core repeating Cys
    columns.append(int(number_of_repeating_CST.split(", ")[0]))
    # Number of core repeating Ser
    columns.append(int(number_of_repeating_CST.split(", ")[1]))
    # Number of core repeating Thr
    columns.append(int(number_of_repeating_CST.split(", ")[2]))
    # Number of blocks of heterocyclizable residues in core
    columns.append(number_of_heteroblocks)
    # Average core heterocycle block length
    if np.isnan(avg_heteroblock_length):
        columns.append(0)
    else:
        columns.append(avg_heteroblock_length)
    # Precursor peptide mass (unmodified)
    columns.append(utils.RobustProteinAnalysis(leader+core, monoisotopic=True).molecular_weight())
    # Unmodified leader peptide mass
    columns.append(utils.RobustProteinAnalysis(leader, monoisotopic=True).molecular_weight())
    # Unmodified core peptide mass
    columns.append(utils.RobustProteinAnalysis(core, monoisotopic=True).molecular_weight())
    # Length of Precursor
    columns.append(len(precursor))
    # Length of Leader
    columns.append(len(leader))
    # Length of Core
    columns.append(len(core))
    # Ratio of length of leader / length of core
    columns.append(float(len(core)) / float(len(leader)))
    # Ratio of heterocyclizable  residues / length of core
    columns.append(float(sum([core.count(aa) for aa in "CST"])) / len(core))
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Aromatics in leader
    columns.append(sum([leader.count(aa) for aa in "FWY"]))
    # Neg charged in leader
    columns.append(sum([leader.count(aa) for aa in "DE"]))
    # Pos charged in leader
    columns.append(sum([leader.count(aa) for aa in "RK"]))
    # Charged in leader
    columns.append(sum([leader.count(aa) for aa in "RKDE"]))
    # Aliphatic in leader
    columns.append(sum([leader.count(aa) for aa in "GAVLMI"]))
    # Hydroxyl in leader
    columns.append(sum([leader.count(aa) for aa in "ST"]))
    # Counts of AAs in core
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Aromatics in core
    columns.append(sum([core.count(aa) for aa in "FWY"]))
    # Neg charged in core
    columns.append(sum([core.count(aa) for aa in "DE"]))
    # Pos charged in core
    columns.append(sum([core.count(aa) for aa in "RK"]))
    # Charged in core
    columns.append(sum([core.count(aa) for aa in "RKDE"]))
    # Aliphatic in core
    columns.append(sum([core.count(aa) for aa in "GAVLMI"]))
    # Hydroxyl in core
    columns.append(sum([core.count(aa) for aa in "ST"]))
    # Counts of AAs in entire precursor (leader+core)
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Aromatics in precursor
    columns.append(sum([precursor.count(aa) for aa in "FWY"]))
    # Neg charged in precursor
    columns.append(sum([precursor.count(aa) for aa in "DE"]))
    # Pos charged in precursor
    columns.append(sum([precursor.count(aa) for aa in "RK"]))
    # Charged in precursor
    columns.append(sum([precursor.count(aa) for aa in "RKDE"]))
    # Aliphatic in precursor
    columns.append(sum([precursor.count(aa) for aa in "GAVLMI"]))
    # Hydroxyl in precursor
    columns.append(sum([precursor.count(aa) for aa in "ST"]))
    return columns


def run_rodeo_svm(csv_columns: List[float]) -> int:
    """Run RODEO SVM"""
    classifier_path = path.get_full_path(__file__, "data", "thiopeptide.classifier.pkl")
    scaler_path = path.get_full_path(__file__, "data", "thiopeptide.scaler.pkl")
    assert os.path.exists(classifier_path) and os.path.exists(scaler_path)
    classifier = joblib.load(classifier_path)
    scaler = joblib.load(scaler_path)
    csv_cols = [[float(i) for i in csv_columns[2:]]]
    scaled = scaler.transform(csv_cols)
    if int(classifier.predict(scaled)[0]) == 1:
        return 10
    return 0


def run_rodeo(leader, core, domains):
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv = acquire_rodeo_heuristics(leader, core, domains)
    rodeo_score += heuristic_score

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, gathered_tabs_for_csv)
    rodeo_score += run_rodeo_svm(csv_columns)

    return rodeo_score >= 20, rodeo_score


def thioscout(core):
    """ThioScout function from Chris Schwalen to count repeat blocks"""
    # rex1 repeating Cys Ser Thr residues
    rex1 = re.compile('[C]{2,}|[S]{2,}|[T]{2,}')

    # rex2 contiguous cyclizable residues
    rex2 = re.compile('[C|S|T]{2,}')

    rexout1 = re.findall(rex1, core)
    number_of_repeat_blocks = len(rexout1)

    temp = "".join(rexout1)
    number_of_repeating_CST = str([temp.count("C"), temp.count("S"), temp.count("T")]).strip("[]")

    rexout2 = re.findall(rex2, core)
    number_of_heteroblocks = len(rexout2)

    if rexout2:
        avg_heteroblock_length = np.mean([len(x) for x in rexout2])
    else:
        avg_heteroblock_length = 0.0

    return number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks
