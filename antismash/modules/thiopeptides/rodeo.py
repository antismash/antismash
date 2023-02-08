# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" RODEO subsections of thiopeptide analysis """

import os
import re
from typing import List, Set, Tuple, Optional

import joblib

from antismash.common import path, utils


def acquire_rodeo_heuristics(leader: str, core: str,  # pylint: disable=too-many-branches,too-many-statements
                             domains: Set[str]) -> Tuple[int, List[float]]:
    """ Calculate heuristic scores for RODEO """
    tabs: List[float] = []
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
    motifs = (('CS{2,}', 1), ('CT{2,}', 1),
              ('S{2,}', 1), ('S{3,}', 1), ('S{4,}', 2),
              ('C{2,}', 1), ('C{3,}', 1), ('C{4,}', 2),
              ('T{2,}', 1), ('T{3,}', 1), ('T{4,}', 2))
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
    stats = ThioStatistics(core)
    if stats.cst_repeats > 4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Avg heterocycle block length > 3
    if stats.average_heteroblock_length > 3:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader net charge < 5
    charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
    leader_charge = sum(charge_dict.get(aa, 0) for aa in leader)
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
    if core[-1] in "CST":
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Core contains >= 2 positive residues
    if sum(core.count(aa) for aa in "RK") >= 2:
        score -= 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Number of heterocyclizable residues to core ratio > 0.4
    if sum(core.count(aa) for aa in "CST") / len(core) >= 0.4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    return score, tabs


def generate_rodeo_svm_csv(leader: str, core: str, previously_gathered_tabs: List[float]) -> List[float]:
    """Generates all the items for one candidate precursor peptide"""
    precursor = leader + core
    columns: List[float] = []
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns.extend(previously_gathered_tabs)
    stats = ThioStatistics(core)
    # Number repeating blocks of heterocyclizable residues in core
    columns.append(stats.block_repeats)
    # Number of core repeating Cys
    columns.append(stats.c_repeats)
    # Number of core repeating Ser
    columns.append(stats.s_repeats)
    # Number of core repeating Thr
    columns.append(stats.t_repeats)
    # Number of blocks of heterocyclizable residues in core
    columns.append(stats.heteroblocks)
    # Average core heterocycle block length
    columns.append(stats.average_heteroblock_length)
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
    columns.append(len(core) / len(leader))
    # Ratio of heterocyclizable residues / length of core
    columns.append(sum(core.count(aa) for aa in "CST") / len(core))
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # groups for aromatics, neg charged, pos charged, aliphatic, and hydroxyl
    groups = ["FWY", "DE", "RK", "RKDE", "GAVLMI", "ST"]
    # leader groups
    for group in groups:
        columns.append(sum(leader.count(aa) for aa in group))
    # Counts of AAs in core
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # core groups
    for group in groups:
        columns.append(sum(core.count(aa) for aa in group))
    # Counts of AAs in entire precursor (leader+core)
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # combined groups
    for group in groups:
        columns.append(sum(precursor.count(aa) for aa in group))
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


def run_rodeo(leader: str, core: str, domains: Set[str]) -> Tuple[bool, int]:
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv = acquire_rodeo_heuristics(leader, core, domains)
    rodeo_score += heuristic_score

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, gathered_tabs_for_csv)
    rodeo_score += run_rodeo_svm(csv_columns)

    return rodeo_score >= 20, rodeo_score


class ThioStatistics:
    """ Contains the same statistics as the old thioscout() function, but
        using lazy evaluation to only calculate when required."""
    def __init__(self, core: str) -> None:
        self._core = core
        self._c_repeats: Optional[int] = None
        self._s_repeats: Optional[int] = None
        self._t_repeats: Optional[int] = None
        self._block_repeats: Optional[int] = None
        self._heteroblocks: Optional[int] = None
        self._average_heteroblock_length: Optional[float] = None

    @property
    def c_repeats(self) -> int:
        """ The number of repeating Cys in the core """
        if self._c_repeats is None:
            self._c_repeats = len("".join(re.findall("C{2,}", self._core)))
        return self._c_repeats

    @property
    def s_repeats(self) -> int:
        """ The number of repeating Ser in the core """
        if self._s_repeats is None:
            self._s_repeats = len("".join(re.findall("S{2,}", self._core)))
        return self._s_repeats

    @property
    def t_repeats(self) -> int:
        """ The number of repeating Thr in the core """
        if self._t_repeats is None:
            self._t_repeats = len("".join(re.findall("T{2,}", self._core)))
        return self._t_repeats

    @property
    def cst_repeats(self) -> int:
        """ The total number of repeating Cys, Ser, and Thr in the core """
        return sum([self.c_repeats, self.s_repeats, self.t_repeats])

    @property
    def block_repeats(self) -> int:
        """ The number of distinct blocks of repeating Cys, Ser, and Thr """
        if self._block_repeats is None:
            self._block_repeats = len(re.findall('C{2,}|S{2,}|T{2,}', self._core))
        return self._block_repeats

    @property
    def heteroblocks(self) -> int:
        """ The number of blocks of heterocyclizable residues in the core """
        if self._heteroblocks is None:
            self._calculate_heteroblocks()
        assert isinstance(self._heteroblocks, int)
        return self._heteroblocks

    @property
    def average_heteroblock_length(self) -> float:
        """ The average length of blocks of heterocyclizable residues in the core """
        if self._average_heteroblock_length is None:
            self._calculate_heteroblocks()
        assert isinstance(self._average_heteroblock_length, float)
        return self._average_heteroblock_length

    def _calculate_heteroblocks(self) -> None:
        results = re.findall('[CST]{2,}', self._core)
        self._heteroblocks = len(results)
        self._average_heteroblock_length = 0.
        if self._heteroblocks:
            self._average_heteroblock_length = sum(len(res) for res in results) / self._heteroblocks
