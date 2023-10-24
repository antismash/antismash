# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
More detailed lassopeptide analysis using HMMer-based leader peptide
cleavage sites prediction as well as prediction of number of disulfide
bridges, molecular mass and macrolactam ring.
"""

from collections import defaultdict
import logging
import re
import os
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from helperlibs.wrappers.io import TemporaryFile
import joblib

from antismash.common import (
    all_orfs,
    comparippson,
    module_results,
    path,
    subprocessing,
    utils,
)
from antismash.common.secmet import Record, CDSFeature, Protocluster, Prepeptide, GeneFunction
from antismash.common.secmet.locations import location_from_string
from antismash.common.secmet.qualifiers.prepeptide_qualifiers import LassoQualifier
from antismash.config import get_config as get_global_config

from .config import get_config as get_lasso_config

# the size limits, in aminos, of a precursor CDS
MIN_PRECURSOR_LENGTH = 20
MAX_PRECURSOR_LENGTH = 100


class LassoResults(module_results.ModuleResults):
    """ Holds the results of lassopeptide analysis for a record

    """
    schema_version = 3

    def __init__(self, record_id: str, comparippson_results: comparippson.MultiDBResults = None) -> None:
        super().__init__(record_id)
        # keep new CDS features
        self._new_cds_features: Set[CDSFeature] = set()
        # keep new CDSMotifs by the gene they match to
        # e.g. self.motifs_by_locus[gene_locus] = [motif1, motif2..]
        self.motifs_by_locus: Dict[str, List[Prepeptide]] = defaultdict(list)
        # keep clusters and which genes in them had precursor hits
        # e.g. self.clusters[cluster_number] = {gene1_locus, gene2_locus}
        self.clusters: Dict[int, Set[str]] = defaultdict(set)
        self.comparippson_results: Optional[comparippson.MultiDBResults] = comparippson_results

    def to_json(self) -> Dict[str, Any]:
        cds_features = [(str(feature.location),
                         feature.get_name()) for feature in self._new_cds_features]
        motifs = {}
        for locus, locus_motifs in self.motifs_by_locus.items():
            motifs[locus] = [motif.to_json() for motif in locus_motifs]
        comparison = None
        if comparison:
            comparison = self.comparippson_results.to_json()
        return {
            "record_id": self.record_id,
            "schema_version": LassoResults.schema_version,
            "motifs": motifs,
            "new_cds_features": cds_features,
            "protoclusters": {key: list(val) for key, val in self.clusters.items()},
            "comparippson_results": self.comparippson_results.to_json() if self.comparippson_results else None,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["LassoResults"]:
        if json.get("schema_version") != LassoResults.schema_version:
            logging.warning("Discarding Lassopeptide results, schema version mismatch")
            return None
        comp_results = comparippson.MultiDBResults.from_json(json["comparippson_results"])
        results = LassoResults(json["record_id"], comparippson_results=comp_results)
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(Prepeptide.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["protoclusters"].items()}
        for location, name in json["new_cds_features"]:
            loc = location_from_string(location)
            cds = all_orfs.create_feature_from_location(record, loc, label=name)
            results.add_cds(cds)
        return results

    def add_to_record(self, record: Record) -> None:
        existing = record.get_cds_name_mapping()
        for feature in self._new_cds_features:
            # since a precursor may be found by other RiPP modules
            if feature.get_name() not in existing:
                record.add_cds_feature(feature)

        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                record.add_cds_motif(motif)

    def add_cds(self, cds: CDSFeature) -> None:
        """ Add a newly found CDS feature that will be added to the record """
        # if already added by another protocluster, don't double up
        for existing in self._new_cds_features:
            if cds.get_name() == existing.get_name():
                return
        self._new_cds_features.add(cds)


class Lassopeptide:
    """ Class to calculate and store lassopeptide information
    """
    def __init__(self, end: int, score: float, rodeo_score: int,
                 leader: str, core: str) -> None:
        self.end = end
        self.score = score
        self.rodeo_score = rodeo_score
        self._leader = str(leader)
        self._lassotype = 'Class II'
        self._core = ''
        self.core = str(core)
        self._weight = -1
        self._monoisotopic_weight = -1
        self.number_bridges = self._calculate_number_of_bridges()
        self._macrolactam = ''
        self._c_cut = ''
        self._cut_weight = -1
        self._cut_mass = -1

    @property
    def core(self) -> str:
        """ The core of the prepeptide """
        return self._core

    @core.setter
    def core(self, seq: str) -> None:
        self._core = str(seq)
        self.number_bridges = self._calculate_number_of_bridges()

    @property
    def leader(self) -> str:
        """ The leader of the prepeptide """
        return self._leader

    @leader.setter
    def leader(self, seq: str) -> None:
        self._leader = str(seq)

    @property
    def c_cut(self) -> str:
        """ The tail of the prepeptide """
        return self._c_cut

    @c_cut.setter
    def c_cut(self, ccut: str) -> None:
        self._c_cut = str(ccut)

    def __repr__(self) -> str:
        return (
            f"Lassopeptide(..{self.end}, {self.score}, {self._lassotype!r}, "
            f"{self._core!r}, {self.number_bridges}, {self._monoisotopic_weight}({self._weight}), "
            f"{self._macrolactam}, {self.c_cut})"
        )

    def _calculate_weight(self, analysis: utils.RobustProteinAnalysis) -> float:
        """ Calculate the molecular weight/monoisotopic mass from the given input
        """
        cc_mass = 2 * self.number_bridges
        mw = analysis.molecular_weight()
        bond = 18.02
        return mw + cc_mass - bond

    @property
    def monoisotopic_mass(self) -> float:
        """ Determines the weight of the core peptide """
        if not self.core:
            raise ValueError("Cannot calculate weights without a core")
        return self._calculate_weight(utils.RobustProteinAnalysis(self.core, monoisotopic=True))

    @property
    def molecular_weight(self) -> float:
        """ Determines the weight of the core peptide """
        if not self.core:
            raise ValueError("Cannot calculate weights without a core")
        return self._calculate_weight(utils.RobustProteinAnalysis(self.core, monoisotopic=False))

    @property
    def cut_mass(self) -> float:
        """ Determines the monoisotopic mass of the core peptide without tail """
        if not self.core:
            raise ValueError("Cannot calculate cut weights without a core")
        if not self.c_cut:
            return self.monoisotopic_mass
        return self._calculate_weight(utils.RobustProteinAnalysis(self.core[:-len(self.c_cut)], monoisotopic=True))

    @property
    def cut_weight(self) -> float:
        """ Determines the weight of the core peptide without tail """
        if not self.core:
            raise ValueError("Cannot calculate cut weights without a core")
        if not self.c_cut:
            return self.molecular_weight
        return self._calculate_weight(utils.RobustProteinAnalysis(self.core[:-len(self.c_cut)], monoisotopic=False))

    @property
    def macrolactam(self) -> str:
        """
        Predict the lassopeptide macrolactam ring
        """

        if not self._core:
            raise ValueError("no core to base macrolactam prediction on")

        seq = self._core[6:9]
        for i, res in enumerate(seq):
            if res in ['E', 'D']:
                self._macrolactam = self._core[:i+7]

        return self._macrolactam

    def _calculate_number_of_bridges(self) -> int:
        """
        Predict the lassopeptide number of disulfide bridges
        """

        aas = utils.RobustProteinAnalysis(self.core, monoisotopic=True).count_amino_acids()
        if aas['C'] >= 4:
            return 2
        if aas['C'] >= 2:
            return 1
        return 0

    @property
    def lasso_class(self) -> str:
        """
        Predict the lassopeptide class based on disulfide bridges
        """
        if not self._core:
            raise ValueError("no core to base class prediction on")

        if self.number_bridges == 1:
            return 'Class III'
        if self.number_bridges == 2:
            return 'Class I'
        return 'Class II'


def predict_cleavage_site(query_hmmfile: str, target_sequence: str, threshold: float
                          ) -> Tuple[Optional[int], float]:
    """
    Function extracts from HMMER the start position, end position and score
    of the HMM alignment
    """
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)
    resvec = (None, 0.)
    for res in hmmer_res:
        for hits in res:
            for hsp in hits:
                # when hmm includes 1st macrolactam residue: end-2
                if hsp.bitscore > threshold:
                    resvec = (hsp.query_end - 1, hsp.bitscore)
                    break
    return resvec


def run_cleavage_site_phmm(fasta: str, hmmer_profile: str, threshold: float) -> Tuple[Optional[int], float]:
    """Try to identify cleavage site using pHMM"""
    profile = path.get_full_path(__file__, 'data', hmmer_profile)
    return predict_cleavage_site(profile, fasta, threshold)


def run_cleavage_site_regex(fasta: str) -> Optional[int]:
    """Try to identify cleavage site using regular expressions"""
    # Regular expressions; try 1 first, then 2, etc.
    rex1 = re.compile('(Y..P.L...G.....T)')
    rex2 = re.compile('(G.....T)')
    rex3 = re.compile('(Y..P.L)')
    rex4 = re.compile('(Y..P)')

    end = -1

    # For each regular expression, check if there is a match that is <10 AA from the end
    if re.search(rex1, fasta) and len(re.split(rex1, fasta)[-1]) > 14:
        _, end = [m.span() for m in rex1.finditer(fasta)][-1]
        end -= 5
    elif re.search(rex2, fasta) and len(re.split(rex2, fasta)[-1]) > 14:
        _, end = [m.span() for m in rex2.finditer(fasta)][-1]
        end -= 5
    elif re.search(rex3, fasta) and len(re.split(rex3, fasta)[-1]) > 14:
        _, end = [m.span() for m in rex3.finditer(fasta)][-1]
        end += 5
    elif re.search(rex4, fasta) and len(re.split(rex4, fasta)[-1]) > 14:
        _, end = [m.span() for m in rex4.finditer(fasta)][-1]
        end += 7

    if end <= 0:
        return None

    return end


def is_on_same_strand_as(cluster: Protocluster, query: CDSFeature, profile_name: str) -> bool:
    """Check if a query CDS is on same strand as gene with pHMM hit"""
    for cds in cluster.cds_children:
        if not cds.sec_met:
            continue
        if query.strand != cds.strand:
            continue
        if profile_name in cds.sec_met.domain_ids:
            return True
    return False


def acquire_rodeo_heuristics(record: Record, cluster: Protocluster, query: CDSFeature,
                             leader: str, core: str) -> Tuple[int, List[Union[float, int]]]:
    """Calculate heuristic scores for RODEO"""
    tabs: List[Union[float, int]] = []
    score = 0
    # Calcd. lasso peptide mass (Da) (with Xs average out)
    core_analysis = utils.RobustProteinAnalysis(core, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(core_analysis.molecular_weight()))

    # Distance to any biosynthetic protein (E, B, C)
    hmmer_profiles = ['PF13471', 'PF00733', 'PF05402']
    distance = utils.distance_to_pfam(record, query, hmmer_profiles)
    tabs.append(distance)
    # Within 500 nucleotides of any biosynthetic protein (E, B, C)	+1
    if distance < 500:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Within 150 nucleotides of any biosynthetic protein (E, B, C)	+1
    if distance < 150:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Greater than 1000 nucleotides from every biosynthetic protein (E, B, C)	-2
    if distance > 1000:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core cluster has 2 or 4 Cys residues	+1
    if core.count("C") in [2, 4]:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader cluster is longer than core cluster	+2
    if len(leader) > len(core):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core has 7 (Glu) or 8(Glu/Asp) or 9 (Asp) membered ring possible	+1
    if 'E' in core[6:8] or 'D' in core[7:9]:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader cluster contains GxxxxxT	+3
    if re.search('(G.....T)', leader):
        score += 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Core starts with G	+2
    if core.startswith("G"):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Peptide and lasso cyclase are on same strand	+1
    if is_on_same_strand_as(cluster, query, 'PF00733'):
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader/core cluster length ratio < 2 and > 0.5	+1
    if 0.5 <= len(leader) / len(core) <= 2:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Core starts with Cys and has an even number of Cys	0
    if core.startswith("C") and core.count("C") % 2 == 0:
        score += 0
        tabs.append(1)
    else:
        tabs.append(0)
    # Core contains no Gly	-4
    if "G" not in core:
        score -= 4
        tabs.append(1)
    else:
        tabs.append(0)
    # Core has at least one aromatic residue	+1
    if set("FWY") & set(core):
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Core has at least 2 aromatic residues	+2
    if sum(core.count(aa) for aa in "FWY") >= 2:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core has odd number of Cys	-2
    if core.count("C") % 2 != 0:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader cluster contains Trp	-1
    if "W" in leader:
        score -= 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader cluster contains Lys	+1
    if "K" in leader:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader cluster has Cys	-2
    if "C" in leader:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster does not contain PF13471	-2
    if utils.distance_to_pfam(record, query, ['PF13471']) == -1 or \
       utils.distance_to_pfam(record, query, ['PF13471']) > 10000:
        score -= 2
    # Peptide utilizes alternate start codon	-1
    if not str(query.extract(record.seq)).startswith("ATG"):
        score -= 1
    return score, tabs


def identify_lasso_motifs(leader: str, core: str) -> Tuple[List[int], int, Dict[int, float]]:
    """Run FIMO to identify lasso peptide-specific motifs"""
    motif_file = path.get_full_path(__file__, 'data', "lasso_motifs_meme.txt")
    with TemporaryFile() as tempfile:
        with open(tempfile.name, "w", encoding="utf-8") as out_file:
            out_file.write(f">query\n{leader}{core}")
        fimo_output = subprocessing.run_fimo_simple(motif_file, tempfile.name)
    fimo_scores = {}
    fimo_motifs = []
    for motif in fimo_output:
        if not motif.pattern_name.isdigit():
            continue
        motif_id = int(motif.pattern_name)
        fimo_motifs.append(motif_id)
        fimo_scores[motif_id] = motif.score
    # Calculate score
    motif_score = 0
    if 2 in fimo_motifs:
        motif_score += 4
    elif fimo_motifs:
        motif_score += 2
    else:
        motif_score += -1
    return fimo_motifs, motif_score, fimo_scores


def generate_rodeo_svm_csv(record: Record, query: CDSFeature, leader: str, core: str,
                           previously_gathered_tabs: List[Union[float, int]], fimo_motifs: List[int],
                           fimo_scores: Dict[int, float]) -> List[Union[float, int]]:
    """Generates all the items for a single precursor peptide candidate"""
    columns: List[Union[float, int]] = []
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns += previously_gathered_tabs
    # cluster has PF00733?
    if utils.distance_to_pfam(record, query, ['PF00733']) == -1 or \
       utils.distance_to_pfam(record, query, ['PF00733']) > 10000:
        columns.append(0)
    else:
        columns.append(1)
    # cluster has PF05402?
    if utils.distance_to_pfam(record, query, ['PF05402']) == -1 or \
       utils.distance_to_pfam(record, query, ['PF05402']) > 10000:
        columns.append(0)
    else:
        columns.append(1)
    # cluster has PF13471?
    if utils.distance_to_pfam(record, query, ['PF13471']) == -1 or \
       utils.distance_to_pfam(record, query, ['PF13471']) > 10000:
        columns.append(0)
    else:
        columns.append(1)
    # Leader has LxxxxxT motif?
    if re.search('(L.....T)', leader):
        columns.append(1)
    else:
        columns.append(0)
    # Core has adjacent identical aas (doubles)?
    if any(core[i] == core[i+1] for i in range(len(core) - 1)):
        columns.append(1)
    else:
        columns.append(0)
    # Core length (aa)
    columns.append(len(core))
    # Leader length (aa)
    columns.append(len(leader))
    # Precursor length (aa)
    columns.append(len(leader) + len(core))
    # Leader/core ratio
    columns.append(len(core) / len(leader))
    # Number of Pro in first 9 aa of core?
    columns.append(core[:9].count("P"))
    # Estimated core charge
    charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
    columns.append(sum(charge_dict.get(aa, 0) for aa in core))
    # Estimated leader charge
    columns.append(sum(charge_dict.get(aa, 0) for aa in leader))
    # Estimated precursor charge
    columns.append(sum(charge_dict.get(aa, 0) for aa in leader + core))
    # Absolute value of core charge
    columns.append(abs(sum(charge_dict.get(aa, 0) for aa in core)))
    # Absolute value of leader charge
    columns.append(abs(sum(charge_dict.get(aa, 0) for aa in leader)))
    # Absolute value of precursor charge
    columns.append(abs(sum(charge_dict.get(aa, 0) for aa in leader + core)))
    # Counts of AAs in leader
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
    # Counts (0 or 1) of amino acids within first AA position of core sequence
    columns += [core[0].count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Counts of AAs in leader+core
    precursor = leader + core
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]  # Temp to work with current training CSV
    # combined groups
    for group in groups:
        columns.append(sum(precursor.count(aa) for aa in group))
    # Motifs
    columns += [1 if motif in fimo_motifs else 0 for motif in range(1, 17)]
    # Total motifs hit
    columns.append(len(fimo_motifs))
    # Motif scores
    columns += [fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)]
    # Sum of MEME scores
    columns.append(sum(fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)))
    # No Motifs?
    if not fimo_motifs:
        columns.append(1)
    else:
        columns.append(0)
    # Alternate Start Codon?
    if not str(query.extract(record.seq)).startswith("ATG"):
        columns.append(1)
    else:
        columns.append(0)
    return columns


def run_rodeo_svm(csv_columns: List[float]) -> int:
    """Run RODEO SVM"""
    classifier_path = path.get_full_path(__file__, "data", "lassopeptide.classifier.pkl")
    scaler_path = path.get_full_path(__file__, "data", "lassopeptide.scaler.pkl")
    assert os.path.exists(classifier_path) and os.path.exists(scaler_path)
    classifier = joblib.load(classifier_path)
    scaler = joblib.load(scaler_path)
    csv_cols = [[float(i) for i in csv_columns[2:]]]
    scaled = scaler.transform(csv_cols)
    if int(classifier.predict(scaled)[0]) == 1:
        return 10
    return 0


def run_rodeo(record: Record, cluster: Protocluster, query: CDSFeature, leader: str, core: str) -> Tuple[bool, int]:
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv = acquire_rodeo_heuristics(record, cluster, query, leader, core)
    rodeo_score += heuristic_score

    fimo_motifs: List[int] = []
    fimo_scores: Dict[int, float] = {}
    motif_score = 0

    if not get_global_config().without_fimo and get_lasso_config().fimo_present:
        # Incorporate motif scores
        fimo_motifs, motif_score, fimo_scores = identify_lasso_motifs(leader, core)
    rodeo_score += motif_score

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(record, query, leader, core, gathered_tabs_for_csv, fimo_motifs, fimo_scores)
    rodeo_score += run_rodeo_svm(csv_columns)

    return rodeo_score >= 15, rodeo_score


def determine_precursor_peptide_candidate(record: Record, cluster: Protocluster,
                                          query: CDSFeature, query_sequence: str) -> Optional[Lassopeptide]:
    """Identify precursor peptide candidates and split into two"""

    if not MIN_PRECURSOR_LENGTH < len(query_sequence) < MAX_PRECURSOR_LENGTH:
        return None

    # Create FASTA sequence for feature under study
    lasso_a_fasta = f">query\n{query_sequence}"

    # Run sequence against pHMM to find the cleavage site position
    end, score = run_cleavage_site_phmm(lasso_a_fasta, 'precursor_2637.hmm', -20.00)

    # If no pHMM hit, try regular expression
    if end is None:
        score = 0
        end = run_cleavage_site_regex(lasso_a_fasta)
        if end is None or end > len(query_sequence) - 3:
            end = len(query_sequence) // 2 - 5

    assert end is not None

    # Run RODEO to assess whether candidate precursor peptide is judged real
    valid, rodeo_score = run_rodeo(record, cluster, query, query_sequence[:end], query_sequence[end:])
    if not valid:
        return None

    # Determine the leader and core peptide
    leader = query_sequence[:end]
    core = query_sequence[end:]
    return Lassopeptide(end + 1, score, rodeo_score, leader, core)


def run_lassopred(record: Record, cluster: Protocluster, query: CDSFeature) -> Optional[Prepeptide]:
    """General function to predict and analyse lasso peptides"""

    # Run checks to determine whether an ORF encodes a precursor peptide
    result = determine_precursor_peptide_candidate(record, cluster, query, query.translation)
    if result is None:
        return None

    # prediction of cleavage in C-terminal based on lasso's core sequence
    c_term_hmmer_profile = 'tail_cut.hmm'
    thresh_c_hit = -7.5

    aux = result.core[(len(result.core) // 2):]
    core_a_fasta = f">query\n{aux}"

    profile = path.get_full_path(__file__, 'data', c_term_hmmer_profile)
    hmmer_res = subprocessing.run_hmmpfam2(profile, core_a_fasta)

    for res in hmmer_res:
        for hits in res:
            for seq in hits:
                if seq.bitscore > thresh_c_hit:
                    result.c_cut = aux[seq.query_start+1:]

    if result is None:
        logging.debug('%r: No C-terminal cleavage site predicted', query.get_name())
        return None

    query.gene_functions.add(GeneFunction.ADDITIONAL, "lassopeptides",
                             "predicted lassopeptide")

    return result_vec_to_motif(query, result)


def result_vec_to_motif(query: CDSFeature, result: Lassopeptide) -> Prepeptide:
    """ Converts a Lassopeptide to a Prepeptide """
    core = result.core
    tail = result.c_cut
    if tail:
        core = result.core[:-len(tail)]
    weight = result.molecular_weight
    cut_mass = result.cut_mass
    cut_weight = result.cut_weight

    product = "lassopeptide"
    name = f"{query.get_name()}_{product}"
    feature = Prepeptide(query.location, product, core, name, "lassopeptides",
                         peptide_subclass=result.lasso_class,
                         score=result.score, monoisotopic_mass=result.monoisotopic_mass,
                         molecular_weight=weight,
                         leader=result.leader, tail=tail)
    feature.detailed_information = LassoQualifier(result.rodeo_score, result.number_bridges,
                                                  result.macrolactam, cut_mass, cut_weight)
    return feature


def specific_analysis(record: Record) -> LassoResults:
    """ Runs the full lassopeptide analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated LassoResults object
    """
    results = LassoResults(record.id)
    motif_count = 0
    existing_motifs: dict[str, set[str]] = defaultdict(set)
    for cluster in record.get_protoclusters():
        if cluster.product != 'lassopeptide':
            continue
        precursor_candidates = list(cluster.cds_children)
        # Find candidate ORFs that are not yet annotated
        extra_orfs = all_orfs.find_all_orfs(record, cluster, min_length=MIN_PRECURSOR_LENGTH * 3)
        precursor_candidates.extend(extra_orfs)

        for candidate in precursor_candidates:
            motif = run_lassopred(record, cluster, candidate)
            if motif is None:
                continue
            # regardless of whether it's a duplication or not, track via protocluster
            results.clusters[cluster.get_protocluster_number()].add(candidate.get_name())

            # but if it is in multiple protoclusters, don't add it more than once
            if motif.get_name() in existing_motifs[candidate.get_name()]:
                continue

            existing_motifs[candidate.get_name()].add(motif.get_name())
            results.motifs_by_locus[candidate.get_name()].append(motif)
            motif_count += 1
            # track new CDSFeatures if found with all_orfs
            if candidate.region is None:
                results.add_cds(candidate)

    cores = {}
    for motifs in results.motifs_by_locus.values():
        for motif in motifs:
            cores[motif.get_name().rsplit("_", 1)[0]] = motif.core
    results.comparippson_results = comparippson.compare_precursor_cores(cores, get_global_config())

    logging.debug("Lassopeptide module marked %d motifs", motif_count)
    return results
