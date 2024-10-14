# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
More detailed lanthipeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of number of lanthionine
bridges and molcular mass.
"""

from collections import defaultdict
import logging
import os
import re
from typing import Any, Dict, List, Iterable, Optional, Set

from antismash.common.signature import HmmSignature
from antismash.common import all_orfs, comparippson, path, subprocessing, module_results, utils
from antismash.common.fasta import get_fasta_from_features
from antismash.common.secmet import CDSFeature, Protocluster, GeneFunction, Prepeptide, Record, Region
from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import location_from_string
from antismash.common.secmet.qualifiers.prepeptide_qualifiers import LanthiQualifier
from antismash.config import get_config

from .rodeo import run_rodeo

KNOWN_PRECURSOR_DOMAINS = set([
    'Antimicr18',
    'Gallidermin',
    'L_biotic_A',
    'lacticin_l',
    'leader_d',
    'leader_abc',
    'leader_eh',
    'mature_ha',
    'lacticin_mat',
    'mature_b',
    'mature_d',
    'mature_ab',
    'mature_a',
    'TIGR03731',
    'LD_lanti_pre',
    'strep_PEQAXS',
])

THRESH_DICT = {'Class-I': -15,
               'Class-II': -7.3,
               'Class-III': -3.5}

# the maximal number of bases in each direction of a core enzyme that a
# precursor can be defined in
MAX_PRECURSOR_DISTANCE = 10000
# the size limits, in aminos, of a precursor CDS
MIN_PRECURSOR_LENGTH = 20
MAX_PRECURSOR_LENGTH = 80


class LanthiResults(module_results.ModuleResults):
    """ Holds the results of lanthipeptide analysis for a record

    """
    schema_version = 5

    def __init__(self, record_id: str, *, comparippson_results: comparippson.MultiDBResults = None) -> None:
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
        cds_features = [(str(feature.location), feature.get_name()) for feature in self._new_cds_features]
        motifs = {}
        for locus, locus_motifs in self.motifs_by_locus.items():
            motifs[locus] = [motif.to_json() for motif in locus_motifs]
        return {
            "record_id": self.record_id,
            "schema_version": LanthiResults.schema_version,
            "motifs": motifs,
            "new_cds_features": cds_features,
            "protoclusters": {key: list(val) for key, val in self.clusters.items()},
            "comparippson": self.comparippson_results.to_json() if self.comparippson_results else None,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["LanthiResults"]:
        if json.get("schema_version") != LanthiResults.schema_version:
            logging.warning("Discarding Lanthipeptide results, schema version mismatch")
            return None
        comparison_results = None
        if json.get("comparippson") is not None:
            comparison_results = comparippson.MultiDBResults.from_json(json["comparippson"])
        results = LanthiResults(json["record_id"], comparippson_results=comparison_results)
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(Prepeptide.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["protoclusters"].items()}
        for location, name in json["new_cds_features"]:
            cds = all_orfs.create_feature_from_location(record, location_from_string(location), label=name)
            results.add_cds(cds)
        return results

    def add_to_record(self, record: Record) -> None:
        existing = record.get_cds_name_mapping()
        for feature in self._new_cds_features:
            # since a precursor may be found by other RiPP modules
            if feature.get_name() not in existing:
                record.add_cds_feature(feature)

        motifs_added: Set[str] = set()
        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                if motif.get_name() not in motifs_added:
                    record.add_cds_motif(motif)
                    motifs_added.add(motif.get_name())

    def get_motifs_for_region(self, region: Region) -> Dict[str, List[Prepeptide]]:
        """ Given a region, return a subset of motifs_by_locus for hits within
            that region
        """
        results = {}
        for cluster in region.get_unique_protoclusters():
            for locus in self.clusters.get(cluster.get_protocluster_number(), []):
                results[locus] = self.motifs_by_locus[locus]
        return results

    def add_cds(self, cds: CDSFeature) -> None:
        """ Add a newly found CDS feature that will be added to the record """
        # if already added by another protocluster, don't double up
        for existing in self._new_cds_features:
            if cds.get_name() == existing.get_name():
                return
        self._new_cds_features.add(cds)


class PrepeptideBase:
    """ A generic prepeptide class for tracking various typical components """
    def __init__(self, end: int, score: float, rodeo_score: int = 0) -> None:
        self.end = int(end)  # cleavage site position
        self.score = float(score)  # of cleavage site
        self.rodeo_score = int(rodeo_score)
        self._leader: Optional[str] = None
        self._core = ''
        self._tail: Optional[str] = None
        self._lan_bridges = -1
        self._weight = -1
        self._monoisotopic_weight = -1
        self._alt_weights: List[float] = []
        self.core_analysis_monoisotopic: Optional[utils.RobustProteinAnalysis] = None
        self.core_analysis: Optional[utils.RobustProteinAnalysis] = None

    @property
    def core(self) -> str:
        """ The sequence of the prepeptide core """
        return self._core

    @core.setter
    def core(self, seq: str) -> None:
        self.core_analysis_monoisotopic = utils.RobustProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = utils.RobustProteinAnalysis(seq, monoisotopic=False)
        self._core = seq
        self._calculate_mw()

    @property
    def leader(self) -> Optional[str]:
        """ The sequence of the prepeptide leader """
        return self._leader

    @leader.setter
    def leader(self, seq: str) -> None:
        self._leader = seq

    def __repr__(self) -> str:
        base = "PrepeptideBase(..%s, %s, %r, %r)"
        return base % (self.end, self.score, self.rodeo_score, self._core)

    @property
    def number_of_lan_bridges(self) -> int:
        """
        function determines the number of lanthionine bridges in the core peptide
        """
        raise NotImplementedError()

    def _calculate_mw(self) -> None:
        """
        (re)calculate the monoisotopic mass and molecular weight
        """
        assert self._core
        raise NotImplementedError()

    @property
    def monoisotopic_mass(self) -> float:
        """ weight of the dehydrated core
        """
        if self._monoisotopic_weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self) -> float:
        """ Weight of the dehydrated core
        """
        if self._weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self) -> List[float]:
        """ The possible alternative weights assuming one or more of the Ser/Thr
            residues aren't dehydrated
        """
        if self._alt_weights is None:
            raise ValueError("No core to calculate weight of")
        return self._alt_weights


class CleavageSiteHit:  # pylint: disable=too-few-public-methods
    """ A simple container for storing cleavage site information """
    def __init__(self, end: int, score: float, lantype: str) -> None:
        self.end = int(end)
        self.score = float(score)
        self.lantype = lantype

    def __repr__(self) -> str:
        return f"CleavageSiteHit(end={self.end}, score={self.score}, lantype={self.lantype!r})"


class Lanthipeptide(PrepeptideBase):
    """ Calculates and stores lanthipeptide information
    """
    def __init__(self, hit: CleavageSiteHit, rodeo_score: int, leader: str, core: str) -> None:
        super().__init__(hit.end, hit.score, rodeo_score)
        self.lantype = hit.lantype
        self._aminovinyl = False
        self._chlorinated = False
        self._oxygenated = False
        self._lac = False
        self.leader = leader
        self.core = core

    def __repr__(self) -> str:
        base = "Lanthipeptide(..%s, %s, %r, %r, %s, %s(%s))"
        return base % (self.end, self.score, self.lantype, self._core,
                       self._lan_bridges, self._monoisotopic_weight, self._weight)

    @property
    def number_of_lan_bridges(self) -> int:
        """ Determines the number of lanthionine bridges in the core peptide
        """
        if not self._core:
            raise ValueError("No core to calculate bridges from")
        assert self.core_analysis is not None

        amino_counts = self.core_analysis.count_amino_acids()
        no_cys = amino_counts['C']
        no_thr_ser = amino_counts['T'] + amino_counts['S']
        self._lan_bridges = min(no_cys, no_thr_ser)
        if self._aminovinyl:
            self._lan_bridges -= 1
        return self._lan_bridges

    def _calculate_mw(self) -> None:
        """ (re)calculates the monoisotopic mass and molecular weight
        """
        assert self._core, "calculating weight without a core"
        assert self.core_analysis is not None
        assert self.core_analysis_monoisotopic is not None

        amino_counts = self.core_analysis.count_amino_acids()
        no_thr_ser = amino_counts['T'] + amino_counts['S']

        mol_mass = self.core_analysis.molecular_weight()
        mods = 18.02 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._weight = mol_mass - mods

        # every unbridged Ser or Thr might not be dehydrated
        self._alt_weights = []
        for i in range(1, no_thr_ser - amino_counts['C'] + 1):
            self._alt_weights.append(self._weight + 18.02 * i)

        monoisotopic_mass = self.core_analysis_monoisotopic.molecular_weight()
        mods = 18 * no_thr_ser
        if self._aminovinyl:
            mods += 46
        if self._chlorinated:
            mods -= 34
        if self._oxygenated:
            mods -= 16
        if self._lac:
            mods -= 2
        self._monoisotopic_weight = monoisotopic_mass - mods

    @property
    def aminovinyl_group(self) -> bool:
        """ Returns True if lanthipeptide contains an aminovinyl group
        """
        return self._aminovinyl

    @aminovinyl_group.setter
    def aminovinyl_group(self, value: bool) -> None:
        """ Sets whether lanthipeptide contains an aminovinyl group and triggers
            recalculation of the molecular weight
        """
        self._aminovinyl = value
        if self._core:
            self._calculate_mw()
            # recalculate the number of lan bridges
            self._lan_bridges = -1

    @property
    def chlorinated(self) -> bool:
        """ Returns True if lanthipeptide is chlorinated
        """
        return self._chlorinated

    @chlorinated.setter
    def chlorinated(self, value: bool) -> None:
        """ Sets whether lanthipeptide is chlorinated and triggers
            recalculation of the molecular weight
        """
        self._chlorinated = value
        if self._core:
            self._calculate_mw()

    @property
    def oxygenated(self) -> bool:
        """ Returns True if lanthipeptide is oxygenated
        """
        return self._oxygenated

    @oxygenated.setter
    def oxygenated(self, value: bool) -> None:
        """ Sets whether lanthipeptide is oxygenated and triggers
            recalculation of the molecular weight
        """
        self._oxygenated = value
        if self._core:
            self._calculate_mw()

    @property
    def lactonated(self) -> bool:
        """ Returns True if lanthipeptide starts with a lactone
        """
        return self._lac

    @lactonated.setter
    def lactonated(self, value: bool) -> None:
        """ Sets whether lanthipeptide has a lactone and triggers
            recalculation of the molecular weight
        """
        self._lac = value
        if self._core:
            self._calculate_mw()


def get_detected_domains(genes: Iterable[CDSFeature]) -> List[str]:
    """ Gathers all detected domains in a cluster, including some not detected
        by hmm_detection.

        Arguments:
            genes: a list of genes to check

        Returns:
            a list of strings, each string being the name of a domain in the
            cluster
    """
    found_domains: List[str] = []
    if not genes:
        return found_domains
    # Gather biosynthetic domains
    for feature in genes:
        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    # Gather non-biosynthetic domains
    cluster_fasta = get_fasta_from_features(genes)
    assert cluster_fasta
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(cluster_fasta)
    non_biosynthetic_hmms_found: List[str] = []
    for hsps_found in non_biosynthetic_hmms_by_id.values():
        for hsp in hsps_found:
            if hsp not in non_biosynthetic_hmms_found:
                non_biosynthetic_hmms_found.append(hsp)
    found_domains += non_biosynthetic_hmms_found

    return found_domains


def run_non_biosynthetic_phmms(fasta: str) -> Dict[str, List[str]]:
    """ Finds lanthipeptide-specific domains in the input fasta

        Arguments:
            fasta: a string containing gene sequences in fasta format

        Returns:
            a dictionary mapping the hit id to a list of matching query ids
    """
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"),
              "r", encoding="utf-8") as handle:
        hmmdetails = [line.strip().split("\t") for line in handle if line.count("\t") == 3]
    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    non_biosynthetic_hmms_by_id: Dict[str, List[str]] = defaultdict(list)
    for sig in signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path.rpartition(os.sep)[2])
        runresults = subprocessing.run_hmmsearch(sig.path, fasta)
        for runresult in runresults:
            # Store result if it is above cut-off
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    non_biosynthetic_hmms_by_id[hsp.hit_id].append(hsp.query_id)
    return non_biosynthetic_hmms_by_id


def predict_cleavage_site(query_hmmfile: str, target_sequence: str,
                          threshold: float = -100.) -> Optional[CleavageSiteHit]:
    """ Extracts from HMMER the start position, end position and score
        of the HMM alignment for a cleavage site

        Arguments:
            query_hmmfile: the path to a HMM file for the cleavage site profile
            target_sequence: the sequence of a CDS feature
            threshold: a minimum bitscore for a HMMer hit, exclusive

        Returns:
            a CleavageSiteHit instance with the information about the hit, or
            None if no hit was above the threshold
    """
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)

    for res in hmmer_res:
        for hits in res:
            lanthi_type = hits.description
            for hsp in hits:
                if hsp.bitscore > threshold:
                    return CleavageSiteHit(hsp.query_end, hsp.bitscore, lanthi_type)
    return None


def predict_class_from_genes(focus: CDSFeature, genes: Iterable[CDSFeature]) -> Optional[str]:
    """ Predict the lanthipeptide class from the gene cluster

        Arguments:
            genes: a list of genes to check

        Returns:
            a string representing the class, or None if no class predicted
    """
    found_domains = set()
    for feature in genes:
        if not feature.sec_met:
            continue
        found_domains.update(set(feature.sec_met.domain_ids))
    if focus.sec_met:
        found_domains.update(set(focus.sec_met.domain_ids))

    if 'Lant_dehydr_N' in found_domains or 'Lant_dehydr_C' in found_domains:
        return 'Class-I'
    if 'DUF4135' in found_domains:
        return 'Class-II'
    if 'Pkinase' in found_domains:
        # this could be class 3 or class 4, but as nobody has seen class 4
        # in vivo yet, we'll ignore that
        return 'Class-III'

    return None


def run_cleavage_site_phmm(fasta: str, hmmer_profile: str, threshold: float) -> Optional[CleavageSiteHit]:
    """ Try to identify cleavage site using pHMM """
    profile = path.get_full_path(__file__, hmmer_profile)
    return predict_cleavage_site(profile, fasta, threshold)


def run_cleavage_site_regex(fasta: str, lan_class: str) -> Optional[CleavageSiteHit]:
    """ Try to identify cleavage site using regular expressions"""
    # Regular expressions; try 1 first, then 2, etc.
    rex1 = re.compile('F?LD')
    rex2 = re.compile('[LF]?LQ')

    # For regular expression, check if there is a match that is <10 AA from the end
    if re.search(rex1, fasta) and len(re.split(rex1, fasta)[-1]) > 10:
        _, end = [m.span() for m in rex1.finditer(fasta)][-1]
        end += 16
    elif re.search(rex2, fasta) and len(re.split(rex2, fasta)[-1]) > 10:
        _, end = [m.span() for m in rex2.finditer(fasta)][-1]
        end += 15
    else:
        return None
    return CleavageSiteHit(end, -100., lan_class)


def determine_precursor_peptide_candidate(record: Record, query: CDSFeature, domains: List[str],
                                          hmmer_profile: str, lant_class: str) -> Optional[Lanthipeptide]:
    """ Identify precursor peptide candidates and split into two,
        only valid for Class-I lanthipeptides
    """

    if not MIN_PRECURSOR_LENGTH <= len(query.translation) <= MAX_PRECURSOR_LENGTH:
        return None

    # Create FASTA sequence for feature under study
    lan_a_fasta = f">query\n{query.translation}"

    # Run sequence against pHMM; if positive, parse into a vector containing START, END and SCORE
    cleavage_result = run_cleavage_site_phmm(lan_a_fasta, hmmer_profile, THRESH_DICT[lant_class])

    if cleavage_result is None or cleavage_result.end > len(query.translation) - 8:
        # If no pHMM hit, try regular expression
        cleavage_result = run_cleavage_site_regex(lan_a_fasta, lant_class)
        if cleavage_result is None or cleavage_result.end > len(query.translation) - 8:
            # still no good, so abort, since RODEO will predict duplicates based
            # only on cluster attributes
            return None

    # if the cleavage results in no core, that's not valid
    if cleavage_result.end == len(query.translation):
        return None

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(record, query, query.translation[:cleavage_result.end],
                             query.translation[cleavage_result.end:], domains)
    if rodeo_result < 14:
        return None

    # Determine the leader and core peptide
    leader = query.translation[:cleavage_result.end]
    core = query.translation[cleavage_result.end:]

    return Lanthipeptide(cleavage_result, rodeo_result, leader, core)


def run_lanthipred(record: Record, query: CDSFeature, lant_class: str,
                   domains: List[str]) -> Optional[Lanthipeptide]:
    """ Determines if a CDS is a predicted lanthipeptide based on the class
        and any contained domains.

        Arguments:
            record: the parent Record of the feature
            query: the CDSFeature to analyse
            lant_class: a string representing the class
            domains: a list of domain names in the current cluster
    """
    hmmer_profiles = {'Class-I': 'data/class1.hmm',
                      'Class-II': 'data/class2.hmm',
                      'Class-III': 'data/class3.hmm', }
    query_sequence = query.translation

    if lant_class in ("Class-II", "Class-III"):
        profile = path.get_full_path(__file__, hmmer_profiles[lant_class])
        lan_a_fasta = f">query\n{query_sequence}"
        cleavage_result = predict_cleavage_site(profile, lan_a_fasta)

        if cleavage_result is None:
            return None

        if THRESH_DICT[lant_class] > cleavage_result.score:
            return None

        # if the cleavage results in no core, that's not valid
        if cleavage_result.end == len(query_sequence):
            return None
        cleavage_result.lantype = lant_class
        leader = query_sequence[:cleavage_result.end]
        core = query_sequence[cleavage_result.end:]
        result = Lanthipeptide(cleavage_result, 0, leader, core)

    else:
        candidate = determine_precursor_peptide_candidate(record, query, domains,
                                                          hmmer_profiles[lant_class], lant_class)
        if candidate is None:
            return None
        result = candidate

    # extract now (that class is known and thus the END component) the core peptide
    if result.number_of_lan_bridges == 0:
        return None

    query.gene_functions.add(GeneFunction.ADDITIONAL, "lanthipeptides",
                             "predicted lanthipeptide")
    return result


def find_lan_a_features(area_feature: CDSCollection) -> List[CDSFeature]:
    """ Finds all lanthipeptide candidate features """
    lan_a_features = []
    for feature in area_feature.cds_children:
        assert feature.is_contained_by(area_feature)

        if len(feature.translation) < MAX_PRECURSOR_LENGTH:
            lan_a_features.append(feature)
            continue
        # if it has known precursor domains, that's also of interest if it's vaguely the right size
        if len(feature.translation) > MAX_PRECURSOR_LENGTH * 1.25:
            continue
        if feature.sec_met and set(feature.sec_met.domain_ids).intersection(KNOWN_PRECURSOR_DOMAINS):
            lan_a_features.append(feature)

    return lan_a_features


def contains_feature_with_domain(genes: List[CDSFeature], domains: Set[str]) -> bool:
    """ Checks for the existence of one of the given domains in the given features

        Arguments:
            genes: a list of genes to check
            domains: the set of domain names allowable

        Returns:
            True if a feature matching the conditions was found, otherwise False
    """
    for feature in genes:
        if not feature.sec_met:
            continue
        if domains.intersection(set(feature.sec_met.domain_ids)):
            return True
    return False


def result_vec_to_feature(orig_feature: CDSFeature, res_vec: Lanthipeptide) -> Prepeptide:
    """ Generates a Prepeptide feature from a CDSFeature and a Lanthipeptide

        Arguments:
            orig_feature: the CDSFeature the lanthipeptide was found in
            res_vec: the Lanthipeptide instance that was calculated

        Returns:
            a Prepeptide instance
    """
    assert res_vec.leader is not None
    product = "lanthipeptide"
    name = f"{orig_feature.get_name()}_{product}"
    feature = Prepeptide(orig_feature.location, product, res_vec.core,
                         name, "lanthipeptides", res_vec.lantype, res_vec.score,
                         res_vec.monoisotopic_mass, res_vec.molecular_weight,
                         res_vec.alternative_weights, res_vec.leader)
    qual = LanthiQualifier(res_vec.number_of_lan_bridges,
                           res_vec.rodeo_score, res_vec.aminovinyl_group,
                           res_vec.chlorinated, res_vec.oxygenated, res_vec.lactonated)
    feature.detailed_information = qual
    return feature


def find_neighbours_in_range(center: CDSFeature,
                             candidates: Iterable[CDSFeature]) -> List[CDSFeature]:
    """ Restrict a set of genes to those within precursor range of a central
        gene.

        Arguments:
            center: the gene to find the neighbours of
            candidates: the genes to filter by range

        Returns:
            a list of genes within range, with the same ordering as the input
    """
    neighbours = []
    for candidate in candidates:
        if candidate < center:
            if center.location.start - candidate.location.start <= MAX_PRECURSOR_DISTANCE:
                neighbours.append(candidate)
        else:
            if candidate.location.end - center.location.end <= MAX_PRECURSOR_DISTANCE:
                neighbours.append(candidate)
            else:
                # skip looking further to the right if the previous one was too far away
                break
    return neighbours


def run_lanthi_on_genes(record: Record, focus: CDSFeature, cluster: Protocluster,
                        genes: List[CDSFeature], results: LanthiResults) -> None:
    """ Runs lanthipeptide around a single focus gene which is a core biosynthetic
        enzyme for lanthipeptides.
        Updates the results object with any precursors found.

        Arguments:
            record: the Record instance containing the genes
            focus: a core lanthipeptide gene
            cluster: the Protocluster being analysed
            genes: a list of candidate precursor genes
            results: a LanthiResults object to update

        Returns:
            None
    """
    if not genes:
        return
    domains = get_detected_domains(cluster.cds_children)
    non_candidate_neighbours = find_neighbours_in_range(focus, cluster.cds_children)
    flavoprotein_found = contains_feature_with_domain(non_candidate_neighbours, {"Flavoprotein"})
    halogenase_found = contains_feature_with_domain(non_candidate_neighbours, {"Trp_halogenase"})
    oxygenase_found = contains_feature_with_domain(non_candidate_neighbours, {"p450"})
    dehydrogenase_found = contains_feature_with_domain(non_candidate_neighbours, {"adh_short", "adh_short_C2"})

    lant_class = predict_class_from_genes(focus, cluster.cds_children)
    if not lant_class:
        return

    for candidate in genes:
        result_vec = run_lanthipred(record, candidate, lant_class, domains)
        if result_vec is None:
            continue
        result_vec.aminovinyl_group = flavoprotein_found
        result_vec.chlorinated = halogenase_found
        result_vec.oxygenated = oxygenase_found
        result_vec.lactonated = dehydrogenase_found and result_vec.core.startswith('S')
        motif = result_vec_to_feature(candidate, result_vec)
        results.motifs_by_locus[focus.get_name()].append(motif)
        results.clusters[cluster.get_protocluster_number()].add(focus.get_name())
        # track new CDSFeatures if found with all_orfs
        if candidate.region is None:
            results.add_cds(candidate)


def run_specific_analysis(record: Record) -> LanthiResults:
    """ Runs the full lanthipeptide analysis over the given record

        Arguments:
            record: the Record instance to analyse

        Returns:
            A populated LanthiResults object
    """
    results = LanthiResults(record.id)
    for cluster in record.get_protoclusters():
        if not cluster.product.startswith("lanthipeptide"):
            continue

        # find core biosynthetic enzyme locations
        core_domain_names = {'Lant_dehydr_N', 'Lant_dehydr_C', 'DUF4135', 'Pkinase'}
        core_genes = []
        for gene in cluster.cds_children:
            if not gene.sec_met:
                continue
            # We seem to hit Lant_dehydr_C on some O-Methyltranferases that also hit PCMT
            if 'PCMT' in gene.sec_met.domain_ids:
                continue
            if core_domain_names.intersection(set(gene.sec_met.domain_ids)):
                core_genes.append(gene)

        precursor_candidates = find_lan_a_features(cluster)
        # Find candidate ORFs that are not yet annotated
        extra_orfs = all_orfs.find_all_orfs(record, cluster, min_length=MIN_PRECURSOR_LENGTH * 3)
        for orf in extra_orfs:
            if len(orf.translation) < MAX_PRECURSOR_LENGTH:
                precursor_candidates.append(orf)

        for gene in core_genes:
            neighbours = find_neighbours_in_range(gene, precursor_candidates)
            if not neighbours:
                continue
            run_lanthi_on_genes(record, gene, cluster, neighbours, results)

    logging.debug("Lanthipeptide module marked %d motifs", sum(map(len, results.motifs_by_locus)))

    cores = {}
    for precursors in results.motifs_by_locus.values():
        for precursor in precursors:
            cores[precursor.get_name().rsplit("_", 1)[0]] = precursor.core
    results.comparippson_results = comparippson.compare_precursor_cores(cores, get_config())
    return results
