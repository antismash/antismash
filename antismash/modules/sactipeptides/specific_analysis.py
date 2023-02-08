# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
More detailed sactipeptide analysis using HMMer-based leader peptide
cleavage sites prediction as well as prediction of number of  dissulfide
bridges, molcular mass and macrolactam ring.
'''

from collections import defaultdict
import logging
import os
import re
from typing import Any, Dict, Set, List, Optional, Tuple, Union

import joblib

from antismash.common import (
    all_orfs,
    comparippson,
    module_results,
    secmet,
    subprocessing,
    path,
    fasta,
    utils,
)
from antismash.common.hmmscan_refinement import HSP
from antismash.common.secmet.qualifiers import SecMetQualifier
from antismash.common.secmet.locations import location_from_string
from antismash.common.signature import HmmSignature
from antismash.config import get_config

# Cys chain limits (for CnnnC style sections)
CHAIN_LOWER = 1  # CnC
CHAIN_UPPER = 6  # CnnnnnnC

# the size limits, in aminos, of a precursor CDS
MIN_PRECURSOR_LENGTH = 20
MAX_PRECURSOR_LENGTH = 100


class SactiResults(module_results.ModuleResults):
    """ Holds the results of sactipeptide analysis for a record

    """
    schema_version = 3

    def __init__(self, record_id: str, comparippson_results: comparippson.MultiDBResults = None) -> None:
        super().__init__(record_id)
        # keep new CDS features
        self._new_cds_features: Set[secmet.CDSFeature] = set()
        # keep new CDSMotifs by the gene they match to
        # e.g. self.motifs_by_locus[gene_locus] = [motif1, motif2..]
        self.motifs_by_locus: Dict[str, List[secmet.Prepeptide]] = defaultdict(list)
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
        return {
            "record_id": self.record_id,
            "schema_version": SactiResults.schema_version,
            "motifs": motifs,
            "new_cds_features": cds_features,
            "protoclusters": {key: list(val) for key, val in self.clusters.items()},
            "comparippson_results": self.comparippson_results.to_json() if self.comparippson_results else None,
        }

    @staticmethod
    def from_json(json: Dict[str, Any], record: secmet.Record) -> Optional["SactiResults"]:
        if json.get("schema_version") != SactiResults.schema_version:
            logging.warning("Discarding Sactipeptide results, schema version mismatch")
            return None
        comp_results = comparippson.MultiDBResults.from_json(json["comparippson_results"])
        results = SactiResults(json["record_id"], comp_results)
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(secmet.Prepeptide.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["protoclusters"].items()}
        for location, name in json["new_cds_features"]:
            loc = location_from_string(location)
            cds = all_orfs.create_feature_from_location(record, loc, label=name)
            results.add_cds(cds)
        return results

    def add_to_record(self, record: secmet.Record) -> None:
        existing = record.get_cds_name_mapping()
        for feature in self._new_cds_features:
            # since a precursor may be found by other RiPP modules
            if feature.get_name() not in existing:
                record.add_cds_feature(feature)

        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                record.add_cds_motif(motif)

    def add_cds(self, cds: secmet.CDSFeature) -> None:
        """ Add a newly found CDS feature that will be added to the record """
        # if already added by another protocluster, don't double up
        for existing in self._new_cds_features:
            if cds.get_name() == existing.get_name():
                return
        self._new_cds_features.add(cds)


def get_detected_domains(cluster: secmet.Protocluster) -> Dict[str, int]:
    """ Gathers all detected domain ids from a cluster. Includes detection of
        some extra HMM profiles specific to sactipeptides.

        Arguments:
            cluster: the Protocluster to gather domains from

        Returns:
            a dictionary mapping domain ids to number of times that domain was found
    """
    found_domains: Dict[str, int] = {}
    # Gather biosynthetic domains
    for feature in cluster.cds_children:
        if not feature.sec_met:
            continue
        for domain_id in feature.sec_met.domain_ids:
            found_domains[domain_id] = found_domains.get(domain_id, 0) + 1

    # Gather non-biosynthetic domains
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(fasta.get_fasta_from_features(cluster.cds_children))
    for hsps_found_for_this_id in non_biosynthetic_hmms_by_id.values():
        for hsp in hsps_found_for_this_id:
            found_domains[hsp.query_id] = found_domains.get(hsp.query_id, 0) + 1

    return found_domains


def run_non_biosynthetic_phmms(cluster_fasta: str) -> Dict[str, List[HSP]]:
    """ Try to identify cleavage site using pHMM """
    if not cluster_fasta:
        return {}
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"),
              "r", encoding="utf-8") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]
    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]

    non_biosynthetic_hmms_by_id: Dict[str, Any] = defaultdict(list)
    for sig in signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path)
        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            # Store result if it is above cut-off
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    non_biosynthetic_hmms_by_id[hsp.hit_id].append(hsp)
    return non_biosynthetic_hmms_by_id


def cds_has_domains(cds: secmet.CDSFeature, domains: Set[str]) -> bool:
    """ Tests whether a cds has any of the given domains

        Arguments:
            cds: the CDSFeature to check
            domains: a set of domain names

        Returns:
            True if any of the domains are present in the CDS, otherwise False
    """
    return bool(cds.sec_met and set(cds.sec_met.domain_ids).intersection(domains))


def acquire_rodeo_heuristics(cluster: secmet.Protocluster, query: secmet.CDSFeature,
                             leader: str, core: str,
                             domains: Dict[str, int]) -> Tuple[int, List[float], List[int]]:
    """Calculate heuristic scores for RODEO"""
    tabs = []
    score = 0
    precursor = leader + core
    # Calcd. precursor peptide mass (Da)
    precursor_analysis = utils.RobustProteinAnalysis(precursor, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(precursor_analysis.molecular_weight()))
    # Calcd. leader peptide mass (Da)
    leader_analysis = utils.RobustProteinAnalysis(leader, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(leader_analysis.molecular_weight()))
    # Calcd. core peptide mass (Da)
    core_analysis = utils.RobustProteinAnalysis(core, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(core_analysis.molecular_weight()))
    # Distance to any biosynthetic protein (E, B, C)
    hmmer_profiles = ['PF04055']
    distance = utils.distance_to_pfam(cluster.parent_record, query, hmmer_profiles)
    tabs.append(distance)
    # rSAM within 500 nt?
    if utils.distance_to_pfam(cluster.parent_record, query, ['PF04055']) < 500:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # rSAM within 150 nt?
    if utils.distance_to_pfam(cluster.parent_record, query, ['PF04055']) < 150:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # rSAM further than 1000 nt?
    if utils.distance_to_pfam(cluster.parent_record, query, ['PF04055']) == -1 or \
       utils.distance_to_pfam(cluster.parent_record, query, ['PF04055']) > 10000:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Ratio of N-term to 1st Cys 0.25<x<0.60; Ratio of N-term to 1st Cys <0.25 or >0.60
    if "C" not in precursor:
        score -= 2
        tabs += [0, 1]
    elif 0.25 <= precursor.find("C") / len(precursor) <= 0.60:
        score += 2
        tabs += [1, 0]
    else:
        score -= 2
        tabs += [0, 1]
    # Three or more Cys; Less than 3 Cys
    if precursor.count("C") >= 3:
        score += 4
        tabs += [1, 0]
    else:
        score -= 4
        tabs += [0, 1]
    # CxC/CxxC/CxxxC/CxxxxxC; # CC/CCC
    motifs = (('C.{5}C', 2), ('C.{3}C', 1),
              ('C.{2}C', 1), ('C.{1}C', 1),
              ('CC', -2), ('CCC', -2))
    for motif in motifs:
        if re.search(motif[0], core):
            score += motif[1]
            tabs.append(1)
        else:
            tabs.append(0)
    # No Cys in last 1/4th?
    quarter_length = -len(precursor) // 4
    if "C" not in precursor[quarter_length:]:
        score += 1
        tabs.append(1)
    else:
        score -= 1
        tabs.append(0)
    # 2 Cys in first 2/3rds of precursor, 1 Cys in last 1/3rd of precursor
    two_thirds = 2 * len(precursor) // 3
    if precursor[:two_thirds].count("C") == 2 and precursor[two_thirds:].count("C") == 1:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Peptide matches SboA hmm
    if cds_has_domains(query, {"Subtilosin_A"}):
        score += 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Peptide matches SkfA hmm
    if cds_has_domains(query, {"TIGR04404"}):
        score += 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Peptide matches SCIFF hmm
    if cds_has_domains(query, {"TIGR03973"}):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has PqqD/RRE (PF05402)
    if "PF05402" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has SPASM domain (PF13186)
    if "SPASM" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # PF04055 (rSAM) domain start > 80
    runresults = subprocessing.run_hmmsearch(path.get_full_path(__file__, "data", "PF04055.hmm"),
                                             fasta.get_fasta_from_features(cluster.cds_children))
    max_start = 0
    hitstarts = []
    hitends = []
    for runresult in runresults:
        # Store result if it is above cut-off
        for hsp in runresult.hsps:
            if hsp.bitscore > 40:
                hitstarts.append(hsp.hit_start)
                max_start = max(hsp.hit_start, max_start)
                hitends.append(hsp.hit_end)
    if hitstarts and max_start > 80:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has peptidase
    peptidase_domains = ["Peptidase_M16_C", "Peptidase_S8", "Peptidase_M16", "Peptidase_S41"]
    no_peptidase = True
    for pepdom in peptidase_domains:
        if pepdom in domains:
            score += 1
            tabs.append(1)
            no_peptidase = False
        else:
            tabs.append(0)
    # cluster has transporter
    transport_domains = ["PF00005", "PF00664"]
    for transpdom in transport_domains:
        if transpdom in domains:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
    # cluster has response regulator (PF00072)
    if "PF00072" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has major facilitator (PF07690)
    if "PF07690" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has ATPase (PF13304)
    if "PF13304" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has Fer4_12 (PF13353)
    if "PF13353" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has rSAM (PF04055)
    if "PF04055" in domains or "TIGR03975" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # cluster has no recognized peptidase
    if no_peptidase:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # C-terminal portion is < 0.35 or > 0.65; C-terminal portion is defined as
    # the part from the last cysteine in the last identified Cx(n)C motif to the C-terminus
    # the binary opposite is also included as the next field
    last_motif_c = 0
    index = -1
    for aa in reversed(precursor):
        if aa == "C" and "C" in precursor[index-6:index]:
            last_motif_c = index + 1
        index -= 1
    if 0.35 <= last_motif_c / len(precursor) <= 0.65:
        score += 3
        tabs += [0, 1]
    else:
        score -= 2
        tabs += [1, 0]
    # SS profile count > 1
    # is there more than one Cx..C structure in the sequence
    cysrex = f"(?=(C.{{{CHAIN_LOWER},{CHAIN_UPPER}}}C))"
    rex4 = re.compile(cysrex)
    if len(rex4.findall(core)) > 1:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    return score, tabs, hitends


def structure_analysis(seq: str) -> Tuple[int, float, float, List[int]]:
    """ Gathers information on the structural properties of the sequence

        Arguments:
            seq: the sequence to analyse

        Returns:
            a tuple of
                the number of [G|A][G|A]C (and reversed) structures found
                the average distance between the GAC structures
                a normalised location of final C..C terminal within the sequence
                a list containing counts of each possible C...C size
    """
    # define Cys-gap structures with a c-terminal Cys
    cysrex = f"(?=(C.{{{CHAIN_LOWER},{CHAIN_UPPER}}}C))"

    core = str(seq)

    rex1 = re.compile(r'([G|A]{1,2}C)')
    rex2 = re.compile(r'(C[G|A]{1,2})')
    rex3 = re.compile(r'(C.{2,5})(?=C)')
    rex4 = re.compile(cysrex)

    locations = []
    matches = []
    for regex in [rex1, rex2]:
        for hit in regex.finditer(core):
            locations.append(hit.start())
            matches.append(hit.group())

    total = 0
    count = 0
    for i in range(len(locations) - 1):
        distance = locations[i + 1] - locations[i]
        if distance > 0:
            total += distance
        else:
            total += 1
        count += 1

    avg = 0.
    if total:
        avg = total / count

    cterm = len(rex3.split(core)[-1]) / len(core)

    profile = []
    chain_lengths = [len(ring) for ring in rex4.findall(core)]
    for size in range(CHAIN_LOWER, CHAIN_UPPER + 1):
        profile.append(chain_lengths.count(size + 2))  # includes the C on each end
    assert len(profile) == CHAIN_UPPER - CHAIN_LOWER + 1

    return len(matches), avg, cterm, profile


def generate_rodeo_svm_csv(leader: str, core: str, previously_gathered_tabs: List[Union[float, int]],
                           hitends: List[int], domain_counts: Dict[str, int]) -> List[float]:
    """Generates all the items for one candidate precursor peptide"""
    columns: List[float] = []
    precursor = leader + core
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns += previously_gathered_tabs
    # Length of leader peptide
    columns.append(len(leader))
    # Length of precursor peptide
    columns.append(len(precursor))
    # Length of core peptide
    columns.append(len(core))
    # Length of core / length of precursor ratio
    columns.append(len(core) / len(precursor))
    # Length of core / length of leader ratio
    columns.append(len(core) / len(leader))
    # Ratio of length of N-terminus to first Cys / length of precursor
    columns.append(precursor.count("C") / len(precursor))
    # Number of occurrences of CxNC motifs
    numrings, avg, cterm, profile = structure_analysis(precursor)
    columns.append(numrings)
    # Average distance between CxNC motifs
    columns.append(avg)
    # Ratio of length from last CxNC to C-terminus / length of core
    columns.append(cterm)
    # Number of instances of CxC where x exists (CHAIN_LOWER..CHAIN_UPPER) times
    columns.extend(profile)
    # Number in entire precursor of each amino acid
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in entire precursor of each amino acid type
    # (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
    amino_groups = ["FWY", "DE", "RK", "RKDE", "GAVLMI", "ST"]
    for group in amino_groups:
        columns.append(sum(precursor.count(aa) for aa in group))
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in leader of each amino acid type
    for group in amino_groups:
        columns.append(sum(leader.count(aa) for aa in group))
    # Number in core of each amino acid
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in core of each amino acid type
    for group in amino_groups:
        columns.append(sum(core.count(aa) for aa in group))
    # Number of each peptidase Pfam hit (PF05193/PF00082/PF03572/PF00675/PF02517/PF02163/PF00326)
    peptidase_domains = ["PF05193", "PF00082", "PF03572", "PF00675", "PF02517", "PF02163", "PF00326"]
    for peptidase in peptidase_domains:
        columns.append(domain_counts.get(peptidase, 0))
    # Number of each ABC transporter Pfam hit (PF00005/PF00664)
    transp_domains = ["PF00005", "PF00664"]
    for transp_dom in transp_domains:
        columns.append(domain_counts.get(transp_dom, 0))
    # Number of each response regulator Pfam hit (PF00072)
    columns.append(domain_counts.get("PF00072", 0))
    # Number of each major facilitator Pfam hit (PF07690)
    columns.append(domain_counts.get("PF07690", 0))
    # Number of each ATPase Pfam hit (PF02518/PF13304)
    atpase_domains = ["PF02518", "PF13304"]
    for atpase_dom in atpase_domains:
        columns.append(domain_counts.get(atpase_dom, 0))
    # Number of each Fer4_12 Pfam hit (PF13353)
    columns.append(domain_counts.get("PF13353", 0))
    # Number of each rSAM Pfam hit (PF04055)
    columns.append(domain_counts.get("PF04055", 0))
    # Length of rSAM including PqqD domain
    if not hitends:
        columns.append(0)
    else:
        columns.append(max(hitends))
    return columns


def run_rodeo_svm(csv_columns: List[float]) -> int:
    """Run RODEO SVM"""
    classifier_path = path.get_full_path(__file__, "data", "sactipeptide.classifier.pkl")
    scaler_path = path.get_full_path(__file__, "data", "sactipeptide.scaler.pkl")
    assert os.path.exists(classifier_path) and os.path.exists(scaler_path)
    classifier = joblib.load(classifier_path)
    scaler = joblib.load(scaler_path)
    csv_cols = [[float(i) for i in csv_columns[2:]]]
    scaled = scaler.transform(csv_cols)
    if int(classifier.predict(scaled)[0]) == 1:
        return 10
    return 0


def run_rodeo(cluster: secmet.Protocluster, query: secmet.CDSFeature,
              leader: str, core: str, domains: Dict[str, int]) -> Tuple[bool, float]:
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, partial_csv, hitends = acquire_rodeo_heuristics(cluster, query, leader, core, domains)
    rodeo_score += heuristic_score

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, partial_csv, hitends, domains)
    rodeo_score += run_rodeo_svm(csv_columns)

    return rodeo_score >= 26, rodeo_score


def determine_precursor_peptide_candidate(cluster: secmet.Protocluster, query: secmet.CDSFeature,
                                          query_sequence: str, domains: Dict[str, int]
                                          ) -> Optional[secmet.Prepeptide]:
    """Identify precursor peptide candidates and split into two"""

    if not MIN_PRECURSOR_LENGTH <= len(query_sequence) <= MAX_PRECURSOR_LENGTH:
        return None

    end = len(query_sequence) // 4  # TODO: this seems very arbitrary

    # Determine the leader and core peptide
    leader = query_sequence[:end]
    core = query_sequence[end:]

    # Run RODEO to assess whether candidate precursor peptide is judged real
    valid, score = run_rodeo(cluster, query, leader, core, domains)
    if not valid:
        return None
    product = "sactipeptide"
    name = f"{query.get_name()}_{product}"
    return secmet.Prepeptide(query.location, product, core, name,
                             tool="sactipeptides", leader=leader, score=score)


def run_sactipred(cluster: secmet.Protocluster, query: secmet.CDSFeature,
                  domains: Dict[str, int]) -> Optional[secmet.Prepeptide]:
    """General function to predict and analyse sacti peptides"""

    # Run checks to determine whether an ORF encodes a precursor peptide
    result = determine_precursor_peptide_candidate(cluster, query,
                                                   query.translation, domains)
    if result is None:
        return None

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "sactipeptides",
                             "predicted sactipeptide")
    return result


def annotate_orfs(cds_features: List[secmet.CDSFeature], hmm_results: Dict[str, List[HSP]]) -> None:
    """ Annotates newly found ORFs with sactipeptide domain information.
        This is only relevant for CDS features that did not exist during
        the cluster detection stage of antiSMASH.
    """

    domains_by_feature: Dict[str, List[SecMetQualifier.Domain]] = defaultdict(list)
    for hit_id, results in hmm_results.items():
        for result in results:
            domain = SecMetQualifier.Domain(result.query_id, result.evalue, result.bitscore, 0, "sactipeptides")
            domains_by_feature[hit_id].append(domain)
    for cds in cds_features:
        domains = domains_by_feature[cds.get_name()]
        if domains:
            cds.sec_met = SecMetQualifier(domains)


def specific_analysis(record: secmet.Record) -> SactiResults:
    """ Analyse each sactipeptide cluster and find precursors within it.
        If an unannotated ORF would contain the precursor, it will be annotated.

        Arguments:
            record: the Record to analyse

        Returns:
            a SactiResults instance holding all found precursors and new ORFs
    """
    results = SactiResults(record.id)
    new_feature_hits = 0
    motif_count = 0
    for cluster in record.get_protoclusters():
        if cluster.product != 'sactipeptide':
            continue

        # Find candidate ORFs that are not yet annotated
        new_orfs = all_orfs.find_all_orfs(record, cluster, min_length=MIN_PRECURSOR_LENGTH * 3)
        hmm_results = run_non_biosynthetic_phmms(fasta.get_fasta_from_features(new_orfs))
        annotate_orfs(new_orfs, hmm_results)

        # Get all CDS features to evaluate for RiPP-likeness
        candidates = list(cluster.cds_children) + new_orfs
        domains = get_detected_domains(cluster)

        # Evaluate each candidate precursor peptide
        for candidate in candidates:
            motif = run_sactipred(cluster, candidate, domains)
            if motif is None:
                continue

            results.motifs_by_locus[candidate.get_name()].append(motif)
            motif_count += 1
            results.clusters[cluster.get_protocluster_number()].add(candidate.get_name())
            # track new CDSFeatures if found with all_orfs
            if candidate.region is None:
                results.add_cds(candidate)
                new_feature_hits += 1

    if not motif_count:
        logging.debug("Found no sactipeptide motifs")
    else:
        verb = "is" if new_feature_hits == 1 else "are"
        logging.debug("Found %d sactipeptide motif(s) in %d feature(s), %d of which %s new",
                      motif_count, len(results.motifs_by_locus), new_feature_hits, verb)

    cores = {}
    for motifs in results.motifs_by_locus.values():
        for motif in motifs:
            cores[motif.get_name().rsplit("_", 1)[0]] = motif.core
    results.comparippson_results = comparippson.compare_precursor_cores(cores, get_config())
    return results
