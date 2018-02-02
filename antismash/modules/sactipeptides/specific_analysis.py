# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
More detailed sactipeptide analysis using HMMer-based leader peptide
cleavage sites prediction as well as prediction of number of  dissulfide
bridges, molcular mass and macrolactam ring.
'''

from collections import defaultdict
import logging
import numpy as np
import os
import re
from typing import Any, Dict, Set, List, Optional

from sklearn.externals import joblib

from antismash.common import utils, all_orfs, module_results, secmet, serialiser, subprocessing, path, fasta
from antismash.common.signature import HmmSignature


class SactiResults(module_results.ModuleResults):
    """ Holds the results of sactipeptide analysis for a record

    """
    schema_version = 1

    def __init__(self, record_id, *args):
        super().__init__(record_id, *args)
        # keep new CDS features
        self.new_cds_features = set()
        # keep new CDSMotifs by the gene they match to
        # e.g. self.motifs_by_locus[gene_locus] = [motif1, motif2..]
        self.motifs_by_locus = defaultdict(list)
        # keep clusters and which genes in them had precursor hits
        # e.g. self.clusters[cluster_number] = {gene1_locus, gene2_locus}
        self.clusters = defaultdict(set)

    def to_json(self):
        cds_features = [(serialiser.location_to_json(feature.location),
                         feature.get_name()) for feature in self.new_cds_features]
        motifs = {}
        for locus, locus_motifs in self.motifs_by_locus.items():
            motifs[locus] = [motif.to_json() for motif in locus_motifs]
        return {"record_id": self.record_id,
                "schema_version": SactiResults.schema_version,
                "motifs": motifs,
                "new_cds_features": cds_features,
                "clusters": {key: list(val) for key, val in self.clusters.items()}}

    @staticmethod
    def from_json(json, record) -> "SactiResults":
        if json.get("schema_version") != SactiResults.schema_version:
            logging.warning("Discarding Sactipeptide results, schema version mismatch")
            return None
        results = SactiResults(json["record_id"])
        for locus, motifs in json["motifs"].items():
            for motif in motifs:
                results.motifs_by_locus[locus].append(SactipeptideMotif.from_json(motif))
        results.clusters = {int(key): set(val) for key, val in json["clusters"].items()}
        for location, name in json["new_cds_features"]:
            cds = all_orfs.create_feature_from_location(record, location, label=name)
            results.new_cds_features.add(cds)
        return results

    def add_to_record(self, record):
        for feature in self.new_cds_features:
            record.add_cds_feature(feature)

        for motifs in self.motifs_by_locus.values():
            for motif in motifs:
                record.add_cds_motif(motif)


class Sactipeptide:
    '''
    Class to calculate and store sactipeptide information
    '''
    def __init__(self, start, end, rodeo_score):
        self.start = start
        self.end = end
        self.rodeo_score = rodeo_score
        self._leader = ''
        self._core = ''

    @property
    def core(self):
        return self._core

    @core.setter
    def core(self, seq):
        seq = seq.replace('X', '')
        self._core = seq

    @property
    def leader(self):
        return self._leader

    @leader.setter
    def leader(self, seq):
        self._leader = seq

    @property
    def c_cut(self):
        return self._c_cut

    @c_cut.setter
    def c_cut(self, ccut):
        self._c_cut = ccut


def get_detected_domains(cluster: secmet.Cluster) -> Set[str]:
    """ Gathers all detected domain ids from a cluster. Includes detection of
        some extra HMM profiles specific to sactipeptides.

        Arguments:
            cluster: the Cluster to gather domains from

        Returns:
            a dictionary mapping domain ids to number of times that domain was found
    """
    found_domains = {}  # type: Set[str]
    # Gather biosynthetic domains
    for feature in cluster.cds_children:
        if not feature.sec_met:
            continue
        for domain_id in feature.sec_met.domain_ids:  # TODO: check domain_ids is not a set/unique only
            found_domains[domain_id] = found_domains.get(domain_id, 0) + 1

    # Gather non-biosynthetic domains
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(fasta.get_fasta_from_features(cluster.cds_children))
    for hsps_found_for_this_id in non_biosynthetic_hmms_by_id.values():
        for hsp in hsps_found_for_this_id:
            found_domains[hsp.query_id] = found_domains.get(hsp.query_id, 0) + 1

    return found_domains


def run_non_biosynthetic_phmms(cluster_fasta: str) -> Dict[str, Any]:
    """ Try to identify cleavage site using pHMM """
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"), "r") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]
    signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]

    non_biosynthetic_hmms_by_id = defaultdict(list)  # type: Dict[str, Any]
    for sig in signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path.rpartition(os.sep)[2])
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
    return cds.sec_met and set(cds.sec_met.domain_ids).intersection(domains)


def acquire_rodeo_heuristics(record, cluster, query, leader, core, domains):
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
    distance = utils.distance_to_pfam(record, query, hmmer_profiles)
    tabs.append(distance)
    # rSAM within 500 nt?
    if utils.distance_to_pfam(record, query, ['PF04055']) < 500:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # rSAM within 150 nt?
    if utils.distance_to_pfam(record, query, ['PF04055']) < 150:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # rSAM further than 1000 nt?
    if utils.distance_to_pfam(record, query, ['PF04055']) == -1 or \
       utils.distance_to_pfam(record, query, ['PF04055']) > 10000:
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
    motifs = (('C[ARNDBCEQZGHILKMFPSTWYV]{5}C', 2), ('C[ARNDBCEQZGHILKMFPSTWYV]{3}C', 1),
              ('C[ARNDBCEQZGHILKMFPSTWYV]{2}C', 1), ('C[ARNDBCEQZGHILKMFPSTWYV]{1}C', 2),
              ('CC', -2), ('CCC', -2))
    for motif in motifs:
        if re.search(motif[0], core):
            score += motif[1]
            tabs.append(1)
        else:
            tabs.append(0)
    # No Cys in last 1/4th?
    quarter_length = int(-len(precursor) / 4)
    if "C" not in precursor[quarter_length:]:
        score += 1
        tabs.append(1)
    else:
        score -= 1
        tabs.append(0)
    # 2 Cys in first 2/3rds of precursor, 1 Cys in last 1/3rd of precursor
    two_thirds = int(2 * len(precursor) / 3)
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
    # Cluster has PqqD/RRE (PF05402)
    if "PF05402" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has SPASM domain (PF13186)
    if "PF13186" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # PF04055 (rSAM) domain start
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
    # Cluster has peptidase
    peptidase_domains = ["Peptidase_M16_C", "Peptidase_S8", "Peptidase_M16", "Peptidase_S41"]
    no_peptidase = True
    for pepdom in peptidase_domains:
        if pepdom in domains:
            score += 1
            tabs.append(1)
            no_peptidase = False
        else:
            tabs.append(0)
    # Cluster has transporter
    transport_domains = ["PF00005", "PF00664"]
    for transpdom in transport_domains:
        if transpdom in domains:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
    # Cluster has response regulator (PF00072)
    if "PF00072" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has major facilitator (PF07690)
    if "PF07690" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has ATPase (PF13304)
    if "PF13304" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has Fer4_12 (PF13353)
    if "PF13353" in domains:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has rSAM (PF04055)
    if "PF04055" in domains or "TIGR03975" in domains:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster has no recognized peptidase
    if no_peptidase:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # C-terminal portion is < 0.35 or > 0.65; C-terminal portion is defined as
    # the part from the last cysteine in the last identified Cx(n)C motif to the C-terminus
    # And criterion "C-terminal portion is > 0.35 and < 0.65"
    if "C" not in precursor:
        score -= 2
        tabs += [1, 0]
    else:
        last_motif_C = 0
        index = -1
        for aa in reversed(precursor):
            if aa == "C" and "C" in precursor[index-6:index]:
                last_motif_C = len(precursor[:index]) + 1
            index -= 1
        if not (0.35 <= last_motif_C / float(len(precursor)) <= 0.65):
            score -= 2
            tabs += [1, 0]
        else:
            score += 3
            tabs += [1, 0]
    # SS profile sum > 1; The "SS profile" variables are calculated by the
    # lanscout function, namely the last regex search done using the rex4.
    # This provides a count for each number in the range
    lanlower = 1
    lanupper = 6
    cysrex = '(?=(C.{%d,%d}C))' % (lanlower, lanupper)
    rex4 = re.compile(cysrex)
    totrings = rex4.findall(core)
    if len(totrings) > 1:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    return score, tabs, hitends


def lanscout(seq):
    lanlower = 1
    lanupper = 6
    # define lanSactinine ring with a c-terminal Cys
    cysrex = '(?=(C.{%d,%d}C))' % (lanlower, lanupper)

    db = []
    avgs = []
    numrings = []
    cterm = []
    for n in range(0, len(seq)):
        core = str(seq[n][:])

        rex1 = re.compile(r'([G|A]{1,2}C)')
        rex2 = re.compile(r'(C[G|A]{1,2})')
        rex3 = re.compile(r'(C.{2,5})(?=C)')
        rex4 = re.compile(cysrex)

        loc1 = []
        match1 = []
        for county in rex1.finditer(core):
            loc1.append(county.start())
            match1.append(county.group())

        for county in rex2.finditer(core):
            loc1.append(county.start())
            match1.append(county.group())

        tempy = []
        for m in range(len(loc1[:-1])):
            if (loc1[m+1]-loc1[m]) > 0:
                tempy.append(loc1[m+1]-loc1[m])
            else:
                tempy.append(1)

        if tempy:
            avgs.append(np.mean(tempy))
        else:
            avgs.append(0)

        numrings.append(len(match1))
        cterm.append(len(rex3.split(core)[-1])/float(len(core)))

        numringlist = []
        totrings = rex4.findall(core)
        size = []
        for i in range(0, len(totrings)):
            size.append(len(totrings[i]))
        db.append(size)
        numringlist.append(numrings)

    profile = []
    for i in range(len(db)):
        temp = []
        for j in range(lanlower + 2, lanupper + 3):
            temp.append(db[i].count(j))
        profile.append(temp)

    for i in range(len(profile)):
        profile[i] = str(profile[i]).strip('[]')

    return numrings, avgs, cterm, profile


def generate_rodeo_svm_csv(leader, core, previously_gathered_tabs, hitends, domain_counts) -> List[float]:
    """Generates all the items for one candidate precursor peptide"""
    columns = []  # type: List[float]
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
    columns.append(float(len(core)) / float(len(precursor)))
    # Length of core / length of leader ratio
    columns.append(float(len(core)) / float(len(leader)))
    # Ratio of length of N-terminus to first Cys / length of precursor
    columns.append(precursor.count("C") / float(len(precursor)))
    # Number of occurrences of CxNC motifs
    numrings, avgs, cterm, profile = lanscout([precursor])
    if np.isnan(avgs[0]):
        avgs = [0]
    columns.append(numrings[0])
    # Average distance between CxNC motifs
    columns.append(avgs[0])
    # Ratio of length from last CxNC to C-terminus / length of core
    columns.append(cterm[0])
    # Number of instances of CxNC where N = 1
    columns.append(profile[0].split(",")[0])
    # Number of instances of CxNC where N = 2
    columns.append(profile[0].split(",")[1])
    # Number of instances of CxNC where N = 3
    columns.append(profile[0].split(",")[2])
    # Number of instances of CxNC where N = 4
    columns.append(profile[0].split(",")[3])
    # Number of instances of CxNC where N = 5
    columns.append(profile[0].split(",")[4])
    # Number of instances of CxNC where N = 6
    columns.append(profile[0].split(",")[5])
    # Number in entire precursor of each amino acid
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in entire precursor of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
    columns.append(sum([precursor.count(aa) for aa in "FWY"]))
    columns.append(sum([precursor.count(aa) for aa in "DE"]))
    columns.append(sum([precursor.count(aa) for aa in "RK"]))
    columns.append(sum([precursor.count(aa) for aa in "RKDE"]))
    columns.append(sum([precursor.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([precursor.count(aa) for aa in "ST"]))
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in leader of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
    columns.append(sum([leader.count(aa) for aa in "FWY"]))
    columns.append(sum([leader.count(aa) for aa in "DE"]))
    columns.append(sum([leader.count(aa) for aa in "RK"]))
    columns.append(sum([leader.count(aa) for aa in "RKDE"]))
    columns.append(sum([leader.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([leader.count(aa) for aa in "ST"]))
    # Number in core of each amino acid
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in core of each amino acid type (Aromatics, Neg charged, Pos charged, Charged, Aliphatic, Hydroxyl)
    columns.append(sum([core.count(aa) for aa in "FWY"]))
    columns.append(sum([core.count(aa) for aa in "DE"]))
    columns.append(sum([core.count(aa) for aa in "RK"]))
    columns.append(sum([core.count(aa) for aa in "RKDE"]))
    columns.append(sum([core.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([core.count(aa) for aa in "ST"]))
    # Number of each peptidase Pfam hit (PF05193/PF00082/PF03572/PF00675/PF02517/PF02163/PF00326)
    peptidase_domains = ["PF05193", "PF00082", "PF03572", "PF00675", "PF02517", "PF02163", "PF00326"]
    for pepdom in peptidase_domains:
        if pepdom in domain_counts:
            columns.append(domain_counts[pepdom])
        else:
            columns.append(0)
    # Number of each ABC transporter Pfam hit (PF00005/PF00664)
    transp_domains = ["PF00005", "PF00664"]
    for transp_dom in transp_domains:
        if transp_dom in domain_counts:
            columns.append(domain_counts[transp_dom])
        else:
            columns.append(0)
    # Number of each response regulator Pfam hit (PF00072)
    if "PF00072" in domain_counts:
        columns.append(domain_counts["PF00072"])
    else:
        columns.append(0)
    # Number of each major facilitator Pfam hit (PF07690)
    if "PF07690" in domain_counts:
        columns.append(domain_counts["PF07690"])
    else:
        columns.append(0)
    # Number of each ATPase Pfam hit (PF02518/PF13304)
    atpase_domains = ["PF02518", "PF13304"]
    for atpase_dom in atpase_domains:
        if atpase_dom in domain_counts:
            columns.append(domain_counts[atpase_dom])
        else:
            columns.append(0)
    # Number of each Fer4_12 Pfam hit (PF13353)
    if "PF13353" in domain_counts:
        columns.append(domain_counts["PF13353"])
    else:
        columns.append(0)
    # Number of each rSAM Pfam hit (PF04055)
    if "PF04055" in domain_counts:
        columns.append(domain_counts["PF04055"])
    else:
        columns.append(0)
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


def run_rodeo(record, cluster, query, leader, core, domains):
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv, hitends = acquire_rodeo_heuristics(record, cluster, query, leader, core, domains)
    rodeo_score += heuristic_score

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, gathered_tabs_for_csv, hitends, domains)
    rodeo_score += run_rodeo_svm(csv_columns)

    if rodeo_score >= 26:
        return True, rodeo_score
    else:
        return False, rodeo_score


def determine_precursor_peptide_candidate(record, cluster, query, query_sequence, domains) -> Optional[Sactipeptide]:
    """Identify precursor peptide candidates and split into two"""

    # Skip sequences with >100 AA
    if not 20 <= len(query_sequence) <= 100:
        return None

    end = int(len(query_sequence)*0.25)

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(record, cluster, query, query_sequence[:end], query_sequence[end:], domains)
    if not rodeo_result[0]:
        return None

    sacti_peptide = Sactipeptide(0, end + 1, rodeo_result[1])

    # Determine the leader and core peptide
    sacti_peptide.leader = query_sequence[:end]
    sacti_peptide.core = query_sequence[end:]

    return sacti_peptide


def run_sactipred(record, cluster, query, domains) -> Optional[Sactipeptide]:
    """General function to predict and analyse sacti peptides"""

    # Run checks to determine whether an ORF encodes a precursor peptide
    result = determine_precursor_peptide_candidate(record, cluster, query,
                                                   query.translation, domains)
    if result is None:
        return None

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "sactipeptides",
                             "predicted sactipeptide")
    return result


class SactipeptideMotif(secmet.Prepeptide):
    def __init__(self, location, name, score, leader, core):
        super().__init__(location, "sactipeptide", core, name, leader=leader, score=score)

    def to_json(self) -> Dict[str, Any]:
        return {"location": serialiser.location_to_json(self.location),
                "name": self.get_name(),
                "score": self.score,
                "core": self.core,
                "leader": self.leader}

    @staticmethod
    def from_json(json) -> "SactipeptideMotif":
        return SactipeptideMotif(serialiser.location_from_json(json["location"]),
                                 json["name"], json["score"], json["leader"], json["core"])


def result_vec_to_features(orig_feature, res_vec):
    return SactipeptideMotif(orig_feature.location, orig_feature.get_name(),
                             res_vec.rodeo_score, res_vec.leader, res_vec.core)


def specific_analysis(record, options) -> SactiResults:
    results = SactiResults(record.id)
    for cluster in record.get_clusters():
        if 'sactipeptide' not in cluster.products:
            continue

        # Find candidate ORFs that are not yet annotated
        new_orfs = all_orfs.find_all_orfs(record, cluster)

        # Get all CDS features to evaluate for RiPP-likeness
        candidates = list(cluster.cds_children) + new_orfs
        domains = get_detected_domains(cluster)

        # Evaluate each candidate precursor peptide
        for candidate in candidates:
            result_vec = run_sactipred(record, cluster, candidate, domains)
            if result_vec is None:
                continue

            motif = result_vec_to_features(candidate, result_vec)

            results.motifs_by_locus[candidate.get_name()].append(motif)
            results.clusters[cluster.get_cluster_number()].add(candidate.get_name())
            # track new CDSFeatures if found with all_orfs
            if candidate.cluster is None:
                results.new_cds_features.add(candidate)
    return results
