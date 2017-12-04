# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
More detailed thiopeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of molcular mass.
'''

from collections import defaultdict
import logging
import re
import os
from typing import List

from Bio.SeqFeature import FeatureLocation
import numpy as np
from sklearn.externals import joblib

from antismash.common import all_orfs, fasta, module_results, path, secmet, serialiser, subprocessing, utils
from antismash.detection.hmm_detection.signatures import HmmSignature


class ThioResults(module_results.ModuleResults):
    schema_version = 1

    def __init__(self, record_id, *args):
        super().__init__(record_id, *args)
        self.clusters_with_motifs = set()
        self.cds_features = defaultdict(list)  # CDSs found with find_all_orfs
        self.motifs = []

    def to_json(self):
        cds_features_by_cluster = {key: [(serialiser.location_to_json(feature.location), feature.get_name()) for feature in features]
                                   for key, features in self.cds_features.items()}
        return {"record_id": self.record_id,
                "schema_version": ThioResults.schema_version,
                "clusters with motifs": [cluster.get_cluster_number() for cluster in self.clusters_with_motifs],
                "motifs": [motif.to_json() for motif in self.motifs],
                "cds_features": cds_features_by_cluster}

    @staticmethod
    def from_json(data, record):
        if data.get("schema_version") != ThioResults.schema_version:
            logging.warning("Discarding Thiopeptide results, schema version mismatch")
            return None
        results = ThioResults(data["record_id"])
        for motif in data["motifs"]:
            results.motifs.append(ThiopeptideMotif.from_json(motif))
        for cluster in data["clusters with motifs"]:
            results.clusters_with_motifs.add(record.get_cluster(cluster))
        for cluster, features in data["cds_features"]:
            for location, name in features:
                cds = all_orfs.create_feature_from_location(record, location, label=name)
                results.cds_features[cluster].append(cds)
        return results

    def add_to_record(self, record):
        for features in self.cds_features.values():
            for cds in features:
                record.add_cds_feature(cds)

        for motif in self.motifs:
            record.add_cds_motif(motif)


class Thiopeptide(object):
    '''
    Class to calculate and store thiopeptide information
    '''
    def __init__(self, start, end, score, rodeo_score):
        self.start = start
        self.end = end
        self.score = score
        self.rodeo_score = rodeo_score
        self.thio_type = ''
        self._leader = ''
        self._core = ''
        self._weight = -1
        self._monoisotopic_weight = -1
        self._alt_weights = []
        self._macrocycle = ''
        self._amidation = False
        self._mature_alt_weights = []
        self._mature_features = ''
        self._c_cut = ''
        self.core_analysis_monoisotopic = None
        self.core_analysis = None

    @property
    def core(self):
        return self._core

    @core.setter
    def core(self, seq):
        self.core_analysis_monoisotopic = utils.RobustProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = utils.RobustProteinAnalysis(seq, monoisotopic=False)
        self._core = seq

    @property
    def leader(self):
        return self._leader

    @leader.setter
    def leader(self, seq):
        self._leader = seq

    @property
    def macrocycle(self):
        self._predict_macrocycle()
        return self._macrocycle

    @macrocycle.setter
    def macrocycle(self, macro):
        self._macrocycle = macro

    def _predict_macrocycle(self):
        '''
        Prediction of macrocycle ring
        '''

        if not self._core:
            raise ValueError()

        aux = self._core
        if len(self._core) == 17 and self._core[0] in "TIV":
            aux = self._core[4:]

        if aux[0] == "S":
            if "SC" in aux[9:12]:
                if aux[9] == "S":
                    self._macrocycle = "26-member"

                if aux[10] == "S":
                    self._macrocycle = "29-member"

            if re.search(r"S{6,8}?", aux[len(aux)//2:]):
                if aux[12] == "S":
                    self._macrocycle = "35-member"

    @property
    def amidation(self):
        return self._amidation

    @amidation.setter
    def amidation(self, value):
        self._amidation = value

    @property
    def mature_features(self):
        self._predict_mature_core_features()
        return self._mature_features

    @mature_features.setter
    def mature_features(self, feats):
        self._mature_features = feats

    def _predict_mature_core_features(self):
        '''
        Prediction of some features of the mature peptide
        '''
        assert self.thio_type

        self._mature_features = ''
        if self.thio_type == 'Type-I':
            self._mature_features = "Central ring: pyridine tetrasubstituted (hydroxyl group present); second macrocycle"

        if self.thio_type == 'Type-II':
            self._mature_features = "Central ring: piperidine; second macrocycle containing a quinaldic acid moiety"

        if self.thio_type == 'Type-III':
            self._mature_features = "Central ring: pyridine trisubstituted"

        return self._mature_features

    @property
    def c_cut(self):
        return self._c_cut

    @c_cut.setter
    def c_cut(self, ccut):
        self._c_cut = ccut

    def __repr__(self):
        return "Thiopeptide(%s..%s, %s, %r, %r, %s(%s), %s, %s, %s, %s)" % (
                        self.start, self.end, self.score, self._core,
                        self.thio_type, self._monoisotopic_weight, self._weight,
                        self.macrocycle, self.amidation, self._mature_features,
                        self.c_cut)

    def _calculate_mw(self):
        '''
        (re)calculate the monoisotopic mass and molecular weight
        '''
        assert self._core, "calculating weight without a core"

        amino_counts = self.core_analysis.count_amino_acids()
        no_thr_ser = amino_counts['T'] + amino_counts['S']
        no_cys = amino_counts['C']

        mol_mass = self.core_analysis.molecular_weight()
        monoisotopic_mass = self.core_analysis_monoisotopic.molecular_weight()
        self._weight = mol_mass
        self._monoisotopic_weight = monoisotopic_mass

        water_weight = 18.02

        dehydrations = water_weight * (no_thr_ser-1)

        # every unbridged Ser or Thr might not be dehydrated
        for i in range(1, no_thr_ser + 1):
            self._alt_weights.append(mol_mass + water_weight * i)

        if self.thio_type != "Type-III":
            # maturation reactions:
            mol_mass -= dehydrations

            # YcaO > cys/ser to azole ~20 Da for ring
            mol_mass -= (20 * (no_cys - 1))

            # cycloaddition of 2 Dha residues
            mol_mass -= 35

            if self.thio_type == 'Type-I':
                # indolic acid (172) + cyclization
                mol_mass += 172 - 3
            elif self.thio_type == 'Type-II':
                # quinaldic acid + cyclization
                mol_mass += 215

            # considering 2 oxidations
            mol_mass += 32

        if self.amidation:
            mol_mass -= 70

        mass_diff = self._weight - mol_mass
        mod_mol_mass = mol_mass
        mod_monoisotopic_mass = monoisotopic_mass - mass_diff

        self._mature_alt_weights.append(mod_mol_mass)
        self._mature_alt_weights.append(mod_monoisotopic_mass)

        # every unbridged Ser or Thr might not be dehydrated
        for i in range(1, no_thr_ser + 1):
            self._mature_alt_weights.append(mod_mol_mass + water_weight * i)

    @property
    def monoisotopic_mass(self):
        '''
        function determines the weight of the core peptide and substracts
        the weight which is reduced, due to dehydratation
        '''
        if self._monoisotopic_weight > -1:
            return self._monoisotopic_weight

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self):
        '''
        function determines the weight of the core peptide
        '''
        if self._weight > -1:
            return self._weight

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self):
        '''
        function determines the possible alternative weights assuming one or
        more of the Ser/Thr residues aren't dehydrated
        '''
        if self._alt_weights != []:
            return self._alt_weights

        self._calculate_mw()
        return self._alt_weights

    @property
    def mature_alt_weights(self):
        '''
        function determines the possible mature alternative weights which includes mw,
        monoisotopic mass and alternative weightss due to different hydrated residues

        '''
        if self._mature_alt_weights != []:
            return self._mature_alt_weights

        self._calculate_mw()
        return self._mature_alt_weights


def predict_amidation(found_domains):
    amidation = False
    # nosA homologs nocA, tpdK nocA berI pbt
    if 'thio_amide' in found_domains:
        amidation = True

    return amidation


def predict_cleavage_site(query_hmmfile, target_sequence, threshold):
    '''
    Function extracts from HMMER the start position, end position and score
    of the HMM alignment
    '''
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)
    resvec = [None, None, None]

    for res in hmmer_res:
        for hits in res:
            for hsp in hits:
                if hsp.bitscore > threshold:
                    resvec = [hsp.query_start, hsp.query_end-14, hsp.bitscore]
                    return resvec

    return resvec


def predict_type_from_gene_cluster(found_domains):
    '''
    Predict the thiopeptide type from the gene cluster
    '''

    if 'PF06968' in found_domains:
        return 'Type-I'

    if 'PF04055' in found_domains or 'PF00733' in found_domains:
        return 'Type-II'

    return 'Type-III'


def get_detected_domains(cluster):

    found_domains = []
    # Gather biosynthetic domains
    for feature in cluster.cds_children:

        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    # Gather non-biosynthetic domains
    cluster_features = cluster.cds_children
    cluster_fasta = fasta.get_fasta_from_features(cluster_features)
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(cluster_fasta)
    non_biosynthetic_hmms_found = []
    for hmm_id, hsps_found_for_this_id in non_biosynthetic_hmms_by_id.items():
        hsps_found_for_this_id = non_biosynthetic_hmms_by_id[hmm_id]
        for hsp in hsps_found_for_this_id:
            if hsp.query_id not in non_biosynthetic_hmms_found:
                non_biosynthetic_hmms_found.append(hsp.query_id)
    found_domains.extend(non_biosynthetic_hmms_found)

    return found_domains


def run_non_biosynthetic_phmms(cluster_fasta):
    """Try to identify cleavage site using pHMM"""
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"), "r") as handle:
        hmmdetails = [line.split("\t") for line in handle.read().splitlines() if line.count("\t") == 3]
    _signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    non_biosynthetic_hmms_by_id = {}
    for sig in _signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path.rpartition(os.sep)[2])
        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            # Store result if it is above cut-off
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    if hsp.hit_id not in non_biosynthetic_hmms_by_id:
                        non_biosynthetic_hmms_by_id[hsp.hit_id] = [hsp]
                    else:
                        non_biosynthetic_hmms_by_id[hsp.hit_id].append(hsp)
    return non_biosynthetic_hmms_by_id


def run_cleavage_site_phmm(input_fasta, hmmer_profile, threshold):
    """Try to identify cleavage site using pHMM"""
    profile = path.get_full_path(__file__, "data", hmmer_profile)
    return predict_cleavage_site(profile, input_fasta, threshold)


def run_cleavage_site_regex(sequence):
    """Try to identify cleavage site using regular expressions"""
    # Regular expressions; try 1 first, then 2, etc.
    rex1 = re.compile('([I|V]AS)')
    rex2 = re.compile('([G|A|S]AS)')

    # For each regular expression, check if there is a match that is <10 AA from the end
    if re.search(rex1, sequence) and len(re.split(rex1, sequence)[-1]) > 10:
        start, end = [m.span() for m in rex1.finditer(sequence)][-1]
        end -= 5
    elif re.search(rex2, sequence) and len(re.split(rex1, sequence)[-1]) > 10:
        start, end = [m.span() for m in rex2.finditer(sequence)][-1]
        end -= 5
    else:
        return None, None, None

    return start, end, 0


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
        rexout2 = np.mean([len(x) for x in rexout2])
        avg_heteroblock_length = float(str(rexout2).strip("[]"))
    else:
        avg_heteroblock_length = 0.0

    return number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks


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
    if avg_heteroblock_length != "nan" and avg_heteroblock_length > 3:
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
    if avg_heteroblock_length == "nan":
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


def determine_precursor_peptide_candidate(query, query_sequence, domains):
    """Identify precursor peptide candidates and split into two"""

    # Skip sequences with >200 AA
    if len(query_sequence) > 200 or len(query_sequence) < 40:
        return

    # Create FASTA sequence for feature under study
    thio_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    # Run sequence against pHMM; if positive, parse into a vector containing START, END and SCORE
    start, end, score = run_cleavage_site_phmm(thio_a_fasta, 'thio_cleave.hmm', -3.00)

    # If no pHMM hit, try regular expression
    if score is None:
        start, end, score = run_cleavage_site_regex(query_sequence)
        if score is None or end > len(query_sequence) - 5:
            start, end, score = 0, int(len(query_sequence)*0.60) - 14, 0

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(query_sequence[:end], query_sequence[end:], domains)
    if rodeo_result[0] is False:
        return

    thiopeptide = Thiopeptide(start, end + 1, score, rodeo_result[1])

    # Determine the leader and core peptide
    thiopeptide.leader = query_sequence[:end]
    thiopeptide.core = query_sequence[end:]

    return thiopeptide


def run_thiopred(query: secmet.CDSFeature, thio_type, domains):

    query_sequence = query.translation

    # Run checks to determine whether an ORF encodes a precursor peptide
    result = determine_precursor_peptide_candidate(query, query_sequence, domains)
    if result is None:
        return

    # Determine thiopeptide type
    result.thio_type = thio_type

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "thiopeptides",
                             "predicted thiopeptide")  # TODO: unsure of this

    # leader cleavage "validation"
    pep_hmmer_profile = 'thiopep2.hmm'
    thresh_pep_hit = -2

    core_a_fasta = ">%s\n%s" % (query.get_name(), result.core)

    profile_pep = path.get_full_path(__file__, "data", pep_hmmer_profile)
    hmmer_res_pep = subprocessing.run_hmmpfam2(profile_pep, core_a_fasta)

    filter_out = True
    for res in hmmer_res_pep:
        for hits in res:
            for seq in hits:
                if seq.bitscore > thresh_pep_hit:
                    filter_out = False

    if filter_out:
        return

    # additional filter(s) for peptide prediction
    aux = ""
    if re.search("[ISTV][SACNTW][STNCVG][ATCSGM][SVTFC][CGSTEAV][TCGVY]", result.core):
        aux = re.search("[ISTV][SACNTW][STNCVG][ATCSGM][SVTFC][CGSTEAV][TCGVY].*",
                        result.core).group()
    else:
        return

    diff = len(result.core)-len(aux)

    if aux != "" and (len(aux) < 20 and len(aux) > 10):
        result.leader = result.leader+result.core[:diff]
        result.core = aux

    # prediction of cleavage in C-terminal based on thiopeptide's core sequence
    # if last core residue != S or T or C > great chance of a tail cut
    if result.core[-1] not in "SCT":
        C_term_hmmer_profile = 'thio_tail.hmm'
        thresh_C_hit = -9

        temp = result.core[-10:]
        core_a_fasta = ">%s\n%s" % (query.get_name(), temp)

        profile_C = path.get_full_path(__file__, "data", C_term_hmmer_profile)
        hmmer_res_C = subprocessing.run_hmmpfam2(profile_C, core_a_fasta)

        for res in hmmer_res_C:
            for hits in res:
                for seq in hits:
                    if seq.bitscore > thresh_C_hit:
                        result.c_cut = temp[seq.query_end-1:]

    if result is None:
        logging.debug('%r: No C-terminal cleavage site predicted', query.get_name())
        return

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "thiopeptides",
                             "predicted thiopeptide")
    return result


class ThiopeptideMotif(secmet.Prepeptide):
    def __init__(self, core_location, core_seq, leader_location, leader_seq,
                 locus_tag, monoisotopic_mass, molecular_weight, alternative_weights,
                 thio_class, score, rodeo_score, macrocycle, cleaved_residues,
                 core_features, mature_weights, amidation):
        super().__init__("thiopeptide", core_location, locus_tag, thio_class,
                         leader=leader_location, leader_seq=leader_seq)
        self.core_seq = core_seq
        self.monoisotopic_mass = monoisotopic_mass
        self.molecular_weight = molecular_weight
        self.alternative_weights = alternative_weights  # list of floats
        self.score = score
        self.rodeo_score = rodeo_score
        self._notes_appended = False
        self.amidation = amidation
        self.macrocycle = macrocycle
        self.cleaved_residues = cleaved_residues
        self.core_features = core_features
        if thio_class == "Type-III":
            assert not mature_weights
        self.mature_weights = mature_weights
        self.tail_reaction = ''
        if self.amidation:
            self.tail_reaction = "dealkylation of C-Terminal residue; amidation"

    def get_modifications(self):
        mods = []
        return mods

    def to_biopython(self):
        if self._notes_appended:  # TODO: could be more clever
            logging.critical("%s already converted: %s, leader type %s",
                             self.location, self._notes_appended, type(self._leader))
            return super().to_biopython()
        self._notes_appended = True
        self.notes.append('monoisotopic mass: %0.1f' % self.monoisotopic_mass)
        self.notes.append('molecular weight: %0.1f' % self.molecular_weight)
        if self.alternative_weights:
            weights = map(lambda x: "%0.1f" % x, self.alternative_weights)
            self.notes.append('alternative weights: %s' % "; ".join(weights))
        self.notes.append('predicted core seq: %s' % self.core_seq)
        self.notes.append('score: %0.2f' % self.score)
        self.notes.append('RODEO score: %s' % str(self.rodeo_score))
        if self.amidation:
            self.notes.append('predicted tail reaction: %s' % self.tail_reaction)
        return super().to_biopython()

    def to_json(self):
        json = dict(vars(self))  # TODO: use a better system that doesn't break encapsulation
        json["core"] = serialiser.location_to_json(self.location)
        json["_leader"] = serialiser.location_to_json(self._leader)
        json["locus_tag"] = self.locus_tag  # not in vars() due to __slots__
        try:
            assert json["locus_tag"]
        except KeyError:
            logging.critical("bad locus tag on motif %s: %s ... %s", self.core, self.locus_tag, json)
        return json

    @staticmethod
    def from_json(data):
        args = []
        args.append(serialiser.location_from_json(data["core"]))
        args.append(data["core_seq"])
        args.append(serialiser.location_from_json(data["_leader"]))
        for arg_name in ["leader_seq", "locus_tag", "monoisotopic_mass",
                         "molecular_weight", "alternative_weights",
                         "peptide_class", "score", "rodeo_score", "macrocycle",
                         "cleaved_residues", "core_features", "mature_weights",
                         "amidation"]:
            args.append(data[arg_name])
        # pylint doesn't do well with the splat op, so don't report errors
        return ThiopeptideMotif(*args)  # pylint: disable=no-value-for-parameter


def result_vec_to_feature(orig_feature, res_vec):
    # resets .. to recalculate weights considering C-terminal putative cut
    res_vec._weight = -1
    res_vec._monoisotopic_weight = -1
    oldcore = res_vec.core
    res_vec.core = oldcore[:(len(res_vec.core) - len(res_vec.c_cut))]

    start = orig_feature.location.start
    end = orig_feature.location.start + (res_vec.end * 3)
    strand = orig_feature.location.strand
    leader_loc = FeatureLocation(start, end, strand=strand)

    start = end
    end = orig_feature.location.end
    core_loc = FeatureLocation(start, end, strand=strand)
    mature_weights = []
    if res_vec.thio_type != "Type-III":
        mature_weights = res_vec.mature_alt_weights
    feature = ThiopeptideMotif(core_loc, res_vec.core, leader_loc, res_vec.leader,
                               orig_feature.get_name(), res_vec.monoisotopic_mass,
                               res_vec.molecular_weight, res_vec.alternative_weights,
                               res_vec.thio_type, res_vec.score, res_vec.rodeo_score,
                               res_vec.macrocycle, res_vec.c_cut, res_vec.mature_features,
                               mature_weights, res_vec.amidation)

    return feature


def specific_analysis(record, options):
    results = ThioResults(record.id)
    for cluster in record.get_clusters():
        if "thiopeptide" not in cluster.products:
            continue

        # Find candidate ORFs that are not yet annotated
        new_orfs = all_orfs.find_all_orfs(record, cluster)

        thio_fs = list(cluster.cds_children) + new_orfs
        domains = get_detected_domains(cluster)
        thio_type = predict_type_from_gene_cluster(domains)

        if thio_type is None:
            return

        amidation = predict_amidation(domains)

        for thio_f in thio_fs:
            result_vec = run_thiopred(thio_f, thio_type, domains)

            if result_vec is None:
                continue

            if amidation:
                result_vec.amidation = True

            new_feature = result_vec_to_feature(thio_f, result_vec)
            if thio_f in new_orfs:
                results.cds_features[cluster.get_cluster_number()].append(thio_f)
            results.motifs.append(new_feature)
            results.clusters_with_motifs.add(cluster)
    logging.debug("Thiopeptides marked %d motifs", len(results.motifs))
    return results
