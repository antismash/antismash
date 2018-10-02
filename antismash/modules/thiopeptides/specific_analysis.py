# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
More detailed thiopeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of molcular mass.
"""

from collections import defaultdict
import logging
import re
import os
from typing import Any, Dict, List, Optional, Set, Tuple

from antismash.common import all_orfs, fasta, module_results, path, secmet, subprocessing, utils
from antismash.common.secmet.qualifiers.prepeptide_qualifiers import ThioQualifier
from antismash.common.signature import HmmSignature

from .rodeo import run_rodeo


class ThioResults(module_results.ModuleResults):
    """ Results container for thiopeptides """
    schema_version = 1

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.clusters_with_motifs = set()  # type: Set[secmet.Cluster]
        # to track CDSs found with find_all_orfs and within which clusters they were found
        self.cds_features = defaultdict(list)  # type: Dict[int, List[secmet.CDSFeature]]
        # to track the motifs created
        self.motifs = []  # type: List[secmet.Prepeptide]

    def to_json(self) -> Dict[str, Any]:
        """ Converts the results to JSON format """
        cds_features_by_cluster = {key: [(str(feature.location), feature.get_name()) for feature in features]
                                   for key, features in self.cds_features.items()}
        return {"record_id": self.record_id,
                "schema_version": ThioResults.schema_version,
                "clusters with motifs": [cluster.get_cluster_number() for cluster in self.clusters_with_motifs],
                "motifs": [motif.to_json() for motif in self.motifs],
                "cds_features": cds_features_by_cluster}

    @staticmethod
    def from_json(json: Dict, record: secmet.Record) -> Optional["ThioResults"]:
        """ Builds a results object from JSON """
        if json.get("schema_version") != ThioResults.schema_version:
            logging.warning("Discarding Thiopeptide results, schema version mismatch")
            return None
        results = ThioResults(json["record_id"])
        for motif in json["motifs"]:
            results.motifs.append(secmet.Prepeptide.from_json(motif))
        for cluster in json["clusters with motifs"]:
            results.clusters_with_motifs.add(record.get_cluster(cluster))
        for cluster, features in json["cds_features"]:
            for location, name in features:
                cds = all_orfs.create_feature_from_location(record, location, label=name)
                results.cds_features[cluster].append(cds)
        return results

    def add_to_record(self, record: secmet.Record) -> None:
        """ Adds any relevant result constructions to the record """
        for features in self.cds_features.values():
            for cds in features:
                record.add_cds_feature(cds)

        for motif in self.motifs:
            record.add_cds_motif(motif)


class Thiopeptide:
    """ Class to calculate and store thiopeptide information
    """
    def __init__(self, end: int, score: float, rodeo_score: int) -> None:
        self.end = end
        self.score = score
        self.rodeo_score = rodeo_score
        self.thio_type = ''
        self._leader = ''
        self._core = ''
        self._weight = -1.
        self._monoisotopic_weight = -1.
        self._alt_weights = []  # type: List[float]
        self._macrocycle = ''
        self.amidation = False
        self._mature_alt_weights = []  # type: List[float]
        self._mature_features = ''
        self._c_cut = ''
        self.core_analysis_monoisotopic = None  # type: Optional[utils.RobustProteinAnalysis]
        self.core_analysis = None  # type: Optional[utils.RobustProteinAnalysis]

    @property
    def core(self) -> str:
        """ The core section of the motif """
        return self._core

    @core.setter
    def core(self, seq: str) -> None:
        assert isinstance(seq, str)
        assert seq
        self._core = seq
        self._weight = -1.
        self._monoisotopic_weight = -1.
        self.core_analysis_monoisotopic = utils.RobustProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = utils.RobustProteinAnalysis(seq, monoisotopic=False)

    @property
    def leader(self) -> str:
        """ The leader section of the motif """
        return self._leader

    @leader.setter
    def leader(self, leader: str) -> None:
        assert isinstance(leader, str)
        self._leader = leader

    @property
    def macrocycle(self) -> str:
        """ Recalculates the macrocycle prediction and returns it """
        self._predict_macrocycle()
        return self._macrocycle

    @macrocycle.setter
    def macrocycle(self, macro: str) -> None:
        assert macro and isinstance(macro, str), macro
        self._macrocycle = macro

    def _predict_macrocycle(self) -> None:
        """ Predict macrocycle ring """
        if not self._core:
            raise ValueError("Predicting macrocycle without a valid core")

        aux = self._core
        if len(self._core) == 17 and self._core[0] in "TIV":
            aux = self._core[4:]

        if aux[0] == "S":
            if "SC" in aux[9:12]:
                if aux[9] == "S":
                    self.macrocycle = "26-member"

                if aux[10] == "S":
                    self.macrocycle = "29-member"

            if re.search(r"S{6,8}?", aux[len(aux)//2:]):
                if aux[12] == "S":
                    self.macrocycle = "35-member"

    @property
    def mature_features(self) -> str:
        """ Recalculates the mature weights of the thiopeptide and returns them """
        self._predict_mature_core_features()
        return self._mature_features

    @mature_features.setter
    def mature_features(self, feats: str) -> None:
        assert isinstance(feats, str)
        self._mature_features = feats

    def _predict_mature_core_features(self) -> str:
        """ Prediction of some features of the mature peptide
        """
        assert self.thio_type

        self.mature_features = ''
        if self.thio_type == 'Type I':
            self.mature_features = "Central ring: pyridine tetrasubstituted (hydroxyl group present); second macrocycle"

        if self.thio_type == 'Type II':
            self.mature_features = "Central ring: piperidine; second macrocycle containing a quinaldic acid moiety"

        if self.thio_type == 'Type III':
            self.mature_features = "Central ring: pyridine trisubstituted"

        return self._mature_features

    @property
    def c_cut(self) -> str:
        """ represents the tail """
        return self._c_cut

    @c_cut.setter
    def c_cut(self, ccut: str) -> None:
        assert isinstance(ccut, str)
        self._c_cut = ccut

    def __repr__(self) -> str:
        return "Thiopeptide(..%s, %s, %r, %r, %s(%s), %s, %s, %s, %s)" % (
                        self.end, self.score, self._core,
                        self.thio_type, self._monoisotopic_weight, self._weight,
                        self.macrocycle, self.amidation, self.mature_features,
                        self.c_cut)

    def _calculate_mw(self) -> None:
        """
        (re)calculate the monoisotopic mass and molecular weight
        """
        assert self._core, "calculating weight without a core"
        assert self.core_analysis and self.core_analysis_monoisotopic, "missing core analyses in weight calculation"

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

        if self.thio_type != "Type III":
            # maturation reactions:
            mol_mass -= dehydrations

            # YcaO > cys/ser to azole ~20 Da for ring
            mol_mass -= (20 * (no_cys - 1))

            # cycloaddition of 2 Dha residues
            mol_mass -= 35

            if self.thio_type == 'Type I':
                # indolic acid (172) + cyclization
                mol_mass += 172 - 3
            elif self.thio_type == 'Type II':
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
    def monoisotopic_mass(self) -> float:
        """ determines the weight of the core peptide and substracts
            the weight which is reduced, due to dehydratation
        """
        if self._monoisotopic_weight > -1:
            return self._monoisotopic_weight

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self) -> float:
        """ determines the weight of the core peptide
        """
        if self._weight > -1:
            return self._weight

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self) -> List[float]:
        """ determines the possible alternative weights assuming one or
            more of the Ser/Thr residues aren't dehydrated
        """
        if self._alt_weights != []:
            return self._alt_weights

        self._calculate_mw()
        return self._alt_weights

    @property
    def mature_alt_weights(self) -> List[float]:
        """ determines the possible mature alternative weights which includes mw,
            monoisotopic mass and alternative weightss due to different hydrated residues
        """
        if self._mature_alt_weights != []:
            return self._mature_alt_weights

        self._calculate_mw()
        return self._mature_alt_weights


def predict_amidation(found_domains: Set[str]) -> bool:
    """ Returns True if amidation likely given the set of domain ids provided
    """
    # nosA homologs nocA, tpdK nocA berI pbt
    return 'thio_amide' in found_domains


def predict_cleavage_site(query_hmmfile: str, target_sequence: str, threshold: float
                          ) -> Tuple[Optional[int], float]:
    """ Extracts the start position, end position and score
        of the HMM alignment from HMMER results.

        Arguments:
            query_hmmfile: the HMM file to search
            target_sequence: the sequence to search
            threshold: the minimum bitscore a hit must have

        Returns:
            a tuple of
                the start of the hit, or None if no hit found
                the end of the hit, or None if no hit found
                the score of the hit, or the best score of all hits if none
                        were above the threshold
    """
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)

    best_score = 0.
    for res in hmmer_res:
        for hits in res:
            for hsp in hits:
                if hsp.bitscore > threshold:
                    return hsp.query_end - 14, hsp.bitscore
                if best_score is None or hsp.bitscore > best_score:
                    best_score = hsp.bitscore

    return None, best_score


def predict_type_from_cluster(found_domains: Set[str]) -> str:
    """ Predict the thiopeptide type from the cluster domains

        Arguments:
            found_domains: the set of domain ids found in the cluster

        Returns:
            the thiopeptide type as a string
    """

    if 'PF06968' in found_domains:
        return 'Type I'

    if 'PF04055' in found_domains or 'PF00733' in found_domains:
        return 'Type II'

    return 'Type III'


def get_detected_domains(cluster: secmet.Cluster) -> Set[str]:
    """ Gathers all detected domain ids from a cluster. Includes detection of
        some extra HMM profiles specific to thiopeptides.

        Arguments:
            cluster: the Cluster to gather domains from

        Return:
            a set of domain ids
    """
    found_domains = []  # type: List[str]
    # Gather biosynthetic domains
    for feature in cluster.cds_children:
        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    # Gather non-biosynthetic domains
    cluster_fasta = fasta.get_fasta_from_features(cluster.cds_children)
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(cluster_fasta)
    non_biosynthetic_hmms_found = []  # type: List[str]
    for hsps_found_for_this_id in non_biosynthetic_hmms_by_id.values():
        for hsp in hsps_found_for_this_id:
            if hsp.query_id not in non_biosynthetic_hmms_found:
                non_biosynthetic_hmms_found.append(hsp.query_id)
    found_domains.extend(non_biosynthetic_hmms_found)

    return set(found_domains)


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


def run_cleavage_site_phmm(input_fasta: str, hmmer_profile: str, threshold: float
                           ) -> Tuple[Optional[int], float]:
    """Try to identify cleavage site using pHMM"""
    profile = path.get_full_path(__file__, "data", hmmer_profile)
    return predict_cleavage_site(profile, input_fasta, threshold)


def run_cleavage_site_regex(sequence: str) -> Optional[int]:
    """Try to identify cleavage site using regular expressions"""
    regex = re.compile('([IVGAS]AS)')

    end = -1

    # For each regular expression, check if there is a match that is <10 AA from the end
    if re.search(regex, sequence) and len(re.split(regex, sequence)[-1]) > 10:
        _, end = [m.span() for m in regex.finditer(sequence)][-1]
        end -= 5

    if end <= 0:
        return None

    return end


def determine_precursor_peptide_candidate(query: secmet.CDSFeature, domains: Set[str]
                                          ) -> Optional[Thiopeptide]:
    """ Identify precursor peptide candidates and split into two

        Arguments:
            query: the CDS feature to check for motifs
            domains: the set of domain ids found in the cluster

        Returns:
            a Thiopeptide instance if a valid precursor found, otherwise None
    """

    query_sequence = query.translation
    # Skip sequences not in the size range desired
    if not 40 < len(query_sequence) < 200:
        return None

    # Create FASTA sequence for feature under study
    thio_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    # Run sequence against pHMM; if positive, parse into a vector containing START, END and SCORE
    end, score = run_cleavage_site_phmm(thio_a_fasta, 'thio_cleave.hmm', -3.00)

    # If no pHMM hit, try regular expression
    if end is None:
        score = 0.
        end = run_cleavage_site_regex(query_sequence)
        if end is None or end > len(query_sequence) - 5:
            end = int(len(query_sequence)*0.60) - 14

    # ensure there's a valid value for end before trying to use it
    assert isinstance(end, int) and end > 0

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(query_sequence[:end], query_sequence[end:], domains)
    if not rodeo_result[0]:
        return Thiopeptide(end + 1, score, 0)

    thiopeptide = Thiopeptide(end + 1, score, rodeo_result[1])

    # Determine the leader and core peptide
    thiopeptide.leader = query_sequence[:end]
    thiopeptide.core = query_sequence[end:]

    return thiopeptide


def find_tail(query: secmet.CDSFeature, core: str) -> str:
    """ Finds the tail of a prepeptide, if it exists

        Arguments:
            query: the CDS feature being checked
            core: the core of the prepeptide as a string

        Returns:
            the translation of the tail, or an empty string if it wasn't found
    """
    # prediction of cleavage in C-terminal based on thiopeptide's core sequence
    # if last core residue != S or T or C > great chance of a tail cut
    tail = ''
    if core[-1] in "SCT":
        return tail
    thresh_c_hit = -9

    temp = core[-10:]
    core_a_fasta = ">%s\n%s" % (query.get_name(), temp)

    c_term_profile = path.get_full_path(__file__, "data", 'thio_tail.hmm')
    c_hmmer_res = subprocessing.run_hmmpfam2(c_term_profile, core_a_fasta)

    for res in c_hmmer_res:
        for hits in res:
            for seq in hits:
                if seq.bitscore > thresh_c_hit:
                    tail = temp[seq.query_end-1:]
    return tail


def run_thiopred(query: secmet.CDSFeature, thio_type: str, domains: Set[str]) -> Optional[Thiopeptide]:
    """ Analyses a CDS feature to determine if it contains a thiopeptide precursor

        Arguments:
            query: the CDS feature to analyse
            thio_type: the suspected type of the thiopeptide
            domains: the set of domains found within the cluster containing the query

        Returns:
            A Thiopeptide instance if a precursor is found, otherwise None
    """
    # Run checks to determine whether an ORF encodes a precursor peptide
    result = determine_precursor_peptide_candidate(query, domains)
    if result is None:
        return None

    # Determine thiopeptide type
    result.thio_type = thio_type

    # leader cleavage "validation"
    profile_pep = path.get_full_path(__file__, "data", 'thiopep2.hmm')
    core_a_fasta = ">%s\n%s" % (query.get_name(), result.core)
    hmmer_res_pep = subprocessing.run_hmmpfam2(profile_pep, core_a_fasta)

    thresh_pep_hit = -2
    filter_out = True
    for res in hmmer_res_pep:
        for hits in res:
            for seq in hits:
                if seq.bitscore > thresh_pep_hit:
                    filter_out = False

    if filter_out:
        return None

    # additional filter(s) for peptide prediction
    search = re.search("[ISTV][SACNTW][STNCVG][ATCSGM][SVTFC][CGSTEAV][TCGVY].*",
                       result.core)
    if not search:
        return None
    aux = search.group()

    if 10 < len(aux) < 20:
        diff = len(result.core) - len(aux)
        result.leader = result.leader + result.core[:diff]
        result.core = aux

    result.c_cut = find_tail(query, result.core)

    query.gene_functions.add(secmet.GeneFunction.ADDITIONAL, "thiopeptides",
                             "predicted thiopeptide")
    return result


def result_vec_to_feature(orig_feature: secmet.CDSFeature, res_vec: Thiopeptide) -> secmet.Prepeptide:
    """ Converts a Thiopeptide object to a Prepeptide, based on an original
        CDSFeature.

        Arguments:
            orig_feature: the original CDS feature that the Motif will attach to
            res_vec: a Thiopeptide object containing results

        Returns:
            a Prepeptide
    """
    if res_vec.c_cut:
        res_vec.core = res_vec.core[:-len(res_vec.c_cut)]

    mature_weights = []  # type: List[float]
    if res_vec.thio_type != "Type III":
        mature_weights = res_vec.mature_alt_weights
    feature = secmet.Prepeptide(orig_feature.location, "thiopeptide", res_vec.core, orig_feature.get_name(),
                                "thiopeptides", res_vec.thio_type, res_vec.score, res_vec.monoisotopic_mass,
                                res_vec.molecular_weight, res_vec.alternative_weights,
                                leader=res_vec.leader, tail=res_vec.c_cut)
    feature.detailed_information = ThioQualifier(res_vec.rodeo_score, res_vec.amidation,
                                                 res_vec.macrocycle, res_vec.mature_features, mature_weights)
    return feature


def specific_analysis(record: secmet.Record) -> ThioResults:
    """ Runs thiopeptide prediction over all cluster features and any extra ORFs
        that are found not overlapping with existing features
    """
    results = ThioResults(record.id)
    for cluster in record.get_clusters():
        if cluster.product != "thiopeptide":
            continue

        # Find candidate ORFs that are not yet annotated
        new_orfs = all_orfs.find_all_orfs(record, cluster)

        thio_features = list(cluster.cds_children) + new_orfs
        domains = get_detected_domains(cluster)
        thio_type = predict_type_from_cluster(domains)

        amidation = predict_amidation(domains)

        for thio_feature in thio_features:
            result_vec = run_thiopred(thio_feature, thio_type, domains)

            if result_vec is None:
                continue

            if amidation:
                result_vec.amidation = True
            new_feature = result_vec_to_feature(thio_feature, result_vec)
            if thio_feature in new_orfs:
                results.cds_features[cluster.get_cluster_number()].append(thio_feature)
            results.motifs.append(new_feature)
            results.clusters_with_motifs.add(cluster)
    logging.debug("Thiopeptides marked %d motifs", len(results.motifs))
    return results
