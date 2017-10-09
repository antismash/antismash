# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
More detailed lanthipeptide analysis using HMMer-based leader peptide
cleavage site prediction as well as prediction of number of lanthionine
bridges and molcular mass.
'''

from collections import defaultdict
import logging
import os
import re

from Bio.SeqFeature import FeatureLocation
from helperlibs.wrappers.io import TemporaryFile
from sklearn.externals import joblib

from antismash.detection.hmm_detection.signatures import HmmSignature
from antismash.common import deprecated, path, subprocessing, secmet, module_results, serialiser
from antismash.config import get_config as get_global_config

from .config import get_config as get_lanthi_config

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
               'Class-III': 3.5}


class LanthiResults(module_results.ModuleResults):
    schema_version = 1

    def __init__(self, record_id, *args):
        super().__init__(record_id, *args)
        self.clusters_with_motifs = set()
        self.motifs = []

    def to_json(self):
        return {"record_id": self.record_id,
                "schema_version": LanthiResults.schema_version,
                "clusters with motifs": [cluster.get_cluster_number() for cluster in self.clusters_with_motifs],
                "motifs": [motif.to_json() for motif in self.motifs]}

    @staticmethod
    def from_json(data, record):
        if data.get("schema_version") != LanthiResults.schema_version:
            logging.warning("Discarding Lanthipeptide results, schema version mismatch")
            return None
        results = LanthiResults(data["record_id"])
        for motif in data["motifs"]:
            results.motifs.append(LanthipeptideMotif.from_json(motif))
        for cluster in data["clusters with motifs"]:
            results.clusters_with_motifs.add(record.get_cluster(cluster))
        return results

    def add_to_record(self, record):
        for motif in self.motifs:
            record.add_cds_motif(motif)


class PrepeptideBase:
    def __init__(self, start, end, score, rodeo_score=None):
        self.start = start  # same as CDS
        self.end = end  # same as CDS
        self.score = score  # of cleavage site
        self.rodeo_score = rodeo_score
        self._leader = None
        self._core = ''
        self._tail = None
        self._lan_bridges = -1
        self._weight = -1
        self._monoisotopic_weight = -1
        self._alt_weights = None
        self.core_analysis_monoisotopic = None
        self.core_analysis = None

    @property
    def core(self):
        return self._core

    @core.setter
    def core(self, seq):
        self.core_analysis_monoisotopic = deprecated.RobustProteinAnalysis(seq, monoisotopic=True)
        self.core_analysis = deprecated.RobustProteinAnalysis(seq, monoisotopic=False)
        self._core = seq
        self._calculate_mw()

    @property
    def leader(self):
        return self._leader

    @leader.setter
    def leader(self, seq):
        self._leader = seq

    def __repr__(self):
        return "LanthipeptideBase(%s..%s, %s, %r, %r)" % (self.start, self.end, self.score, self.rodeo_score, self._core)

    @property
    def number_of_lan_bridges(self):
        '''
        function determines the number of lanthionine bridges in the core peptide
        '''
        raise NotImplementedError()

    def _calculate_mw(self):
        '''
        (re)calculate the monoisotopic mass and molecular weight
        '''
        assert self._core
        raise NotImplementedError()

    @property
    def monoisotopic_mass(self):
        ''' weight of the dehydrated core
        '''
        if self._monoisotopic_weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._monoisotopic_weight

    @property
    def molecular_weight(self):
        ''' Weight of the dehydrated core
        '''
        if self._weight is None:
            raise ValueError("No core to calculate weight of")

        self._calculate_mw()
        return self._weight

    @property
    def alternative_weights(self):
        ''' The possible alternative weights assuming one or more of the Ser/Thr
            residues aren't dehydrated
        '''
        if self._alt_weights is None:
            raise ValueError("No core to calculate weight of")
        return self._alt_weights


class Lanthipeptide(PrepeptideBase):
    '''
    Class to calculate and store lanthipeptide information
    '''
    def __init__(self, start, end, score, rodeo_score, lantype):
        super().__init__(start, end, score, rodeo_score)
        self.lantype = lantype
        self._aminovinyl = False
        self._chlorinated = False
        self._oxygenated = False
        self._lac = False

    def __repr__(self):
        return "Lanthipeptide(%s..%s, %s, %r, %r, %s, %s(%s))" % (self.start,
                    self.end, self.score, self.lantype, self._core,
                    self._lan_bridges, self._monoisotopic_weight, self._weight)

    @property
    def number_of_lan_bridges(self):
        '''
        function determines the number of lanthionine bridges in the core peptide
        '''
        if not self._core:
            raise ValueError("No core to calculate bridges from")

        amino_counts = self.core_analysis.count_amino_acids()
        no_cys = amino_counts['C']
        no_thr_ser = amino_counts['T'] + amino_counts['S']
        self._lan_bridges = min(no_cys, no_thr_ser)
        if self._aminovinyl:
            self._lan_bridges -= 1
        return self._lan_bridges

    def _calculate_mw(self):
        '''
        (re)calculate the monoisotopic mass and molecular weight
        '''
        if not self._core:
            raise ValueError()

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
    def aminovinyl_group(self):
        '''
        Check if lanthipeptide contains an aminovinyl group
        '''
        return self._aminovinyl

    @aminovinyl_group.setter
    def aminovinyl_group(self, value):
        '''
        Define if lanthipeptide contains an aminovinyl group and trigger
        recalculation of the molecular weight if needed
        '''
        self._aminovinyl = value
        if self._core:
            self._calculate_mw()
            # recalculate the number of lan bridges
            self._lan_bridges = -1

    @property
    def chlorinated(self):
        '''
        Check if lanthipeptide is chlorinated
        '''
        return self._chlorinated

    @chlorinated.setter
    def chlorinated(self, value):
        '''
        Define if lanthipeptide is chlorinated and trigger
        recalculation of the molecular weight if needed
        '''
        self._chlorinated = value
        if self._core:
            self._calculate_mw()

    @property
    def oxygenated(self):
        '''
        Check if lanthipeptide is oxygenated
        '''
        return self._oxygenated

    @oxygenated.setter
    def oxygenated(self, value):
        '''
        Define if lanthipeptide is oxygenated and trigger
        recalculation of the molecular weight if needed
        '''
        self._oxygenated = value
        if self._core:
            self._calculate_mw()

    @property
    def lactonated(self):
        '''
        Check if lanthipeptide starts with a lactone
        '''
        return self._lac

    @lactonated.setter
    def lactonated(self, value):
        self._lac = value
        if self._core:
            self._calculate_mw()


class CleavageSiteHit(object):
    def __init__(self, start, end, score, lantype):
        self.start = start
        self.end = end
        self.score = score
        self.lantype = lantype


def get_detected_domains(cluster):
    found_domains = []
    # Gather biosynthetic domains
    for feature in cluster.cds_children:
        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    # Gather non-biosynthetic domains
    cluster_features = cluster.cds_children
    cluster_fasta = deprecated.get_specific_multifasta(cluster_features)
    non_biosynthetic_hmms_by_id = run_non_biosynthetic_phmms(cluster_fasta)
    non_biosynthetic_hmms_found = []
    for hmm in non_biosynthetic_hmms_by_id:
        hsps_found_for_this_id = non_biosynthetic_hmms_by_id[hmm]
        for hsp in hsps_found_for_this_id:
            if hsp.query_id not in non_biosynthetic_hmms_found:
                non_biosynthetic_hmms_found.append(hsp.query_id)
    found_domains += non_biosynthetic_hmms_found

    return found_domains


def run_non_biosynthetic_phmms(fasta):
    """Try to identify cleavage site using pHMM"""
    with open(path.get_full_path(__file__, "data", "non_biosyn_hmms", "hmmdetails.txt"), "r") as handle:
        hmmdetails = [line.strip().split("\t") for line in handle if line.count("\t") == 3]
    _signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]
    non_biosynthetic_hmms_by_id = defaultdict(list)
    for sig in _signature_profiles:
        sig.path = path.get_full_path(__file__, "data", "non_biosyn_hmms", sig.path.rpartition(os.sep)[2])
        runresults = subprocessing.run_hmmsearch(sig.path, fasta)
        for runresult in runresults:
            # Store result if it is above cut-off
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    non_biosynthetic_hmms_by_id[hsp.hit_id].append(hsp)
    return non_biosynthetic_hmms_by_id


def cds_has_domain(cds, domain):
    """Function to test whether a cds has a certain domain"""
    return cds.sec_met and domain in cds.sec_met.domain_ids


def predict_cleavage_site(query_hmmfile, target_sequence, threshold=-100):
    '''
    Function extracts from HMMER the start position, end position and score
    of the HMM alignment
    '''
    hmmer_res = subprocessing.run_hmmpfam2(query_hmmfile, target_sequence)

    for res in hmmer_res:
        for hits in res:
            lanthi_type = hits.description
            for hsp in hits:
                if hsp.bitscore > threshold:
                    return CleavageSiteHit(hsp.query_start - 1, hsp.query_end, hsp.bitscore, lanthi_type)
    return None


def predict_class_from_gene_cluster(cluster):
    '''
    Predict the lanthipeptide class from the gene cluster
    '''
    found_domains = []
    for feature in cluster.cds_children:
        if not feature.sec_met:
            continue
        found_domains.extend(feature.sec_met.domain_ids)

    if 'Lant_dehyd_N' in found_domains or 'Lant_dehyd_C' in found_domains:
        return 'Class-I'
    if 'DUF4135' in found_domains:
        return 'Class-II'
    if 'Pkinase' in found_domains:
        # this could be class 3 or class 4, but as nobody has seen class 4
        # in vivo yet, we'll ignore that
        return 'Class-III'

    # Ok, no biosynthetic enzymes found, let's try the prepeptide
    if 'Gallidermin' in found_domains:
        return 'Class-I'

    return None


def run_cleavage_site_phmm(fasta, hmmer_profile, threshold):
    """Try to identify cleavage site using pHMM"""
    profile = path.get_full_path(__file__, hmmer_profile)
    return predict_cleavage_site(profile, fasta, threshold)


def identify_lanthi_motifs(leader, core):
    """Run FIMO to identify lanthipeptide-specific motifs"""
    motifs_file = path.get_full_path(__file__, "data", "lanthi_motifs_meme.txt")
    with TemporaryFile() as tempfile:
        out_file = open(tempfile.name, "w")
        out_file.write(">query\n%s%s" % (leader, core))
        out_file.close()
        fimo_output = subprocessing.run_fimo_simple(motifs_file, tempfile.name)
    fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
    fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[5]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
    return fimo_motifs, fimo_scores


def run_cleavage_site_regex(fasta):
    """Try to identify cleavage site using regular expressions"""
    # Regular expressions; try 1 first, then 2, etc.
    rex1 = re.compile('F?LD')
    rex2 = re.compile('[LF]?LQ')

    # For regular expression, check if there is a match that is <10 AA from the end
    if re.search(rex1, fasta) and len(re.split(rex1, fasta)[-1]) > 10:
        start, end = [m.span() for m in rex1.finditer(fasta)][-1]
        end += 16
    elif re.search(rex2, fasta) and len(re.split(rex2, fasta)[-1]) > 10:
        start, end = [m.span() for m in rex2.finditer(fasta)][-1]
        end += 15
    else:
        return [None, None, None]
    return start, end, 0


def run_rodeo_svm(csv_columns):
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


def run_rodeo(seq_record, query, leader, core, domains):
    """Run RODEO heuristics + SVM to assess precursor peptide candidate"""
    rodeo_score = 0

    # Incorporate heuristic scores
    heuristic_score, gathered_tabs_for_csv = acquire_rodeo_heuristics(seq_record, query, leader, core, domains)
    rodeo_score += heuristic_score

    fimo_motifs = []
    fimo_scores = {}

    if not get_global_config().without_fimo and get_lanthi_config().fimo_present:
        # Find motifs
        fimo_motifs, fimo_scores = identify_lanthi_motifs(leader, core)

    # Incorporate SVM scores
    csv_columns = generate_rodeo_svm_csv(leader, core, gathered_tabs_for_csv, fimo_motifs, fimo_scores)
    svm_score = run_rodeo_svm(csv_columns)
    rodeo_score += svm_score
    return rodeo_score >= 14, rodeo_score


def lanscout(seq):
    """ define lanthionine ring with a c-terminal Cys """
    lanlower = 2
    lanupper = 6
    lan = '(?=([T|S].{%d,%d}C))' % (lanlower, lanupper)
    sizes = []
    numringlist = []
    seq = [num for elem in seq for num in elem]
    for sub_seq in seq:
        core = str(sub_seq[:])
        totrings = re.compile(lan, re.I).findall(core)
        size = []
        for i in range(0, len(totrings)):
            size.append(len(totrings[i]))
        sizes.append(size)
        numrings = len(totrings)
        numringlist.append(numrings)

    profile = []
    for size in sizes:
        temp = []
        for j in range(lanlower + 2, lanupper + 3):
            temp.append(size.count(j))
        profile.append(temp)

    for i in range(0, len(profile)):
        profile[i] = str(profile[i]).strip('[]')
    return numringlist, profile


def acquire_rodeo_heuristics(seq_record, query, leader, core, domains):
    """Calculate heuristic scores for RODEO"""
    tabs = []
    score = 0
    precursor = leader + core
    # Leader peptide contains FxLD motif
    if re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', leader):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core residue position of Sx4C motif
    if re.search('S[ARNDBCEQZGHILKMFPSTWYV]{4}C', core):
        tabs.append(re.search('S[ARNDBCEQZGHILKMFPSTWYV]{4}C', core).span()[0])
    else:
        tabs.append(0)
    # Core residue position of Tx4C motif
    if re.search('T[ARNDBCEQZGHILKMFPSTWYV]{4}C', core):
        tabs.append(re.search('T[ARNDBCEQZGHILKMFPSTWYV]{4}C', core).span()[0])
    else:
        tabs.append(0)
    # Core residue position of Sx5C motif
    if re.search('S[ARNDBCEQZGHILKMFPSTWYV]{5}C', core):
        tabs.append(re.search('S[ARNDBCEQZGHILKMFPSTWYV]{5}C', core).span()[0])
    else:
        tabs.append(0)
    # Core residue position of Tx5C motif
    if re.search('T[ARNDBCEQZGHILKMFPSTWYV]{5}C', core):
        tabs.append(re.search('T[ARNDBCEQZGHILKMFPSTWYV]{5}C', core).span()[0])
    else:
        tabs.append(0)
    # Precursor is within 500 nt?
    hmmer_profiles = ['LANC_like', 'Lant_dehyd_C']
    distance = deprecated.distance_to_pfam(seq_record, query, hmmer_profiles)
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
    if "Gallidermin" in domains or "mature_a" in domains or "mature_b" in domains or "matura_ab" in domains:
        tabs.append(1)
    else:
        tabs.append(0)
    # Cluster contains PF8130
    if "Antimicr18" in domains:
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide mass < 4000 Da
    precursor_analysis = deprecated.RobustProteinAnalysis(precursor,
                                                          monoisotopic=True,
                                                          ignore_invalid=True)
    if precursor_analysis.molecular_weight() < 4000:
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Core peptide mass < 2000 Da
    core_analysis = deprecated.RobustProteinAnalysis(core, monoisotopic=True,
                                                     ignore_invalid=True)
    if core_analysis.molecular_weight() < 2000:
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide pHMMs below:
    precursor_hit = False
    # Precursor peptide hits gallidermin superfamily (cl03420) HMM
    if cds_has_domain(query, "TIGR03731") or cds_has_domain(query, "Gallidermin"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits lantibio_gallid (TIGR03731) HMM
    if cds_has_domain(query, "TIGR03731"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits lanti_SCO0268 superfamily (cl22812) HMM
    if cds_has_domain(query, "TIGR04451") or cds_has_domain(query, "strep_PEQAXS"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits LD_lanti_pre (TIGR04363) HMM
    if cds_has_domain(query, "LD_lanti_pre"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits Antimicrobial18 (cl06940) HMM
    if cds_has_domain(query, "Antimicr18"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # Precursor peptide hits gallidermin (PF02052) HMM
    if cds_has_domain(query, "Gallidermin") or cds_has_domain(query, "mature_a") or cds_has_domain(query, "mature_ab") or cds_has_domain(query, "mature_b"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)
    # precursor peptide hits Antimicrobial18 (PF08130) HMM
    if cds_has_domain(query, "Antimicr18"):
        precursor_hit = True
        tabs.append(1)
    else:
        tabs.append(0)

    if precursor_hit:
        score += 3

    # Precursor peptide mass (unmodified)
    precursor_analysis = deprecated.RobustProteinAnalysis(precursor, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(precursor_analysis.molecular_weight()))

    # Unmodified leader peptide mass
    leader_analysis = deprecated.RobustProteinAnalysis(leader, monoisotopic=True, ignore_invalid=False)
    tabs.append(float(leader_analysis.molecular_weight()))

    # Unmodified core peptide mass
    core_analysis = deprecated.RobustProteinAnalysis(core, monoisotopic=True, ignore_invalid=False)
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
    if re.search('CC', core[:-3]):
        score -= 3
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader peptide has > 4 negatively charge motifs
    if sum([leader.count(aa) for aa in "DE"]) > 4:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader peptide has net negative charge
    charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
    if sum([charge_dict[aa] for aa in leader if aa in charge_dict]) < 0:
        score += 1
        tabs.append(1)
    else:
        tabs.append(0)
    # Leader residue position of FxLD motif
    if re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', leader):
        tabs.append(re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', leader).span()[0])
    else:
        tabs.append(0)
    # Core peptide contains C-terminal CC (within last 3 residues)
    if re.search('CC', core[-3:]):
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Core peptide contains DGCGxTC / SFNS / SxxLC / CTxGC / TPGC / SFNSxC motifs
    motifs = (('DGCG[ARNDBCEQZGHILKMFPSTWYV]TC', 2), ('SFNS', 2),
              ('S[ARNDBCEQZGHILKMFPSTWYV]{2}LC', 2), ('CT[ARNDBCEQZGHILKMFPSTWYV]{1}GC', 1),
              ('TPGC', 1), ('SFNS[ARNDBCEQZGHILKMFPSTWYV]C', 1))
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
    for aa, penalty in [("C", -10), ("S", -4), ("T", -4)]:
        if aa not in core:
            score += penalty
            tabs.append(1)
        else:
            tabs.append(0)
    # Lanthionine regex maximum ring number > 4
    numringlist, profile = lanscout([[core]])
    if numringlist[0] > 4:
        score += 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Lanthionine regex maximum ring number < 3
    if numringlist[0] < 3:
        score -= 2
        tabs.append(1)
    else:
        tabs.append(0)
    # Lanthionine regex 4-membered ring/5-membered ring/6-membered ring/7-membered ring/8-membered ring
    scores = [2, 2, 2, 2, 1]
    scorepos = 0
    for ringsize in profile[0].split(", ")[:2]:
        if ringsize != "0" and ringsize != "1" and ringsize != "2":
            score += scores[scorepos]
            tabs.append(1)
        else:
            tabs.append(0)
        scorepos += 1
    for ringsize in profile[0].split(", ")[2:]:
        if ringsize != "0":
            score += scores[scorepos]
            tabs.append(1)
        else:
            tabs.append(0)
        scorepos += 1
    return score, tabs


def generate_rodeo_svm_csv(leader, core, previously_gathered_tabs, fimo_motifs, fimo_scores):
    """Generates all the items for one candidate precursor peptide"""
    precursor = leader + core
    columns = []
    # Precursor Index
    columns.append(1)
    # classification
    columns.append(0)
    columns += previously_gathered_tabs
    # Lanthionine regex maximum ring number
    numringlist, profile = lanscout([[core]])
    columns.append(numringlist[0])
    # Lanthionine regex 4-membered ring count
    columns.append(int(profile[0].split(", ")[0]))
    # Lanthionine regex 5-membered ring count
    columns.append(int(profile[0].split(", ")[1]))
    # Lanthionine regex 6-membered ring count
    columns.append(int(profile[0].split(", ")[2]))
    # Lanthionine regex 7-membered ring count
    columns.append(int(profile[0].split(", ")[3]))
    # Lanthionine regex 8-membered ring count
    columns.append(int(profile[0].split(", ")[4]))
    # Ratio of number of Cys in core peptide to sum of Ser/Thr in core peptide
    if "S" in core or "T" in core:
        columns.append(core.count("C") / float(core.count("S") + core.count("T")))
    else:
        columns.append(1.0)
    # Ratio of number of Cys/Ser/Thr to length of core peptide
    columns.append(float(core.count("S") + core.count("T") + core.count("C")) / len(core))
    # log10 p-value MEME motif 1
    if 1 in fimo_motifs:
        columns.append(fimo_scores[1])
    else:
        columns.append(0)
    # log10 p-value MEME motif 2
    if 2 in fimo_motifs:
        columns.append(fimo_scores[2])
    else:
        columns.append(0)
    # log10 p-value MEME motif 3
    if 3 in fimo_motifs:
        columns.append(fimo_scores[3])
    else:
        columns.append(0)
    # log10 p-value MEME motif 4
    if 4 in fimo_motifs:
        columns.append(fimo_scores[4])
    else:
        columns.append(0)
    # log10 p-value MEME motif 5
    if 5 in fimo_motifs:
        columns.append(fimo_scores[5])
    else:
        columns.append(0)
    # Number in leader of each amino acid
    columns += [leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in leader of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    columns.append(sum([leader.count(aa) for aa in "FWY"]))
    columns.append(sum([leader.count(aa) for aa in "DE"]))
    columns.append(sum([leader.count(aa) for aa in "RK"]))
    columns.append(sum([leader.count(aa) for aa in "RKDE"]))
    columns.append(sum([leader.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([leader.count(aa) for aa in "ST"]))
    # Number in core of each amino acid
    columns += [core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in core of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    columns.append(sum([core.count(aa) for aa in "FWY"]))
    columns.append(sum([core.count(aa) for aa in "DE"]))
    columns.append(sum([core.count(aa) for aa in "RK"]))
    columns.append(sum([core.count(aa) for aa in "RKDE"]))
    columns.append(sum([core.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([core.count(aa) for aa in "ST"]))
    # Number in entire precursor of each amino acid
    columns += [precursor.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
    # Number in entire precursor of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
    columns.append(sum([precursor.count(aa) for aa in "FWY"]))
    columns.append(sum([precursor.count(aa) for aa in "DE"]))
    columns.append(sum([precursor.count(aa) for aa in "RK"]))
    columns.append(sum([precursor.count(aa) for aa in "RKDE"]))
    columns.append(sum([precursor.count(aa) for aa in "GAVLMI"]))
    columns.append(sum([precursor.count(aa) for aa in "ST"]))
    return columns


def determine_precursor_peptide_candidate(seq_record, query, query_sequence, domains, hmmer_profile):
    """ Identify precursor peptide candidates and split into two,
        only valid for Class-I lanthipeptides
    """

    # Skip sequences with >200 AA
    if len(query_sequence) > 200 or len(query_sequence) < 20:
        return

    # Create FASTA sequence for feature under study
    lan_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    # Run sequence against pHMM; if positive, parse into a vector containing START, END and SCORE
    cleavage_result = run_cleavage_site_phmm(lan_a_fasta, hmmer_profile, THRESH_DICT["Class-I"])

    if cleavage_result is not None and cleavage_result.end <= len(query_sequence) - 8:
        start = cleavage_result.start
        end = cleavage_result.end
        score = cleavage_result.score
        lanthi_type = cleavage_result.lantype
    else:
        # If no pHMM hit, try regular expression
        start, end, score = run_cleavage_site_regex(lan_a_fasta)
        if score is None or end > len(query_sequence) - 8:
            start, end, score = 0, int(len(query_sequence)*0.50), 0
        lanthi_type = "lanthipeptide"

    # Run RODEO to assess whether candidate precursor peptide is judged real
    rodeo_result = run_rodeo(seq_record, query, query_sequence[:end], query_sequence[end:], domains)
    if rodeo_result[0] is False:
        return
    else:
        lanthipeptide = Lanthipeptide(start, end, score, rodeo_result[1], lanthi_type)

    # Determine the leader and core peptide
    lanthipeptide.leader = query_sequence[:end]
    lanthipeptide.core = query_sequence[end:]

    return lanthipeptide


def run_lanthipred(seq_record, query, lant_class, domains):
    hmmer_profiles = {'Class-I': 'data/class1.hmm',
                      'Class-II': 'data/class2.hmm',
                      'Class-III': 'data/class3.hmm', }
    query_sequence = query.get_aa_sequence(to_stop=True)
    lan_a_fasta = ">%s\n%s" % (query.get_name(), query_sequence)

    if lant_class in ("Class-II", "Class-III"):
        profile = path.get_full_path(__file__, hmmer_profiles[lant_class])
        cleavage_result = predict_cleavage_site(profile, lan_a_fasta)

        if cleavage_result is None:
            #logging.debug('%s: No cleavage site predicted.', query.get_name())
            return None

        if THRESH_DICT[lant_class] > cleavage_result.score:
            #logging.debug('%r: Score %0.2f below threshold %0.2f for class %r',
            #              query.get_name(), cleavage_result.score,
            #              THRESH_DICT[lant_class], lant_class)
            return None

        result = Lanthipeptide(cleavage_result.start, cleavage_result.end, cleavage_result.score, "N/A", lant_class)
        result.leader = query_sequence[:result.end]
        result.core = query_sequence[result.end:]

    else:
        result = determine_precursor_peptide_candidate(seq_record, query, query_sequence, domains, hmmer_profiles[lant_class])
        if result is None:
            return None

    # extract now (that class is known and thus the END component) the core peptide
    if result.core.find('C') < 0:
        #logging.debug('%r: No Cysteine residues found in core, false positive',
        #              utils.get_gene_id(query))
        return None

    if not query.gene_function:
        query.gene_function = secmet.GeneFunction.ADDITIONAL

    return result


def find_lan_a_features(cluster):
    lan_a_features = []
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue

        aa_seq = feature.get_aa_sequence()
        if len(aa_seq) < 80:
            lan_a_features.append(feature)
            continue
        if feature.sec_met and set(feature.sec_met.domain_ids).intersection(KNOWN_PRECURSOR_DOMAINS):
            lan_a_features.append(feature)

    return lan_a_features


def has_only_domain(feature, domain):
    return feature.sec_met and feature.sec_met.domain_ids == [domain]


def find_flavoprotein(cluster):
    "Look for an epiD-like flavoprotein responsible for aminovinylcystein"
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue
        if has_only_domain(feature, "Flavoprotein"):
            return True
    return False


def find_halogenase(cluster):
    "Look for a halogenase"
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue
        if has_only_domain(feature, "Trp_halogenase"):
            return True
    return False


def find_p450_oxygenase(cluster):
    "Look for a p450 oxygenase"
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue
        if has_only_domain(feature, 'p450'):
            return True
    return False


def find_short_chain_dehydrogenase(cluster):
    "Look for an eciO-like short-chain dehydrogenase responsible for N-terminal lactone"
    for feature in cluster.cds_children:
        if not feature.is_contained_by(cluster):
            continue
        if has_only_domain(feature, 'adh_short') or has_only_domain(feature, 'adh_short_C2'):
            return True
    return False


class LanthipeptideMotif(secmet.Prepeptide):
    def __init__(self, core_location, core_seq, leader_location, leader_seq,
                 locus_tag, monoisotopic_mass, molecular_weight, alternative_weights,
                 lan_bridges, lanthi_class, score, rodeo_score, aminovinyl,
                 chlorinated, oxygenated, lactonated):
        super().__init__("lanthipeptide", core_location, locus_tag, lanthi_class,
                         leader=leader_location, leader_seq=leader_seq)
        self.core_seq = core_seq
        self.monoisotopic_mass = monoisotopic_mass
        self.molecular_weight = molecular_weight
        self.alternative_weights = alternative_weights  # list of floats
        self.lan_bridges = lan_bridges
        self.score = score
        self.rodeo_score = rodeo_score
        self.aminovinyl_group = aminovinyl  # bool
        self.chlorinated = chlorinated  # bool
        self.oxygenated = oxygenated  # bool
        self.lactonated = lactonated  # bool
        self._notes_appended = False

    def get_modifications(self):
        mods = []
        if self.aminovinyl_group:
            mods.append("AviCys")
        if self.chlorinated:
            mods.append("Cl")
        if self.oxygenated:
            mods.append("OH")
        if self.lactonated:
            mods.append("Lac")
        return mods

    def to_biopython(self):
        logging.critical("%s already converted: %s, leader type %s", self.location, self._notes_appended, type(self._leader))
        if self._notes_appended:  # TODO: could be more clever
            return super().to_biopython()
        self._notes_appended = True
        self.notes.append('monoisotopic mass: %0.1f' % self.monoisotopic_mass)
        self.notes.append('molecular weight: %0.1f' % self.molecular_weight)
        if self.alternative_weights:
            weights = map(lambda x: "%0.1f" % x, self.alternative_weights)
            self.notes.append('alternative weights: %s' % "; ".join(weights))
        self.notes.append('number of bridges: %s' % self.lan_bridges)
        self.notes.append('predicted core seq: %s' % self.core_seq)
        self.notes.append('score: %0.2f' % self.score)
        self.notes.append('RODEO score: %s' % str(self.rodeo_score))
        if self.aminovinyl_group:
            self.notes.append('predicted additional modification: AviCys')
        if self.chlorinated:
            self.notes.append('predicted additional modification: Cl')
        if self.oxygenated:
            self.notes.append('predicted additional modification: OH')
        if self.lactonated:
            self.notes.append('predicted additional modification: Lac')
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
                         "molecular_weight", "alternative_weights", "lan_bridges",
                         "peptide_class", "score", "rodeo_score", "aminovinyl_group",
                         "chlorinated", "oxygenated", "lactonated"]:
            args.append(data[arg_name])
# pylint doesn't do well with the splat op, so don't report errors
# pylint: disable=no-value-for-parameter
        return LanthipeptideMotif(*args)
# pylint: enable=no-value-for-parameter


def result_vec_to_feature(orig_feature, res_vec):
    start = orig_feature.location.start
    end = orig_feature.location.start + (res_vec.end * 3)
    strand = orig_feature.location.strand
    leader_loc = FeatureLocation(start, end, strand=strand)

    start = end
    end = orig_feature.location.end
    core_loc = FeatureLocation(start, end, strand=strand)
    feature = LanthipeptideMotif(core_loc, res_vec.core, leader_loc, res_vec.leader,
                                 orig_feature.get_name(), res_vec.monoisotopic_mass,
                                 res_vec.molecular_weight, res_vec.alternative_weights,
                                 res_vec.number_of_lan_bridges, res_vec.lantype,
                                 res_vec.score, res_vec.rodeo_score, res_vec.aminovinyl_group,
                                 res_vec.chlorinated, res_vec.oxygenated, res_vec.lactonated)
    return feature


def specific_analysis(seq_record):
    results = LanthiResults(seq_record.id)
    for cluster in seq_record.get_clusters():
        if 'lanthipeptide' not in cluster.products:
            continue

        lan_as = find_lan_a_features(cluster)

        # Find candidate ORFs that are not yet annotated
        for orf in deprecated.find_all_orfs(seq_record, cluster):
            aa_seq = orf.get_aa_sequence()
            if len(aa_seq) < 80:
                lan_as.append(orf)

        domains = get_detected_domains(cluster)
        flavoprotein_found = find_flavoprotein(cluster)
        halogenase_found = find_halogenase(cluster)
        oxygenase_found = find_p450_oxygenase(cluster)
        dehydrogenase_found = find_short_chain_dehydrogenase(cluster)

        lant_class = predict_class_from_gene_cluster(cluster)
        if not lant_class:
            continue
        for lan_a in lan_as:
            result_vec = run_lanthipred(seq_record, lan_a, lant_class, domains)
            if result_vec is None:
                continue
            result_vec.aminovinyl_group = flavoprotein_found
            result_vec.chlorinated = halogenase_found
            result_vec.oxygenated = oxygenase_found
            result_vec.lactonated = dehydrogenase_found and result_vec.core.startswith('S')
            motif = result_vec_to_feature(lan_a, result_vec)
            results.motifs.append(motif)
            results.clusters_with_motifs.add(cluster)
            lan_a.gene_function = secmet.GeneFunction.ADDITIONAL
            if "allorf" in lan_a.get_name():
                seq_record.add_cds_feature(lan_a)  # TODO shift to add_to_record?
                if lan_a.location.start < cluster.location.start:
                    logging.critical("Cluster location being altered in lanthipeptides")
                    cluster.location = FeatureLocation(lan_a.location.start, cluster.location.end)
                if lan_a.location.end > cluster.location.end:
                    logging.critical("Cluster location being altered in lanthipeptides")
                    cluster.location = FeatureLocation(cluster.location.start, lan_a.location.end)
    logging.debug("Lanthipeptide module marked %d motifs", len(results.motifs))
    return results
