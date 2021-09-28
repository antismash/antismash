# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions to find, classify, and annotate NRPS and PKS domains within CDS
    features.
"""

from collections import defaultdict
import logging
from typing import Any, Dict, List, Optional

from antismash.common import module_results, path, subprocessing, utils
from antismash.common.fasta import get_fasta_from_features
from antismash.common.hmmscan_refinement import refine_hmmscan_results, HMMResult
from antismash.common.secmet.record import Record
from antismash.common.secmet.features import (
    AntismashDomain,
    CDSFeature,
    CDSMotif,
    Module as ModuleFeature,
)
from antismash.common.secmet.locations import FeatureLocation

from .module_identification import (
    build_modules_for_cds,
    CDSModuleInfo,
    combine_modules,
    IncompatibleComponentError,
    Module,
)
from .modular_domain import ModularDomain


class CDSResult:
    """ Stores and enables reconstruction of all results for a single CDS """
    def __init__(self, domain_hmms: List[HMMResult], motif_hmms: List[HMMResult],
                 feature_type: str, modules: List[Module], ks_subtypes: List[str]) -> None:
        self.domain_hmms = domain_hmms
        self.motif_hmms = motif_hmms
        self.type = feature_type
        self.modules = modules
        self.ks_subtypes = ks_subtypes
        self.domain_features: Dict[HMMResult, ModularDomain] = {}
        ks_count = sum(dom.hit_id == "PKS_KS" for dom in self.domain_hmms)
        if len(ks_subtypes) != ks_count:
            raise ValueError("mismatching KS subtypes and PKS_KS counts: %d  %d" % (len(ks_subtypes), ks_count))

    def to_json(self) -> Dict[str, Any]:
        """ Create a JSON representation """
        return {
            "domain_hmms": [hmm.to_json() for hmm in self.domain_hmms],
            "motif_hmms": [hmm.to_json() for hmm in self.motif_hmms],
            "modules": [module.to_json() for module in self.modules],
            "type": self.type,
            "ks_subtypes": self.ks_subtypes,
        }

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "CDSResult":
        """ Reconstruct from a JSON representation """
        domain_hmms = [HMMResult.from_json(hmm) for hmm in data["domain_hmms"]]
        motif_hmms = [HMMResult.from_json(hmm) for hmm in data["motif_hmms"]]
        modules = [Module.from_json(module) for module in data["modules"]]
        return CDSResult(domain_hmms, motif_hmms, data["type"], modules, data["ks_subtypes"])

    def annotate_domains(self, record: Record, cds: CDSFeature) -> None:
        """ Adds domain annotations to CDSFeatures and creates ModularDomain
            features for all domains found
        """
        if not self.domain_hmms:
            return

        cds.nrps_pks.type = self.type

        # generate domain features
        self.domain_features = generate_domain_features(cds, self.domain_hmms)
        ks_sub = iter(self.ks_subtypes)
        for domain, domain_feature in self.domain_features.items():
            if domain.hit_id == "PKS_KS":
                sub = next(ks_sub)
                domain_feature.domain_subtype = sub
            else:
                sub = ""
            record.add_antismash_domain(domain_feature)
            # update the CDS' NRPS_PKS qualifier
            cds.nrps_pks.add_domain(domain, domain_feature.get_name(), sub)

        # construct CDSMotif features
        if not self.motif_hmms:
            return

        motif_features = generate_motif_features(cds, self.motif_hmms)

        for motif in motif_features:
            record.add_cds_motif(motif)
        cds.motifs.extend(motif_features)


class NRPSPKSDomains(module_results.DetectionResults):
    """ Results tracking for NRPS and PKS domains """
    schema_version = 3

    def __init__(self, record_id: str, cds_results: Dict[CDSFeature, CDSResult] = None) -> None:
        super().__init__(record_id)
        if cds_results is None:
            cds_results = {}
        self.cds_results = cds_results

    def to_json(self) -> Dict[str, Any]:
        return {"cds_results": {cds.get_name(): cds_result.to_json() for cds, cds_result in self.cds_results.items()},
                "schema_version": NRPSPKSDomains.schema_version,
                "record_id": self.record_id}

    def add_to_record(self, record: Record) -> None:
        # track multi-CDS modules to avoid duplication
        added_modules = set()
        for cds, result in self.cds_results.items():
            for module in result.modules:
                if module in added_modules:
                    continue
                added_modules.add(module)
                domains: List[AntismashDomain] = []
                for component in module:
                    if component.locus == cds.get_name():
                        domain = result.domain_features[component.domain]
                    else:
                        other_cds_results = self.cds_results[record.get_cds_by_name(component.locus)]
                        domain = other_cds_results.domain_features[component.domain]
                    domains.append(domain)
                mod_type = ModuleFeature.types.UNKNOWN
                if module.is_nrps():
                    mod_type = ModuleFeature.types.NRPS
                elif module.is_pks():
                    mod_type = ModuleFeature.types.PKS
                feature = ModuleFeature(domains, mod_type, complete=module.is_complete(),
                                        starter=module.is_starter_module(),
                                        final=module.is_termination_module(),
                                        iterative=module.is_iterative())
                record.add_module(feature)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["NRPSPKSDomains"]:
        if NRPSPKSDomains.schema_version != json.get("schema_version"):
            logging.warning("Schema version mismatch, discarding NRPS PKS domain results")
            return None
        if record.id != json.get("record_id"):
            logging.warning("Record identifier mismatch, discarding NRPS PKS domain results")
            return None

        cds_results = {}
        for cds_name, cds_result in json["cds_results"].items():
            cds = record.get_cds_by_name(cds_name)
            cds_result = CDSResult.from_json(cds_result)
            cds_result.annotate_domains(record, cds)
            cds_results[cds] = cds_result

        return NRPSPKSDomains(record.id, cds_results)


def match_subtypes_to_ks_domains(domains: List[HMMResult], subtypes: List[HMMResult]) -> List[str]:
    """ Returns a subtype name for each PKS_KS domain in the domains provided.
        If no subtype hit overlaps with a PKS_KS domain, the subtype name will be the
        empty string.

        Arguments:
            domains: a list of HMMResults, non-PKS_KS domains will be ignored
            subtypes: a list of PKS_KS subtypes

        Returns:
            a list of strings, one for each PKS_KS domain
    """
    domains = [dom for dom in domains if dom.hit_id == "PKS_KS"]
    if not domains:
        return []
    subs = []
    for domain in domains:
        sub = ""
        for subtype in subtypes:
            if domain.query_end >= subtype.query_start and subtype.query_end >= domain.query_start:
                sub = subtype.hit_id
                break
        subs.append(sub)
    assert len(domains) == len(subs)
    return subs


def generate_domains(record: Record) -> NRPSPKSDomains:
    """ Annotates NRPS/PKS domains on CDS features. The `nrps_pks` member of
        each feature will be updated, along with creating CDSMotif features
        when relevant.

        Arguments:
            record: the secmet.Record of which to annotate CDS features

        Returns:
            a NRPSPKSDomains instance containing all found motifs and domain HMMs for each CDS
    """
    results = NRPSPKSDomains(record.id)

    cds_within_regions = record.get_cds_features_within_regions()
    assert cds_within_regions  # because every cluster should have genes

    fasta = get_fasta_from_features(cds_within_regions)
    cds_domains = find_domains(fasta, record)
    cds_ks_subtypes = find_ks_domains(fasta)
    cds_motifs = find_ab_motifs(fasta)

    prev: Optional[CDSModuleInfo] = None
    for cds in cds_within_regions:
        domains = cds_domains.get(cds.get_name(), [])
        motifs = cds_motifs.get(cds.get_name(), [])
        if not (domains or motifs):
            continue
        subtype_names = match_subtypes_to_ks_domains(domains, cds_ks_subtypes.get(cds.get_name(), []))
        domain_type = classify_cds([domain.hit_id for domain in domains], subtype_names)
        modules = build_modules_for_cds(domains, subtype_names, cds.get_name())
        results.cds_results[cds] = CDSResult(domains, motifs, domain_type, modules, subtype_names)

        # combine modules that cross CDS boundaries, if possible and relevant
        info = CDSModuleInfo(cds, modules)
        if prev and prev.modules and info.modules:
            combine_modules(info, prev)  # modifies the lists of modules linked in each CDSResult
        prev = info

    for cds, cds_result in results.cds_results.items():
        # filter out modules with only a single component, they're just noise at this stage
        cds_result.modules = list(filter(lambda mod: len(mod.components) > 1, cds_result.modules))
        cds_result.annotate_domains(record, cds)
    return results


def filter_nonterminal_docking_domains(record: Record, cds_domains: Dict[str, List[HMMResult]]
                                       ) -> Dict[str, List[HMMResult]]:
    """ For multiprotein domains, remove all docking terminal predictions that
        aren't overlapping with the first or last 50 amino acids of the protein.
    """
    dockingdomains = {'NRPS-COM_Nterm', 'NRPS-COM_Cterm',
                      'PKS_Docking_Cterm', 'PKS_Docking_Nterm'}
    feature_by_id = record.get_cds_name_mapping()
    results = {}
    for cds_name in list(cds_domains):
        new = []
        cds_length = len(feature_by_id[cds_name].translation)
        for hit in cds_domains[cds_name]:
            if hit.hit_id in dockingdomains and \
                    not (cds_length - max(hit.query_start, hit.query_end) < 50
                         or min(hit.query_start, hit.query_end) < 50):
                continue
            new.append(hit)
        if new:
            results[cds_name] = new
    return results


def find_ab_motifs(fasta: str) -> Dict[str, List[HMMResult]]:
    """ Analyse for abMotifs

        Arguments:
            fasta: a group of features in fasta format

        Returns:
            a dictionary mapping feature name to a list of motif results for that feature
    """
    opts = ["-E", "0.25"]
    motif_file = path.get_full_path(__file__, "data", "abmotifs.hmm")
    abmotif_results = subprocessing.run_hmmscan(motif_file, fasta, opts)
    lengths = utils.get_hmm_lengths(motif_file)
    return refine_hmmscan_results(abmotif_results, lengths, neighbour_mode=True)


def find_domains(fasta: str, record: Record) -> Dict[str, List[HMMResult]]:
    """ Analyse for C/A/PCP/E/KS/AT/ATd/DH/KR/ER/ACP/TE/TD/COM/Docking/MT/CAL domains

        Arguments:
            fasta: a group of features in fasta format
            record: the Record that contains all the features

        Returns:
            a dictionary mapping feature name to a list of domain results for that feature
    """
    opts = ["--cut_tc"]
    nrpspks_file = path.get_full_path(__file__, "data", "nrpspksdomains.hmm")
    nrpspksdomain_results = subprocessing.run_hmmscan(nrpspks_file, fasta, opts)
    lengths = utils.get_hmm_lengths(nrpspks_file)
    domains = refine_hmmscan_results(nrpspksdomain_results, lengths, neighbour_mode=True)
    return filter_nonterminal_docking_domains(record, domains)


def find_ks_domains(fasta: str) -> Dict[str, List[HMMResult]]:
    """ Analyse KS domains & PKS/NRPS protein domain composition to detect NRPS/PKS types

        Arguments:
            fasta: a group of features in fasta format

        Returns:
            a dictionary mapping feature name to a list of KS domain results for that feature
    """
    opts = ["--cut_tc"]
    ks_file = path.get_full_path(__file__, "data", "ksdomains.hmm")
    lengths = utils.get_hmm_lengths(ks_file)
    domains = subprocessing.run_hmmscan(ks_file, fasta, opts)
    return refine_hmmscan_results(domains, lengths, neighbour_mode=True)


class KetosynthaseCounter:
    """ Keeps track of the counts of various KS domains and simplifies
        finding the largest value.
    """
    def __init__(self, domain_names: List[str]) -> None:
        """
            Arguments:
                domain_names: a collection of domain names
        """
        self.modular = 0
        self.trans_at = 0
        self.enediyne = 0
        self.iterative = 0
        self.pks = 0

        for domain in domain_names:
            if domain == "PKS_KS":
                self.pks += 1
            elif domain == "Trans-AT-KS":
                self.trans_at += 1
            elif domain == "Modular-KS":
                self.modular += 1
            elif domain == "Enediyne-KS":
                self.enediyne += 1
            elif domain == "Iterative-KS":
                self.iterative += 1

    def trans_is_greatest(self) -> bool:
        """ Returns true if the trans_at count is strictly greater than others """
        return self.trans_at > max([self.modular, self.enediyne, self.iterative])

    def ene_is_greatest(self) -> bool:
        """ Returns true if the enediyne count is strictly greater than others """
        return self.enediyne > max([self.modular, self.trans_at, self.iterative])

    def modular_is_greatest(self) -> bool:
        """ Returns true if the modular count is strictly greater than others """
        return self.modular > max([self.enediyne, self.trans_at, self.iterative])

    def iterative_is_greatest(self) -> bool:
        """ Returns true if the iterative count is strictly greater than others """
        return self.iterative > max([self.enediyne, self.trans_at, self.modular])


def classify_cds(domain_names: List[str], ks_domain_subtypes: List[str]) -> str:
    """ Classifies a CDS based on the type and counts of domains present.

        Arguments:
            domain_names: a list of domain names present in the CDS

        Returns:
            a string of the classification (e.g. 'NRPS-like protein')
    """
    # get the set of domains and count the relevant types
    counter = KetosynthaseCounter(domain_names + ks_domain_subtypes)
    domains = set(domain_names)

    # which rule does it match
    pks_domains = domains.intersection({"PKS_KS", "PKS_AT"})
    nrps_domains = domains.intersection({"Condensation_LCL", "Condensation_DCL",
                                         "Condensation_Starter", "Cglyc",
                                         "Condensation_Dual", "AMP-binding"})
    if not pks_domains and not nrps_domains:
        classification = "other"
    elif {"Cglyc", "Epimerization", "AMP-binding"}.issubset(domains) and not pks_domains:
        classification = "Glycopeptide NRPS"
    elif len(nrps_domains) >= 2 and "AMP-binding" in domains:
        if pks_domains:
            classification = "Hybrid PKS-NRPS"
        else:
            classification = "NRPS"
    elif not nrps_domains:
        if {"PKS_KS", "Trans-AT_docking"}.issubset(domains) and "PKS_AT" not in domains and counter.trans_is_greatest():
            classification = "Type I Trans-AT PKS"
        elif len(pks_domains) == 2:  # are both KS and AT domains present
            if counter.iterative_is_greatest() and counter.pks < 3:
                classification = "Type I Iterative PKS"
            elif counter.ene_is_greatest() and counter.pks < 3:
                classification = "Type I Enediyne PKS"
            elif counter.modular_is_greatest() or counter.pks > 3:
                classification = "Type I Modular PKS"
            else:
                classification = "PKS-like protein"
        else:
            classification = "PKS/NRPS-like protein"
    elif not pks_domains:
        classification = "NRPS-like protein"
    else:
        classification = "PKS/NRPS-like protein"
    return classification


def generate_domain_features(gene: CDSFeature, domains: List[HMMResult]) -> Dict[HMMResult, ModularDomain]:
    """ Generates ModularDomain features for each provided HMMResult

        Arguments:
            gene: the CDSFeature the domains were found in
            domains: a list of HMMResults found in the CDSFeature

        Returns:
            a dictionary mapping the HMMResult used to the matching ModularDomain
    """
    new_features = {}
    domain_counts: Dict[str, int] = defaultdict(int)
    for domain in domains:
        loc = gene.get_sub_location_from_protein_coordinates(domain.query_start, domain.query_end)
        prot_loc = FeatureLocation(domain.query_start, domain.query_end)

        # set up new feature
        new_feature = ModularDomain(loc, protein_location=prot_loc,
                                    locus_tag=gene.get_name())
        new_feature.domain = domain.hit_id
        new_feature.locus_tag = gene.locus_tag or gene.get_name()
        new_feature.detection = "hmmscan"
        new_feature.database = "nrpspksdomains.hmm"
        new_feature.evalue = domain.evalue
        new_feature.score = domain.bitscore

        new_feature.translation = gene.translation[domain.query_start:domain.query_end]

        domain_counts[domain.hit_id] += 1  # 1-indexed, so increment before use
        domain_name = "{}_{}.{}".format(gene.get_name(), domain.hit_id, domain_counts[domain.hit_id])

        new_feature.domain_id = "nrpspksdomains_" + domain_name
        new_feature.label = domain_name

        new_features[domain] = new_feature
    return new_features


def generate_motif_features(feature: CDSFeature, motifs: List[HMMResult]) -> List[CDSMotif]:
    """ Convert a list of HMMResult to a list of CDSMotif features """
    # use a locus tag if one exists
    locus_tag = feature.get_name()
    if feature.locus_tag:
        locus_tag = feature.locus_tag

    motif_features = []
    for i, motif in enumerate(motifs):
        i += 1  # user facing, so 1-indexed
        loc = feature.get_sub_location_from_protein_coordinates(motif.query_start, motif.query_end)
        prot_loc = FeatureLocation(motif.query_start, motif.query_end)
        new_motif = CDSMotif(loc, feature.get_name(), prot_loc, tool="nrps_pks_domains")
        new_motif.label = motif.hit_id
        new_motif.domain_id = 'nrpspksmotif_{}_{:04d}'.format(locus_tag, i)
        new_motif.evalue = motif.evalue
        new_motif.score = motif.bitscore
        new_motif.detection = "hmmscan"
        new_motif.database = "abmotifs"
        new_motif.locus_tag = locus_tag

        new_motif.translation = feature.translation[motif.query_start:motif.query_end]

        motif_features.append(new_motif)
    return motif_features
