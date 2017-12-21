# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions to find, classify, and annotate NRPS and PKS domains within CDS
    features.
"""

from typing import Dict, List

from Bio.SeqFeature import FeatureLocation

from antismash.common import path, subprocessing, utils
from antismash.common.fasta import get_fasta_from_features
from antismash.common.hmmscan_refinement import refine_hmmscan_results, HMMResult
from antismash.common.secmet.record import Record
from antismash.common.secmet.feature import AntismashDomain, CDSFeature, CDSMotif


def annotate_domains(record: Record) -> None:
    """ Annotates NRPS/PKS domains on CDS features. The `nrps_pks` member of
        each feature will be updated, along with creating CDSMotif features
        when relevant.

        Arguments:
            record: the secmet.Record of which to annotate CDS features

        Returns:
            None
    """
    cds_within_clusters = record.get_cds_features_within_clusters()
    assert cds_within_clusters  # because every cluster should have genes

    fasta = get_fasta_from_features(cds_within_clusters)
    cds_domains = find_domains(fasta, record)
    cds_motifs = find_ab_motifs(fasta)

    for cds in cds_within_clusters:
        cds_name = cds.get_name()
        # gather domains and classify
        domains = cds_domains.get(cds_name)
        if not domains:
            continue

        # generate domain features
        domain_features = generate_domain_features(record, cds, domains)
        for domain_feature in domain_features:
            record.add_antismash_domain(domain_feature)

        domain_type = classify_feature([domain.hit_id for domain in domains])
        cds.nrps_pks.type = domain_type

        for domain in domains:
            cds.nrps_pks.add_domain(domain)

        # construct motif features
        motifs = cds_motifs.get(cds_name)
        if not motifs:
            continue
        motif_features = generate_motif_features(record, cds, motifs)

        for motif in motif_features:
            record.add_cds_motif(motif)
        cds.motifs.extend(motif_features)


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
    refine_hmmscan_results(domains, lengths, neighbour_mode=True)
    raise NotImplementedError("no return value used from refine_hmmscan_results")


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


def classify_feature(domain_names: List[str]) -> str:
    """ Classifies a CDS based on the type and counts of domains present.

        Arguments:
            domain_names: a list of domain names present in the CDS

        Returns:
            a string of the classification (e.g. 'NRPS-like protein')
    """
    # get the set of domains and count the relevant types
    counter = KetosynthaseCounter(domain_names)
    domains = set(domain_names)

    # which rule does it match
    pks_domains = domains.intersection({"PKS_KS", "PKS_AT"})
    nrps_domains = domains.intersection({"Condensation_LCL", "Condensation_DCL",
                                         "Condensation_Starter", "Cglyc",
                                         "Condensation_Dual", "AMP-binding"})
    if not pks_domains and not nrps_domains:
        classification = "other"
    elif {"Cglyc", "Epimerization", "AMP-binding"} in domains and not pks_domains:
        classification = "Glycopeptide NRPS"
    elif len(nrps_domains) >= 2 and "AMP-binding" in domains:
        if pks_domains:
            classification = "Hybrid PKS-NRPS"
        else:
            classification = "NRPS"
    elif not nrps_domains:
        if {"PKS_KS", "Trans-AT_docking"}.issubset(domains) and "PKS_AT" not in domains and counter.trans_is_greatest():
            classification = "Type I Trans-AT PKS"
        elif len(pks_domains) == 2:
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


def generate_domain_features(record, gene: CDSFeature, domains: List[HMMResult]) -> List[AntismashDomain]:
    new_features = []
    nrat = 0
    nra = 0
    nrcal = 0
    nrkr = 0
    nrXdom = 0
    for domain in domains:
        # calculate respective positions based on aa coordinates
        if gene.location.strand == 1:
            start = gene.location.start + 3 * domain.query_start
            end = gene.location.start + 3 * domain.query_end
        else:
            end = gene.location.end - 3 * domain.query_start
            start = gene.location.end - 3 * domain.query_end
        loc = FeatureLocation(start, end, strand=gene.location.strand)

        # set up new feature
        new_feature = AntismashDomain(loc)
        new_feature.domain = domain.hit_id
        new_feature.locus_tag = gene.locus_tag
        new_feature.detection = "hmmscan"
        new_feature.database = "nrpspksdomains.hmm"
        new_feature.evalue = domain.evalue
        new_feature.score = domain.bitscore

        transl_table = gene.transl_table or 1
        new_feature.translation = str(new_feature.extract(record.seq).translate(table=transl_table))

        if domain.hit_id == "AMP-binding":
            nra += 1
            domainname = "{}_A{}".format(gene.get_name(), nra)
            new_feature.label = domainname
            new_feature.domain_id = "nrpspksdomains_" + domainname
        elif domain.hit_id == "PKS_AT":
            nrat += 1
            domainname = "{}_AT{}".format(gene.get_name(), nrat)
            new_feature.label = domainname
            new_feature.domain_id = "nrpspksdomains_" + domainname
        elif domain.hit_id == "CAL_domain":
            nrcal += 1
            domainname = gene.get_name() + "_CAL" + str(nrcal)
            new_feature.label = domainname
            new_feature.domain_id = "nrpspksdomains_" + domainname
        elif domain.hit_id == "PKS_KR":
            nrkr += 1
            domainname = gene.get_name() + "_KR" + str(nrkr)
            new_feature.label = domainname
            new_feature.domain_id = "nrpspksdomains_" + domainname
        else:
            nrXdom += 1
            new_feature.domain_id = "nrpspksdomains_" + gene.get_name().partition(".")[0] + "_Xdom"+'{:02d}'.format(nrXdom)
        new_features.append(new_feature)
    return new_features


def generate_motif_features(record: Record, feature: CDSFeature, motifs) -> List[CDSMotif]:
    """ Convert a list of HMMResult to a list of CDSMotif features """
    # use a locus tag if one exists
    locus_tag = feature.get_name()
    if feature.locus_tag:
        locus_tag = feature.locus_tag
    # grab the translation table if it's there
    if feature.transl_table:
        transl_table = feature.transl_table
    else:
        transl_table = 1

    motif_features = []
    for i, motif in enumerate(motifs):
        i += 1  # user facing, so 1-indexed
        if feature.location.strand == 1:
            start = feature.location.start + 3 * motif.query_start
            end = feature.location.start + 3 * motif.query_end
        else:
            end = feature.location.end - 3 * motif.query_start
            start = feature.location.end - 3 * motif.query_end
        loc = FeatureLocation(start, end, strand=feature.strand)
        new_motif = CDSMotif(loc)
        new_motif.label = motif.hit_id
        new_motif.motif = motif.hit_id  # TODO: why both label AND motif?
        new_motif.domain_id = 'nrpspksmotif_{}_{:04d}'.format(locus_tag, i)
        new_motif.evalue = motif.evalue
        new_motif.score = motif.bitscore
        new_motif.tool = "pksnrpsmotif"
        new_motif.detection = "hmmscan"
        new_motif.database = "abmotifs"
        new_motif.locus_tag = locus_tag

        new_motif.translation = str(new_motif.extract(record.seq).translate(table=transl_table))
        new_motif.notes.append("NRPS/PKS Motif: %s (e-value: %s, bit-score: %s)" % (
                               motif.hit_id, motif.evalue, motif.bitscore))  # TODO move to CDSMotif

        motif_features.append(new_motif)
    return motif_features
