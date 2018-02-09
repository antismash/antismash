# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

from Bio.SeqFeature import FeatureLocation

from antismash.common import path, subprocessing, utils
from antismash.common.fasta import get_fasta_from_features
from antismash.common.hmmscan_refinement import refine_hmmscan_results
from antismash.common.secmet.feature import CDSMotif


def annotate_domains(record) -> None:
    """ Annotates NRPS/PKS domains on CDS features. The `nrps_pks` member of
        each feature will be updated, along with creating CDSMotif features
        when relevant.

        Arguments:
            record: the secmet.Record of which to annotate CDS features

        Returns:
            None
    """
    genes_within_clusters = record.get_cds_features_within_clusters()
    assert genes_within_clusters  # because every cluster should have genes

    fasta = get_fasta_from_features(genes_within_clusters)
    gene_domains = find_domains(fasta, record)
    gene_motifs = find_ab_motifs(fasta)

    for gene in genes_within_clusters:
        gene_name = gene.get_name()
        # gather domains and classify
        domains = gene_domains.get(gene_name)
        if not domains:
            continue
        domain_type = classify_feature([domain.hit_id for domain in domains])
        gene.nrps_pks.type = domain_type

        for domain in domains:
            gene.nrps_pks.add_domain(domain)

        # construct motif features
        motifs = gene_motifs.get(gene_name)
        if not motifs:
            continue
        motif_features = generate_motif_features(record, gene, motifs)

        for motif in motif_features:
            record.add_cds_motif(motif)
        gene.motifs.extend(motif_features)


def filter_nonterminal_docking_domains(record, gene_domains):
    """ For multiprotein domains, remove all docking terminal predictions that
        aren't overlapping with the first or last 50 amino acids of the protein.
    """
    dockingdomains = {'NRPS-COM_Nterm', 'NRPS-COM_Cterm',
                      'PKS_Docking_Cterm', 'PKS_Docking_Nterm'}
    feature_by_id = record.get_cds_name_mapping()
    results = {}
    for gene_name in list(gene_domains):
        new = []
        gene_length = len(feature_by_id[gene_name].translation)
        for hit in gene_domains[gene_name]:
            if hit.hit_id in dockingdomains and \
                    not (gene_length - max(hit.query_start, hit.query_end) < 50
                         or min(hit.query_start, hit.query_end) < 50):
                continue
            new.append(hit)
        if new:
            results[gene_name] = new
    return results


def find_ab_motifs(fasta):
    # Analyse for abMotifs
    opts = ["-E", "0.25"]
    motif_file = path.get_full_path(__file__, "data", "abmotifs.hmm")
    abmotif_results = subprocessing.run_hmmscan(motif_file, fasta, opts)
    lengths = utils.get_hmm_lengths(motif_file)
    return refine_hmmscan_results(abmotif_results, lengths, neighbour_mode=True)


def find_domains(fasta, record):
    # Analyse for C/A/PCP/E/KS/AT/ATd/DH/KR/ER/ACP/TE/TD/COM/Docking/MT/CAL domains
    opts = ["--cut_tc"]
    nrpspks_file = path.get_full_path(__file__, "data", "nrpspksdomains.hmm")
    nrpspksdomain_results = subprocessing.run_hmmscan(nrpspks_file, fasta, opts)
    lengths = utils.get_hmm_lengths(nrpspks_file)
    domains = refine_hmmscan_results(nrpspksdomain_results, lengths, neighbour_mode=True)
    return filter_nonterminal_docking_domains(record, domains)


def find_ks_domains(fasta):
    # Analyse KS domains & PKS/NRPS protein domain composition to detect NRPS/PKS types
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
    def __init__(self, domain_names):
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

    def trans_is_greatest(self):
        """ Returns true if the trans_at count is strictly greater than others """
        return self.trans_at > max([self.modular, self.enediyne, self.iterative])

    def ene_is_greatest(self):
        """ Returns true if the enediyne count is strictly greater than others """
        return self.enediyne > max([self.modular, self.trans_at, self.iterative])

    def modular_is_greatest(self):
        """ Returns true if the modular count is strictly greater than others """
        return self.modular > max([self.enediyne, self.trans_at, self.iterative])

    def iterative_is_greatest(self):
        """ Returns true if the iterative count is strictly greater than others """
        return self.iterative > max([self.enediyne, self.trans_at, self.modular])


def classify_feature(domain_names) -> str:
    """ """
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


def generate_motif_features(record, feature, motifs):
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
        new_motif.notes.append("NRPS/PKS Motif: " + motif.hit_id + " (e-value: " + str(motif.evalue) + ", bit-score: " + str(motif.bitscore) + ")")  # TODO move to CDSMotif

        motif_features.append(new_motif)
    return motif_features
