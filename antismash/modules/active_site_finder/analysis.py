# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Analysis methods for a variety of active site types.

    Each method should include a reference to its data, the source of the data,
    and a reference to a paper from which the method is based.
"""

from collections import defaultdict
from typing import Dict, List

from antismash.common import secmet

from .common import ActiveSiteAnalysis

# these are type aliases, not constants
SinglePairing = tuple[secmet.features.Domain, str]
MultiplePairing = tuple[secmet.features.Domain, list[str]]


def run_all_analyses(record: secmet.Record) -> List[MultiplePairing]:
    """ Runs all AFS analyses at once and returns their aggregated results """
    hits_by_feature: Dict[secmet.features.Domain, List[str]] = defaultdict(list)
    for analysis in [asp_ks, asp_ks_c, asp_at, acp_type, asp_acp, asp_pksi_dh, asp_pksi_kr,
                     asp_thioesterase, pksi_er_stereo, pksi_kr_stereo, pksi_at_spec, asp_p450_oxy]:
        for feature, hit in analysis(record):
            hits_by_feature[feature].append(hit)
    return list(hits_by_feature.items())


def pksi_kr_stereo(record: secmet.Record) -> List[SinglePairing]:
    """ Pediction of PKSI KR specificities according to Reid et al., Biochemistry 2003, 42, 72-79

        Data: PKSI-KR.hmm2
        Data source: CLUSEAN
        Reference: R. Reid, M. Piagentini, E. Rodriguez, G. Ashley, N. Viswanathan,
                   J. Carney, D. V. Santi, C. R. Hutchinson, and R. McDaniel. 2003.
                   A model of structure and catalysis for ketoreductase domains
                   in modular polyketide synthases. Biochemistry 42:72-79.
    """
    # method_type = "prediction"
    positions = [149, 162, 166]
    expected = ["S", "Y", "N"]
    target_domain = "PKS_KR"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")

    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-KR.hmm2", positions, expected)

    label = "KR domain putatively catalyzing {}-configuration product formation"
    results = []
    for alignment in analyser.get_alignments():
        value = alignment.get_signature([102])
        if value == "D":
            results.append((alignment.domain, label.format("D")))
        else:
            results.append((alignment.domain, label.format("L")))
    return results


def asp_ks(record: secmet.Record) -> List[SinglePairing]:
    """
        KS active site cysteine according to Huang et al. Microbiology 147 (2001), 631-642

        Data: PKSI-KS_N.hmm
        Data source: PFAM
        Reference: Guozhong Huang, Lianhui Zhang and Robert G. Birch.
                   A multifunctional polyketide–peptide synthetase essential for
                   albicidin biosynthesis in Xanthomonas albilineans.
                   Microbiology 147 (2001), 631-642
    """
    # method_type = "active_site_toCorrect"
    positions = [176, 186, 187, 188]
    expected = ['G', 'S', 'S', 'S']
    emissions = [0.99, 0.9, 0.81, 0.81]
    target_domain = "PKS_KS"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")

    results = []

    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-KS_N.hmm2",
                                  positions, expected, emissions=emissions)

    label = "found active site cysteine: {}, scaffold matched {}: {}"
    for alignment in analyser.get_alignments():
        value = alignment.get_signature([185])
        new_label = label.format(value == "C", "".join(expected), analyser.scaffold_matches(alignment))
        results.append((alignment.domain, new_label))
    return results


def asp_ks_c(record: secmet.Record) -> List[SinglePairing]:
    """ KS active site Histidine according to Huang et al. Microbiology 147 (2001), 631-642

        Data: PKSI-KS_C.hmm
        Data source: PFAM
        Reference: Guozhong Huang, Lianhui Zhang and Robert G. Birch.
                   A multifunctional polyketide–peptide synthetase essential for
                   albicidin biosynthesis in Xanthomonas albilineans.
                   Microbiology 147 (2001), 631-642
    """
    # method_type = "active_site_toCorrect"
    positions = [47, 50, 51, 52, 53, 56, 57, 60, 106, 110, 117, 123]
    expected = ['E', 'G', 'T', 'G', 'T', 'G', 'D', 'E', 'K', 'G', 'G', 'K']
    emissions = [0.96, 0.92, 0.95, 0.94, 0.95, 0.96, 0.99, 0.99, 1., 0.99, 0.96, 0.93]
    target_domain = "PKS_KS"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")

    results = []

    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-KS_C.hmm2",
                                  positions, expected, emissions=emissions)

    for alignment in analyser.get_alignments():
        value = alignment.get_signature([49, 111])  # emissions [0.95, 1.]
        results.append((alignment.domain, f"found active site histidines: {value == 'HH'}"))

    return results


def asp_at(record: secmet.Record) -> List[SinglePairing]:
    """
        Data: PKSI-AT.hmm2
        Data source: CLUSEAN
        Reference: unknown
    """
    # method_type = "active_site"
    positions = [93, 94, 97, 98]
    expected = ['G', 'H', 'G', 'E']
    emissions = [1., 1., 1., 0.94]
    target_domain = "PKS_AT"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-AT.hmm2",
                                  positions, expected, emissions=emissions)

    results = []

    label = "found active site serine: %s, scaffold matched GHGE: %s"
    for alignment in analyser.get_alignments():
        value = alignment.get_signature([95])  # emission 1.
        results.append((alignment.domain, label % (value == "S", analyser.scaffold_matches(alignment))))

    return results


def acp_type(record: secmet.Record) -> List[SinglePairing]:
    """ Beta-branching motif of ACPs according to Haines et al., Nature Chemical Biology 9:685-692

        Data: PKSI-ACP.hmm2
        Data source: CLUSEAN
        Reference: Haines et al.
                   A conserved motif flags acyl carrier proteins for beta-branching
                   in polyketide synthesis. Nature Chemical Biology 9:685-692
    """
    # method_type = "prediction"
    positions = [28, 30, 31]
    expected = ['G', 'D', 'S']
    emissions = [0.97, 0.9, 0.98]
    target_domain = "ACP"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-ACP.hmm2",
                                  positions, expected, emissions=emissions)
    results = []

    for alignment in analyser.get_alignments():
        value = alignment.get_signature([37])  # emission 5.68e-5 for W
        if value == "W":
            label = "beta-branching ACP"
        else:
            label = "non-beta-branching ACP"
        results.append((alignment.domain, label))

    return results


def asp_acp(record: secmet.Record) -> List[SinglePairing]:
    """ ACP active site serine according to Huang et al. Microbiology 147 (2001), 631-642

        Data: PKSI-ACP.hmm2
        Data source: CLUSEAN
        Reference: Guozhong Huang, Lianhui Zhang and Robert G. Birch.
                   A albicidin biosynthesis in Xanthomonas albilineans.
                   Microbiology 147 (2001), 631-642
    """
    # method_type = "active_site"
    positions = [28, 30]
    expected = ['G', 'D']
    emissions = [0.97, 0.9]
    target_domain = "ACP"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-ACP.hmm2",
                                  positions, expected, emissions=emissions)

    results = []

    for alignment in analyser.get_alignments():
        value = alignment.get_signature([31])  # emission 0.98
        results.append((alignment.domain, f"found active site serine: {value == 'S'}"))

    return results


def asp_pksi_dh(record: secmet.Record) -> List[SinglePairing]:
    """
        Data: PKSI-DH.hmm2
        Data source: CLUSEAN
        Reference: unknown
    """
    # method_type = "active_site"
    positions = [4, 5, 6, 30, 32, 34]
    expected = ['L', 'L', 'G', 'P', 'L', 'D']
    emissions = [0.9, 0.78, 0.86, 0.7, 0.85, 0.78]
    target_domain = "PKS_DH"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-DH.hmm2",
                                  positions, expected, emissions=emissions)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "catalytic triad H,G,P inconclusive"))
            continue
        value = alignment.get_signature([5, 39, 44])
        found = value == "HGP"  # emissions 1., 0.62, 097
        results.append((alignment.domain, f"catalytic triad H,G,P found: {found}"))

    return results


def asp_pksi_kr(record: secmet.Record) -> List[SinglePairing]:
    """ Active site prediction according to Reid et al. Biochemistry 42:72-79 (2003)

        Data: PKSI-KR.hmm2
        Data source: CLUSEAN
        Reference: R. Reid, M. Piagentini, E. Rodriguez, G. Ashley, N. Viswanathan,
                   J. Carney, D. V. Santi, C. R. Hutchinson, and R. McDaniel. 2003.
                   A model of structure and catalysis for ketoreductase domains
                   in modular polyketide synthases. Biochemistry 42:72-79.
    """
    # method_type = "active_site"
    positions = [33, 65, 92, 97, 120, 148]
    expected = ['R', 'D', 'G', 'A', 'K', 'S']
    emissions = [0.93, 0.88, 0.93, 0.9, 0.99, 0.87]
    target_domain = "PKS_KR"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-KR.hmm2",
                                  positions, expected, emissions=emissions)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "catalytic triad S,Y,N inconclusive"))
            continue
        value = alignment.get_signature([149, 162, 166])
        found = value == "SYN"
        results.append((alignment.domain, f"catalytic triad S,Y,N found: {found}"))

    return results


def asp_thioesterase(record: secmet.Record) -> List[SinglePairing]:
    """ database: Thioesterase.hmm2
        database_source: CLUSEAN
        reference: unknown
    """
    # method_type = "active_site"
    positions = [73, 79, 83, 87, 108]
    expected = ['G', 'G', 'G', 'A', 'D']
    emissions = [0.93, 1., 0.99, 0.9, 0.9]
    target_domain = "Thioesterase"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "Thioesterase.hmm2",
                                  positions, expected, emissions=emissions)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "active site serine inconclusive"))
            continue
        value = alignment.get_signature([81])
        found = value == "S"
        results.append((alignment.domain, f"active site serine present: {found}"))

    return results


def pksi_er_stereo(record: secmet.Record) -> List[SinglePairing]:
    """ Pediction of PKSI ER specificities according to Kwan et all. ChemBiol 15:1231-1240, 2008

        Data: PKSI-ER.hmm2
        Data source: CLUSEAN
        Reference: Kwan, D. H., Y. Sun, F. Schulz, H. Hong, B. Popovic, J. C.,
                   Sim-Stark, S. F. Haydock, and P. F. Leadlay. 2008.
                   Prediction and  Manipulation of the Stereochemistry of
                   Enoylreduction in Modular Polyketide Synthases. Chem.Biol 15:1231-1240.
    """
    # method_type = "prediction"
    positions = [31, 135, 137, 144, 146, 227]
    expected = ['D', 'L', 'H', 'G', 'A', 'G']
    target_domain = "PKS_ER"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-ER.hmm2", positions, expected)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "ER configuration inconclusive"))
            continue
        value = alignment.get_signature([39])
        label = "ER domain putatively catalyzing {}-configuration product formation"
        if value == "Y":
            results.append((alignment.domain, label.format("2S")))
        else:
            results.append((alignment.domain, label.format("2R")))

    return results


def pksi_at_spec(record: secmet.Record) -> List[SinglePairing]:
    """ Pediction of PKSI AT specificities according to DelVecchio et al., 1999, J Ind. Microbiol Biotechnol 30:489-494.

        Data: PKSI-AT.hmm2
        Data source: CLUSEAN
        Reference: Del Vecchio, F., H. Petkovic, S. G. Kendrew, L. Low, B. Wilkinson,
                   R. Lill, J. Cortes, B. A. Rudd, J. Staunton, and P. F. Leadlay. 2003.
                   Active-site residue, domain and module swaps in modular
                   polyketide synthases. J Ind. Microbiol Biotechnol 30:489-494.
    """
    # method_type = "prediction"
    positions = [93, 94, 95, 120, 196, 198, 199, 227, 244, 245]
    expected = ['G', 'H', 'S', 'R', 'A', 'H', 'S', 'S', 'Y', 'W']
    target_domain = "PKS_AT"
    candidates = record.get_antismash_domains_by_tool("nrps_pks_domains")
    analyser = ActiveSiteAnalysis(target_domain, candidates, "PKSI-AT.hmm2", positions, expected)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "(Methyl)Malonyl-CoA specificity inconclusive"))
            continue
        value = alignment.get_signature([195, 197])
        if value == "HF":
            results.append((alignment.domain, "Malonyl-CoA specific"))
        elif value == "YS":
            results.append((alignment.domain, "Methylmalonyl-CoA specific"))
        else:
            results.append((alignment.domain, "Neither malonyl-CoA or methylmalonyl-CoA specific"))

    return results


def asp_p450_oxy(record: secmet.Record) -> List[SinglePairing]:
    """ Prediction of cytochrome P450 active site cysteine

        database: p450.hmm2
        database_source: PFAM31
        reference: Del Vecchio, F., H. Petkovic, S. G. Kendrew, L. Low, B. Wilkinson,
                   R. Lill, J. Cortes, B. A. Rudd, J. Staunton, and P. F. Leadlay. 2003
                   Active-site residue, domain and module swaps in modular
                   polyketide synthases. J Ind. Microbiol Biotechnol 30:489-494.
    """
    # method_type = "active_site"
    positions = [327, 330, 400, 403, 409]
    expected = ['E', 'R', 'F', 'G', 'G']
    target_domain = "p450"
    candidates = record.get_pfam_domains()
    analyser = ActiveSiteAnalysis(target_domain, candidates, "p450.hmm2", positions, expected)

    results = []

    for alignment in analyser.get_alignments():
        if not analyser.scaffold_matches(alignment):
            results.append((alignment.domain, "active site cysteine inconclusive"))
            continue
        value = alignment.get_signature([407])
        found = value == "C"
        results.append((alignment.domain, f"active site cysteine present: {found}"))
    return results
