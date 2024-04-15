# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import logging
import re
from collections import defaultdict
from typing import Union, Optional

from antismash.common.secmet import (
    CDSFeature,
    Record)

from antismash.common import (
    subprocessing,
    fasta,
    utils)

from antismash.detection.genefunctions.halogenases.halogenases import (
    FlavinDependentHalogenases,
    HalogenaseHmmResult
    )


from antismash.common.subprocessing import hmmscan
from antismash.detection.genefunctions.halogenases.flavin_dependent import substrates
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
    phenolic,
    pyrrolic
)

# pHMM ids and the submodules they belong to for use in calling submodule-specific functions
FDH_SUBGROUPS = {"trp_5_FDH":indolic,
                 "trp_6_7_FDH":indolic,
                 "cycline_orsellinic_FDH":phenolic,
                 "tyrosine-like_hpg_FDH":phenolic,
                 "pyrrole_FDH":pyrrolic}

def _get_substrate_specific_profiles() -> list:
    """ Collects the substrate-specific pHMM profiles from the substrate-specific submodules"""
    profiles = []
    submodules = list(set(FDH_SUBGROUPS.values()))
    for submodule in submodules:
        for profile in submodule.SPECIFIC_PROFILES:
            profiles.append(profile)
    return profiles


def retrieve_fdh_signature_residues(translation: str, hmm_result: HalogenaseHmmResult,
                                    signatures: Union[list[list[int]], list[int]], enzyme_substrates: list = None
                                    ) -> dict[str, Optional[str]]:
    """ Get signature residues for an enzyme from each pHMM

        Arguments:
            sequence: protein sequence
            hmm_result: instance of HmmResult class,
                        which contains information about the hit in a pHMM
            signatures: list of the positions that defines the signature residues in a pHMM

        Returns:
            signature residues which were retrieved from a certain pHMM
    """
    signature_residues: dict[str, Optional[str]] = {}
    if not enzyme_substrates:
        residue = search_residues(translation, signatures, hmm_result)
        signature_residues[hmm_result.query_id] = residue
    else:
        substrates_signatures = dict(zip(enzyme_substrates,signatures))
        for substrate, signature in substrates_signatures.items():
            signature_residues[substrate] = search_residues(translation, signature, hmm_result)
    return signature_residues

def search_residues(sequence: str, positions: Union[list[int], list[list[int]]],
                    hmm_result: HalogenaseHmmResult,
                    max_evalue: float = 0.1) -> Optional[str]:
    """ Get the signature residues from the pHMM for the searched protein sequence

        Arguments:
            sequence: protein sequence
            positions: list of position numbers in the pHMM
            hmm_result: properties of the halogenase pHMM hit
            max_evalue: maximum e-value

        Returns:
            residues that are present in the given positions
    """
    args = ["-E", str(max_evalue)]
    results = subprocessing.hmmpfam.run_hmmpfam2(hmm_result.profile,
                                                 f">query\n{sequence}", extra_args=args)
    if not (results and results[0].hsps):
        logging.debug("no hits for query %s")
        return None

    found = False
    hit = None
    for hit in results[0].hsps:
        if hit.hit_id == hmm_result.query_id:
            found = True
            break

    if not found:
        logging.debug(
            "no hits for the enzyme in %s", hmm_result.query_id)
        return None

    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    sites = utils.extract_by_reference_positions(query, profile,
                                                 [p - offset for p in positions if offset < p])
    return sites

def search_conserved_motif(cds: CDSFeature, motif_positions: list[int],
                           hmm_result: HalogenaseHmmResult,
                           motif_pattern: str) -> str:
    """ Looks for WxWxIP and Fx.Px.Sx.G conserved motifs, characteristic to FDHs

        Arguments:
            cds: gene/CDS and its properties
            motif_positions: positions of the conserved motifs in the pHMM
            hmm_result: details of the hit (e.g. bitscore, name of the profile, etc.)
            motif_pattern: pattern of the motif in regex

        Returns:
            if the conserved motifs are present, then it returns those,
            if they are not present it returns an empty string
    """

    categorized = ""
    signature_residues = search_residues(cds.translation,
                                         motif_positions,
                                         hmm_result)
    if not signature_residues:
        return categorized

    motif = re.search(motif_pattern, signature_residues)
    if not motif:
        return categorized
    categorized = (re.search(f"{motif_pattern}", signature_residues) or [])[0]
    return categorized

def run_halogenase_phmms(cluster_fasta: str, profiles: list
                         ) -> dict[str, list[HalogenaseHmmResult]]:
    """ Check if protein sequences hit any pHMM

        Arguments:
            cluster_fasta: string of protein sequences in a fasta format

        Returns:
            if there is a hit, it returns the properties of that hit
    """
    halogenase_hmms_by_id: dict = defaultdict(list)
    for sig in profiles:
        runresults = subprocessing.run_hmmsearch(sig.hmm_file, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    hit = HalogenaseHmmResult(hsp.hit_id, hsp.bitscore,
                                              hsp.query_id, sig.name, sig.path)
                    halogenase_hmms_by_id[hsp.hit_id].append(hit)
    return halogenase_hmms_by_id



def categorize_on_substrate_level(cds: CDSFeature, halogenase_match: FlavinDependentHalogenases,
                                  hmm_results: list[HalogenaseHmmResult]
                                  ) -> Optional[FlavinDependentHalogenases]:
    """ Check if protein could be categorized as a Flavin-dependent enzyme

        Arguments:
            cds: gene/CDS and its properties
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            Modifies the HalogenasesResults by adding the label of the enzyme family,
            and possibly the signature and position of halogenation determined.
            If it doesn't meet the requirements, it doesn't make changes on the original instance.
    """
    if not hmm_results:
        logging.debug("Hmmsearch did not return any hit.")
        return None

    for hit in hmm_results:
        if hit.query_id in FDH_SUBGROUPS:
            signature_residues = FDH_SUBGROUPS[hit.query_id].get_consensus_signature(cds, hit)
            specific_signature_residues = signature_residues[hit.query_id]
            if not specific_signature_residues:
                return None
            FDH_SUBGROUPS[hit.query_id].update_match(hit.query_id,
                                                    specific_signature_residues,
                                                    halogenase_match, hit)
    if not halogenase_match.potential_matches:
        return None

    return halogenase_match

def fdh_specific_analysis(record: Record) -> Union[list, list[FlavinDependentHalogenases]]:
    """ Categorization of enzyme, categorizes any halogenase in a cds in regions

        Arguments: record instance,
                   which holds information of the identified clusters

        Returns: list of HalogenasesResults instances representing halogenase enzymes,
                 if there is a clear best match for a given enzyme, then the information
                 about the position of halogenation, the confidence of the categorization,
                 and the characteristic residues is provided.
                 If the enzyme can be categorized into several groups, with the same confidence,
                 then the above mentioned informations are not defined,
                 and the information about the catogries is in the potential_enzymes attribute.
    """

    potential_enzymes = []
    enzymes_with_hits = []
    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    hits = hmmscan.run_hmmscan(substrates.ALL_FDH_PROFILES,
                               hmmsearch_fasta)
    for query_result in hits:
        if query_result.hits:
            found = re.search(f'>{query_result.id}\n.*\n',
                              hmmsearch_fasta)
            if not found:
                continue
            enzymes_with_hits.append(found.group())
    if enzymes_with_hits:
        hit_enzyme_fasta = "".join(enzymes_with_hits)
        general_hmm_hits = run_halogenase_phmms(hit_enzyme_fasta,
                                                substrates.GENERAL_FDH_PROFILES)
        if general_hmm_hits:
            specific_profiles = _get_substrate_specific_profiles()
            specific_hmm_hits = run_halogenase_phmms(hit_enzyme_fasta,
                                                     specific_profiles)
        for protein in general_hmm_hits.keys():
            cds = record.get_cds_by_name(protein)
            enzyme = FlavinDependentHalogenases(protein, cofactor="flavin", family="FDH")

            if protein in specific_hmm_hits.keys():
                enzyme = (categorize_on_substrate_level(cds, enzyme, specific_hmm_hits[protein]) or
                          enzyme)

            # if it's not in the specific hits or couldn't be categorized further
            elif not enzyme.potential_matches:
                conserved_motifs = {}
                for hit in general_hmm_hits[protein]:
                    for motif, positions in substrates.GENERAL_FDH_MOTIFS.items():
                        conserved_motif = search_conserved_motif(cds, positions,
                                                                 hit, motif)
                        if conserved_motif:
                            conserved_motifs[motif] = conserved_motif

                enzyme.consensus_residues = conserved_motifs

            potential_enzymes.append(enzyme)

        for enzyme in potential_enzymes:
            enzyme.finalize_enzyme()

    return potential_enzymes
