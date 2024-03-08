# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import re
import logging
from collections import defaultdict
from typing import List, Optional, Union

from antismash.common.secmet import Record, CDSFeature
from antismash.common import (
    subprocessing,
    path,
    fasta
)

from antismash.detection.genefunctions import halogenases
from antismash.detection.genefunctions.halogenases import (
    FlavinDependents,
    TryptophanSubstrate,
    PyrroleSubstrate,
    PhenolicSubstrate,
    HalogenaseHmmResult,
    HalogenaseResult
)

FDH_SUBGROUPS = {"trp_5_FDH": TryptophanSubstrate,
                 "trp_6_7_FDH": TryptophanSubstrate,
                 "pyrrole_FDH": PyrroleSubstrate,
                 "cycline_orsellinic_FDH": PhenolicSubstrate,
                 "tyrosine-like_hpg_FDH": PhenolicSubstrate,}

def get_profile_list() -> List[str]:
    """get the subclasses for a substrate in the FDH family"""
    return [profile.name for profile in halogenases.SPECIFIC_FDH_PROFILES]

def check_conserved_motif(cds: CDSFeature, motif_positions: list[int],
                           hmm_result: HalogenaseHmmResult,
                           motif_pattern):

    categorized = ""

    siganture_residues = halogenases.search_signature_residues(cds.translation,
                                                               motif_positions,
                                                               hmm_result)
    if not siganture_residues:
        return categorized

    motif = re.search(motif_pattern, siganture_residues)
    if not motif:
        return categorized
    else:
        categorized = re.search(f"{motif_pattern}", siganture_residues)[0]
    return categorized

def run_halogenase_phmms(cluster_fasta: str, profiles: list) \
    -> dict[str, list[HalogenaseHmmResult]]:
    """ Check if protein sequences hit any pHMM

        Arguments:
            cluster_fasta: string of protein sequences in a fasta format

        Returns:
            if there is a hit, it returns the properties of that hit
    """
    halogenase_hmms_by_id: dict = defaultdict(list)
    for sig in profiles:
        sig.path = path.get_full_path(f'{sig.hmm_file}/{sig.name}')
        runresults = subprocessing.run_hmmsearch(sig.path, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    hit = HalogenaseHmmResult(hsp.hit_id, hsp.bitscore,
                                              hsp.query_id, sig.name, sig.path)
                    halogenase_hmms_by_id[hsp.hit_id].append(hit)
    return halogenase_hmms_by_id

def check_for_fdh(cds: CDSFeature, halogenase: HalogenaseResult, hit: HalogenaseHmmResult) -> None:
    """ Check if protein could be categorized as a Flavin-dependent enzyme

        Arguments:
            cds: gene/CDS and its properties
            halogenase: enzyme categorized as halogenase
            hit: details of the hit (e.g. bitscore, name of the profile, etc.)

        Returns:
            Modifies the HalogenasesResults by adding the label of the enzyme family,
            and possibly the signature and position of halogenation determined.
            If it doesn't meet the requirements, it doesn't make changes on the original instance.
    """
    for profile_name in get_profile_list():
        if profile_name == hit.query_id:
            signature_residues = FDH_SUBGROUPS[profile_name].get_residues(cds, hit)
            specific_signature_residues = signature_residues[hit.query_id]
            if not specific_signature_residues:
                return
            halogenase_category = FDH_SUBGROUPS[profile_name](hit.query_id,
                                                            specific_signature_residues)
            halogenase_category.update_match(halogenase, hit)

def check_for_halogenases(cds: CDSFeature, hmm_results: list[HalogenaseHmmResult]
                          ) -> Optional[HalogenaseResult]:
    """ Categorizes enzymes based on wether they hit the pHMM
        and have the required signature residues

        Arguments:
            cds: the query cds
            hmm_results: list of the halogenase hits and their properties
        Returns:
            categorized halogenase
    """
    if not hmm_results:
        logging.debug("Hmmsearch did not return any hit.")
        return None

    halogenase_match = HalogenaseResult(cds.get_name())
    for hit in hmm_results:
        check_for_fdh(cds, halogenase_match, hit)
    return halogenase_match

def specific_analysis(record: Record) -> Optional[list[HalogenaseResult]]:
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
    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    general_hmm_hits = run_halogenase_phmms(hmmsearch_fasta,
                                    halogenases.GENERAL_FDH_PROFILES)

    if general_hmm_hits:
        hmm_hits = run_halogenase_phmms(hmmsearch_fasta,
                                        halogenases.SPECIFIC_FDH_PROFILES)
        for protein, hits in hmm_hits.items():
            cds = record.get_cds_by_name(protein)
            potential_enzyme = check_for_halogenases(cds, hits)
            if potential_enzyme:
                potential_enzyme.family = "Flavin-dependent halogenases"
                potential_enzymes.append(potential_enzyme)
        for enzyme in potential_enzymes:
            enzyme.finalize_enzyme()
    if not potential_enzymes and general_hmm_hits:
        conserved_motifs = {}
        for protein, hits in general_hmm_hits.items():
            cds = record.get_cds_by_name(protein)
            for hit in hits:
                for motif, positions in FlavinDependents.GENERAL_FDH_MOTIFS.items():
                    conserved_motif = check_conserved_motif(cds, positions,
                                                hit, motif)
                    if conserved_motif: 
                        conserved_motifs[motif] = conserved_motif 
             
            potential_enzyme = FlavinDependents(protein, conserved_motifs)
            
            if potential_enzyme:
                potential_enzymes.append(potential_enzyme)

    print(potential_enzymes)
    return potential_enzymes
