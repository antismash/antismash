import logging
import importlib
import pkgutil
import re
from collections import defaultdict
from typing import cast, Union, Optional, List
from pathlib import Path

from antismash.common.path import get_full_path
from antismash.custom_typing import AntismashModule
from antismash.common.secmet import (
    CDSFeature,
    Record)

from antismash.common import (
    subprocessing,
    fasta,
    utils)

from antismash.detection.genefunctions.halogenases.halogenases import (
    Match,
    TailoringEnzymes,
    HalogenaseHmmResult
    )

from antismash.detection.genefunctions.halogenases.flavin_dependent import substrates
from antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups import (
    indolic,
    phenolic,
    pyrrolic
)

FDH_SUBGROUPS = {"trp_5_FDH":indolic,
                 "trp_6_7_FDH":indolic,
                 "cycline_orsellinic_FDH":phenolic,
                 "tyrosine-like_hpg_FDH":phenolic,
                 "pyrrole_FDH":pyrrolic}

def _gather_fdh_substrate_modules() -> List[AntismashModule]:
    modules = []
    for module_data in pkgutil.walk_packages([get_full_path(Path(__file__).parents[0],
                                                            "flavin_dependent", "subgroups")]):
        module = importlib.import_module(f"antismash.detection.genefunctions.halogenases.flavin_dependent.subgroups.{module_data.name}")
        modules.append(cast(AntismashModule, module))
    return modules

_ANALYSIS_MODULES = _gather_fdh_substrate_modules()

def _get_analysis_modules() -> List[AntismashModule]:
    """ Return a list of default analysis modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return list(_ANALYSIS_MODULES)

def _get_substrate_specific_profiles():
    profiles = []
    submodules = _get_analysis_modules()
    for submodule in submodules:
        for profile in submodule.SPECIFIC_PROFILES:
            profiles.append(profile)

    return profiles

def check_for_match(name, residues, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult,
                    position: Union[int, List[int]], cutoffs: List[float], *,
                    check_residues: bool = True, sig_residues: Union[str, dict[str,str]] = "",
                    confidence: float = 1, targets: bool = False) -> bool:
    """in case the predefined signatures (self.signatures) is a dict,
        so there are several substrates, the searched signatures (sig_residues)
        also have to be a dict"""

    found = False
    if hit.query_id != name:
        return False
    cutoffs.sort(reverse=True)
    modifier = 1.
    for cutoff in cutoffs:
        if not isinstance(residues, dict):
            if hit.bitscore >= cutoff and (residues == sig_residues or not check_residues):
                halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", position,
                                                    confidence * modifier, residues,""))
                return True
            modifier = 0.5
        elif isinstance(residues, dict) and isinstance(sig_residues, dict):
            for subs, sig_res in residues.items():
                if hit.bitscore >= cutoff and (sig_res == sig_residues[subs] or
                                                not check_residues):
                    if not targets:
                        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", position,
                                                                confidence * modifier,
                                                                sig_res, subs))
                    else:
                        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", subs,
                                                                confidence * modifier,
                                                                sig_res))
                    return True
            modifier = 0.5

        elif isinstance(residues, dict) and not isinstance(sig_residues, dict):
            for subs, sig_res in residues.items():
                if hit.bitscore >= cutoff and (sig_res == sig_residues or
                                                not check_residues):
                    if not targets:
                        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", position,
                                                                confidence * modifier,
                                                                sig_res, subs))
                    else:
                        halogenase.add_potential_matches(Match(hit.query_id, "flavin", "FDH", subs,
                                                                confidence * modifier,
                                                                sig_res))
                    return True
            modifier = 0.5
    return found

def get_residues(translation: str, hmm_result: HalogenaseHmmResult,
                 signatures: Union[list[int], list[list[int]]], substrates: list = None)\
                 -> dict[str, Optional[str]]:
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

    if not substrates:
        residue = search_signature_residues(translation, signatures, hmm_result)
        signature_residues[hmm_result.query_id] = residue
    else:
        if len(signatures) == 1:
            substrates_signatures = dict(zip(substrates,(3*signatures)))
        else:
            substrates_signatures = dict(zip(substrates,signatures))
        for substrate, signature in substrates_signatures.items():
            signature_residues[substrate] = search_signature_residues(translation,
                                                                      signature, hmm_result)

    return signature_residues

def search_signature_residues(sequence: str, positions: list[int],
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

def check_conserved_motif(cds: CDSFeature, motif_positions: list[int],
                           hmm_result: HalogenaseHmmResult,
                           motif_pattern):

    categorized = ""

    siganture_residues = search_signature_residues(cds.translation,
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
        runresults = subprocessing.run_hmmsearch(sig.hmm_file, cluster_fasta)
        for runresult in runresults:
            for hsp in runresult.hsps:
                if hsp.bitscore > sig.cutoff:
                    hit = HalogenaseHmmResult(hsp.hit_id, hsp.bitscore,
                                              hsp.query_id, sig.name, sig.path)
                    halogenase_hmms_by_id[hsp.hit_id].append(hit)
    return halogenase_hmms_by_id



def check_for_fdh(cds: CDSFeature, halogenase: TailoringEnzymes, hit: HalogenaseHmmResult) -> None:
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
    signature_residues = FDH_SUBGROUPS[hit.query_id].get_consensus_signature(cds, hit)
    specific_signature_residues = signature_residues[hit.query_id]
    if not specific_signature_residues:
        return

    FDH_SUBGROUPS[hit.query_id].update_match(hit.query_id,
                                             specific_signature_residues,
                                             halogenase, hit)

def check_for_halogenases(cds: CDSFeature, hmm_results: list[HalogenaseHmmResult]
                          ) -> Optional[TailoringEnzymes]:
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

    halogenase_match = TailoringEnzymes(cds.get_name())
    for hit in hmm_results:
        check_for_fdh(cds, halogenase_match, hit)
    return halogenase_match

def fdh_specific_analysis(record: Record) -> Optional[list[TailoringEnzymes]]:
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
    profiles = _get_substrate_specific_profiles()

    potential_enzymes = []
    features = record.get_cds_features_within_regions()
    hmmsearch_fasta = fasta.get_fasta_from_features(features)

    general_hmm_hits = run_halogenase_phmms(hmmsearch_fasta,
                                    substrates.GENERAL_FDH_PROFILES)

    if general_hmm_hits:
        hmm_hits = run_halogenase_phmms(hmmsearch_fasta,
                                        profiles)
        for protein, hits in hmm_hits.items():
            cds = record.get_cds_by_name(protein)
            potential_enzyme = check_for_halogenases(cds, hits)
            if potential_enzyme:
                potential_enzymes.append(potential_enzyme)
        for enzyme in potential_enzymes:
            enzyme.finalize_enzyme()
    if not potential_enzymes and general_hmm_hits:
        conserved_motifs = {}
        for protein, hits in general_hmm_hits.items():
            cds = record.get_cds_by_name(protein)
            for hit in hits:
                for motif, positions in substrates.GENERAL_FDH_MOTIFS.items():
                    conserved_motif = check_conserved_motif(cds, positions,
                                                hit, motif)
                    if conserved_motif:
                        conserved_motifs[motif] = conserved_motif

            potential_enzyme = TailoringEnzymes(protein, cofactor="flavin", family="FDH", consensus_residues=conserved_motifs)

            if potential_enzyme:
                potential_enzymes.append(potential_enzyme)
    print(potential_enzymes)
    return potential_enzymes
