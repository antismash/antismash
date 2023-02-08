# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the type II PKS module
    The analysis determines starter unit, chain length, molecular weight of the product and
    product class. It does so by identifying cyclases, ketoreductases, oxygenases,
    methylases and glycosyltransferases, etc.
    If applicable, the regioselectivity and function of identified auxillary and tailoring
    protiens will be predicted.

"""

import json
from typing import Dict, Iterable, List, Tuple, Set
from collections import defaultdict

from antismash.common import fasta, path, subprocessing
from antismash.common.hmmscan_refinement import refine_hmmscan_results, HMMResult
from antismash.common.utils import get_hmm_lengths, get_fasta_lengths
from antismash.common.secmet import Protocluster, CDSFeature, Record

from .results import ProtoclusterPrediction, CDSPrediction, Prediction

with open(path.get_full_path(__file__, "data", "weights.json"), "r", encoding="utf-8") as handle:
    WEIGHTS = json.load(handle)

with open(path.get_full_path(__file__, "data", "classification.json"), "r", encoding="utf-8") as handle:
    CLASSIFICATIONS = json.load(handle)


def run_t2pks_hmmscan(cds_features: Iterable[CDSFeature]) -> Dict[str, List[HMMResult]]:
    """ Runs hmmscan for type II PKS proteins on the given CDSFeatures

        Arguments:
            cluster: Protocluster on which to run the type II PKS hmmscan

        Returns:
            a dictionary of key: cds and value: list of HMMResults, for hmmscan results of the cluster
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "t2pks.hmm")
    hmm_results = subprocessing.run_hmmscan(hmm_file, cluster_fasta, opts=['--cut_tc'])
    hmm_lengths = get_hmm_lengths(hmm_file)
    return refine_hmmscan_results(hmm_results, hmm_lengths)


def run_starter_unit_blastp(cds_hmm_hits: Dict[CDSFeature, List[HMMResult]]
                            ) -> Dict[str, List[HMMResult]]:
    """ Runs blastp on starter unit coding sequences in given cluster

        Arguments:
            cds_hmm_hits: HMMResults by cds from type II PKS hmmscan

        Returns:
            a dictionary mapping CDS name to a list of HMMresults
    """
    blastp_results = []
    blastp_fasta_files = set()
    for cds, hmm_hits in cds_hmm_hits.items():
        query_sequence = fasta.get_fasta_from_features([cds])
        for hit in hmm_hits:
            if hit.hit_id not in ['KSIII', 'AT', 'AMID', 'LIG']:
                continue
            blast_database = path.get_full_path(__file__, 'data', hit.hit_id)
            blastp_results.extend(subprocessing.run_blastp(blast_database, query_sequence))
            blastp_fasta_files.add(path.get_full_path(__file__, "data", f"{hit.hit_id}.fasta"))

    if not blastp_results:
        return {}

    fasta_lengths = {}
    for fasta_file in blastp_fasta_files:
        fasta_lengths.update(get_fasta_lengths(fasta_file))

    results = refine_hmmscan_results(blastp_results, fasta_lengths)
    for hits in results.values():
        for i, hit in enumerate(hits):
            if not hit.hit_id.endswith("-CoA"):
                hits[i] = HMMResult(hit.hit_id + "-CoA", hit.query_start, hit.query_end, hit.evalue, hit.bitscore)
    return results


def get_cds_predictions_by_protein_type(cds_predictions: Dict[str, List[CDSPrediction]],
                                        protein_types: List[str]
                                        ) -> Dict[str, List[CDSPrediction]]:
    """ Generates a new dictionary of CDS name to CDSPredictions, removing all
        CDS entries that do not have at least one prediction with a desired
        protein type. All predictions with protein types not in the given list
        will be removed.

        Arguments:
            cds_predictions: a dictionary mapping CDS name to a
                 list of CDSPredictions for that CDS
            protein_types: a list of desired protein types

        Returns:
            a dictionary mapping CDS name to a
                 list of CDSPredictions for that CDS
    """
    ret = {}
    desired = set(protein_types)
    for cds, predictions in cds_predictions.items():
        applicable = []
        for pred in predictions:
            if pred.ptype in desired:
                applicable.append(pred)
        if applicable:
            ret[cds] = applicable
    return ret


def get_predictions_by_protein_type(protein_predictions: Dict[str, List[CDSPrediction]],
                                    protein_types: List[str]) -> List[CDSPrediction]:
    """ Generates a new mapping of protein type to list of CDSPredictions, limiting
        predictions to only those with a desired protein type.

        Arguments:
            cds_predictions: a dictionary mapping CDS name to a
                 list of CDSPredictions for that CDS
            protein_types: a list of desired protein types

        Returns:
            a dictionary mapping protein type to a
                 list of CDSPredictions with that protein type
    """
    results = []
    for protein in protein_types:
        results.extend(protein_predictions[protein])
    return results


def sum_predictions(predictions: List[CDSPrediction]) -> List[Prediction]:
    """ Sum CDS prediction scores by proposed function

        Arguments:
            predictions: a list of CDSPredictions to be summed

        Returns:
            a list of (function, score) tuples sorted by highest score
    """
    if not predictions:
        return []
    summed_preds: Dict[str, Prediction] = {}
    for pred in predictions:
        if pred.pfunc is None:
            continue
        keys = pred.pfunc.split('_')[-1].strip().split('/')
        for key in keys:
            score = pred.bitscore / len(keys)
            evalue = pred.evalue

            if key not in summed_preds:
                summed_pred = Prediction(key, 0., evalue)
                summed_preds[key] = summed_pred

            summed_pred = summed_preds[key]

            summed_pred.score += score
            if evalue < summed_pred.evalue:
                summed_pred.evalue = evalue

    return sorted(summed_preds.values(), key=lambda pred: pred.score, reverse=True)


def predict_product_class(cds_predictions: Dict[str, List[CDSPrediction]]) -> Set[str]:
    """ Predicts potential product classes for a cluster containing the given
        CDS predictions. Any returned product class will be a potential classification
        for all given CDS predictions.

        Only predictions of protein types CLF and CYC are used for determining
        product class.

        Arguments:
            cds_predictions: a dictionary mapping CDS name to
                 a list of CDSPredictions for that CDS
    """
    # TODO: include images for these product classes, either from Rasmus' Fig 1.2
    # or his ref 13: Fritzsche, Ishida, Hertweck, 2008
    # 'Orchestration of Discoid Polyketide Cyclization in the Resistomycin Pathway'

    cds_predictions = get_cds_predictions_by_protein_type(cds_predictions, ['CLF', 'CYC'])

    potential_classes = set()
    for i, predictions in enumerate(cds_predictions.values()):
        for pred in predictions:
            if pred.pfunc is None:
                continue
            new = set(CLASSIFICATIONS[pred.ptype][pred.pfunc])
            # start the set with something, since we're about to do intersections
            if i == 0:
                potential_classes.update(new)
            else:
                potential_classes = potential_classes.intersection(new)

    return potential_classes


def predict_molecular_weight(preds_by_protein: Dict[str, List[CDSPrediction]],
                             starter_names: List[str],
                             elongation_names: List[str]) -> Dict[str, float]:
    """ Predicts the molecular weights of a cluster, one for each combination of
        starter unit and elongation.

        Arguments:
            preds_by_protein: a mapping of protein type to a list of
                CDSPredictions with that protein type
            starter_names: a list of names of any starters found
            elongation_names: a list of lengths as strings, with multiple options
                joined by the pipe (|) character, e.g. "7|8"

        Returns:
            a dictionary mapping starter/elongation combination to the weight of
                that combination
    """

    tailoring_mw = 0.
    for ptype, predictions in preds_by_protein.items():
        if ptype not in WEIGHTS:
            continue
        if ptype == "HAL":
            tailoring_mw += WEIGHTS["HAL"]['cl'] * len(predictions)
        elif ptype == "CYC":
            cyc_weight = WEIGHTS["CYC"]
            for pred in predictions:
                if pred.pfunc and '/' in pred.pfunc:
                    tailoring_mw += 2 * cyc_weight
                else:
                    tailoring_mw += cyc_weight
        else:
            tailoring_mw += WEIGHTS[ptype] * len(predictions)

    mws = {}
    for unit in starter_names:
        starter_weight = WEIGHTS['starter_unit'][unit]
        for elongations in elongation_names:
            for elongation in elongations.split('|'):
                combo = f"{unit}_{elongation}"
                mws[combo] = (starter_weight
                              + int(elongation) * WEIGHTS['malonyl']
                              + tailoring_mw)
    return mws


def make_cds_predictions(cds_hmm_hits: Dict[str, List[HMMResult]],
                         cds_blastp_hits: Dict[str, List[HMMResult]]) -> Tuple[Dict[str, List[CDSPrediction]],
                                                                               Dict[str, List[CDSPrediction]]]:
    """ Converts HMMer and BLASTp hits into a uniform collection of CDSPredictions

        Arguments:
            cds_hmm_hits: a mapping of CDS name to a list of HMMResults found by hmmscan
            cds_blastp_hits: a mapping of CDS name to a list of HMMResults found by blastp

        Returns:
            a tuple of
                a dictionary mapping CDS name to a list of CDSPredictions for that CDS
                and a dictionary mapping protein type to a list of CDSPredictions with
                    that protein type
    """

    preds_by_cds: Dict[str, List[CDSPrediction]] = defaultdict(list)
    preds_by_protein: Dict[str, List[CDSPrediction]] = defaultdict(list)

    for cds_name, hmm_hits in cds_hmm_hits.items():
        # combine blast and hmmscan results
        hmm_hits.extend(cds_blastp_hits.get(cds_name, []))

        for hit in hmm_hits:
            info = hit.hit_id.split('_', 1)
            protein_type = info[0]
            protein_function = info[1] if len(info) == 2 else None
            prediction = CDSPrediction(protein_type, protein_function, hit.bitscore, hit.evalue)
            preds_by_cds[cds_name].append(prediction)
            preds_by_protein[protein_type].append(prediction)

    return preds_by_cds, preds_by_protein


def analyse_cluster(cluster: Protocluster, record: Record) -> ProtoclusterPrediction:
    """ Analyse a type II PKS cluster

        Arguments:
            cluster: the Protocluster to analyse
            record: the Record the given cluster belongs to

        Returns:
            a single ProtoclusterPrediction instance with analysis results
    """
    assert cluster.product == "T2PKS"
    hmm_results_by_name = run_t2pks_hmmscan(cluster.cds_children)

    hmm_results_by_cds = {record.get_cds_by_name(name): hits for name, hits in hmm_results_by_name.items()}
    t2pks_blastp_results = run_starter_unit_blastp(hmm_results_by_cds)
    preds_by_cds, preds_by_protein = make_cds_predictions(hmm_results_by_name, t2pks_blastp_results)

    starter_units = sum_predictions(get_predictions_by_protein_type(preds_by_protein, ['KSIII', 'AT', 'AMID', 'LIG']))
    if not starter_units:
        starter_units = [Prediction('acetyl-CoA', 0., 0.)]
    malonyl_elongations = sum_predictions(preds_by_protein['CLF'])
    product_classes = predict_product_class(preds_by_cds)
    weights = {}
    if malonyl_elongations:
        starter_names = [starter.name for starter in starter_units]
        elong_names = [elong.name for elong in malonyl_elongations]
        weights.update(predict_molecular_weight(preds_by_protein, starter_names, elong_names))

    return ProtoclusterPrediction(preds_by_cds, starter_units, malonyl_elongations,
                                  product_classes, weights, cluster.location.start, cluster.location.end)
