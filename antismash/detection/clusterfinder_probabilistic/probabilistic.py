# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The probabilistic/statistical section of clusterfinder.

    The data values used here are pulled from the original stand-alone ClusterFinder.
"""

import logging
import pickle
from typing import Dict, List, Tuple

import numpy as np

from antismash.common import path
from antismash.common.secmet import Record
from antismash.common.secmet.features import FeatureLocation
from antismash.config import ConfigType


# these pickled and unpickled numpy arrays are direct from the ClusterFinder source:
# https://github.com/petercim/ClusterFinder/tree/master/freqdata
with open(path.get_full_path(__file__, 'data', 'NewTS_all_B_reduced_6filter.pkl'), 'rb') as handle:
    EMISSION_PROBABILITIES = pickle.load(handle, encoding='latin1')
with open(path.get_full_path(__file__, 'data', 'NewTS_all_B_index.pkl'), 'rb') as handle:
    PFAM_INDICES = pickle.load(handle, encoding='latin1')
TRANSITION_PROBABILITIES = np.array([[0.94389093, 0.05610907],
                                     [0.01629607, 0.98370393]])

BIOSYNTHETIC_PFAMS = {"PF00109", "PF02801", "PF08659", "PF00378", "PF08541", "PF08545", "PF02803",
                      "PF00108", "PF02706", "PF03364", "PF08990", "PF00501", "PF00668", "PF08415",
                      "PF00975", "PF03061", "PF00432", "PF00494", "PF03936", "PF01397", "PF00432",
                      "PF04275", "PF00348", "PF02401", "PF04551", "PF00368", "PF00534", "PF00535",
                      "PF02922", "PF01041", "PF00128", "PF00908", "PF02719", "PF04321", "PF01943",
                      "PF02806", "PF02350", "PF02397", "PF04932", "PF01075", "PF00953", "PF01050",
                      "PF03033", "PF01501", "PF05159", "PF04101", "PF02563", "PF08437", "PF02585",
                      "PF01721", "PF02052", "PF02674", "PF03515", "PF04369", "PF08109", "PF08129",
                      "PF09221", "PF09683", "PF10439", "PF11420", "PF11632", "PF11758", "PF12173",
                      "PF04738", "PF04737", "PF04604", "PF05147", "PF08109", "PF08129", "PF08130",
                      "PF00155", "PF00202", "PF00702", "PF06339", "PF04183", "PF10331", "PF03756",
                      "PF00106", "PF01370", "PF00107", "PF08240", "PF00441", "PF02770", "PF02771",
                      "PF08028", "PF01408", "PF02894", "PF00984", "PF00725", "PF03720", "PF03721",
                      "PF07993", "PF02737", "PF00903", "PF00037", "PF04055", "PF00171", "PF00067",
                      "PF01266", "PF01118", "PF02668", "PF00248", "PF01494", "PF01593", "PF03992",
                      "PF00355", "PF01243", "PF00384", "PF01488", "PF00857", "PF04879", "PF08241",
                      "PF08242", "PF00698", "PF00483", "PF00561", "PF00583", "PF01636", "PF01039",
                      "PF00288", "PF00289", "PF02786", "PF01757", "PF02785", "PF02409", "PF01553",
                      "PF02348", "PF00891", "PF01596", "PF04820", "PF02522", "PF08484", "PF08421"}

REPEATS = {'PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130', 'PF01839', 'PF01816',
           'PF07720', 'PF00400', 'PF05594', 'PF07661', 'PF02985', 'PF06049', 'PF08238', 'PF06696',
           'PF00353', 'PF02412', 'PF00023', 'PF02071', 'PF03991', 'PF01469', 'PF07676', 'PF00514',
           'PF00904', 'PF07634', 'PF02370', 'PF03335', 'PF01851', 'PF04728', 'PF06715', 'PF03373',
           'PF04680', 'PF00805', 'PF04508', 'PF07918', 'PF01535', 'PF01011', 'PF05017', 'PF06671',
           'PF00818', 'PF03406', 'PF00399', 'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906',
           'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981', 'PF02095', 'PF00414',
           'PF00560', 'PF05001', 'PF02162', 'PF01473', 'PF05465', 'PF02493', 'PF03578', 'PF08043',
           'PF06392', 'PF07142', 'PF08309', 'PF02184'}


class ClusterFinderHit:
    "A putative cluster identified by ClusterFinder"
    def __init__(self, positions: Tuple[int, int], probability: float) -> None:
        self.location = FeatureLocation(positions[0], positions[1])
        self.probability = probability

    def __str__(self) -> str:
        return "CFHit({}, {})".format(self.location, self.probability)


def find_probabilistic_clusters(record: Record, options: ConfigType) -> List[ClusterFinderHit]:
    """ Find clusters based on ClusterFinder probabilities of PFAM features """
    pfam_features = record.get_pfam_domains()
    cf_clusters = []
    state = "seed"
    cluster_position = (0, 0)
    pfam_ids = []  # type: List[str]
    loop_index = 1
    probabilities = [0.]
    for feature in pfam_features:
        feature_position = (int(feature.location.start), int(feature.location.end))
        if feature.probability is None:
            loop_index += 1
            continue
        if feature.probability >= 0.3:
            if state == "seed":
                state = "extend"
                probabilities = [feature.probability]
                cluster_position = (min(feature_position), max(feature_position))
                pfam_ids.clear()
            else:
                probabilities.append(feature.probability)
                if max(feature_position) > cluster_position[1]:
                    cluster_position = (cluster_position[0], max(feature_position))
            pfam_ids.append(feature.identifier)
        else:
            if state == "extend":
                state = "seed"
                cluster_position, cds_count = find_nr_cds(cluster_position, record)
                if is_good_cluster_hit(cds_count, probabilities, pfam_ids, options):
                    cf_clusters.append(ClusterFinderHit(cluster_position, np.mean(probabilities)))
                cluster_position = (0, 0)
                pfam_ids = []
        if loop_index == len(pfam_features):
            if cluster_position != (0, 0):
                cluster_position, cds_count = find_nr_cds(cluster_position, record)
                if is_good_cluster_hit(cds_count, probabilities, pfam_ids, options):
                    cf_clusters.append(ClusterFinderHit(cluster_position, np.mean(probabilities)))
            cluster_position = (0, 0)
            pfam_ids = []
        loop_index += 1
    logging.debug("ClusterFinder detected %d probabilistic clusters", len(cf_clusters))
    return cf_clusters


def is_good_cluster_hit(num_cds: int, probabilities: List[float], pfam_ids: List[str], options: ConfigType) -> bool:
    "Check if the current cluster is a good hit"
    if num_cds < options.cf_min_cds_features:
        return False
    if np.mean(probabilities) < options.cf_threshold:
        return False
    unique_pfams = len(BIOSYNTHETIC_PFAMS.intersection(set(pfam_ids)))
    return bool(unique_pfams >= options.cf_min_pfams)


def find_nr_cds(cluster_position: Tuple[int, int], record: Record) -> Tuple[Tuple[int, int], int]:
    """ Find the number of CDSs in candidate cluster and adjust the cluster starts
        and ends to match the CDS starts and ends """
    area = FeatureLocation(cluster_position[0], cluster_position[1])
    cds_features = record.get_cds_features_within_location(area, with_overlapping=True)

    if not cds_features:
        return cluster_position, 0

    startlocations = [int(cds.location.start) for cds in cds_features]
    endlocations = [int(cds.location.end) for cds in cds_features]
    # avoid getting the complete genome as cluster if one CDS
    # starts at end and finishes at start of genome
    if not (0 in startlocations and len(record.seq) in endlocations):
        cluster_position = (min(startlocations), max(endlocations))
    return cluster_position, len(cds_features)


def get_pfam_probabilities(observations: List[str]) -> List[float]:
    """ Converts a list of PFAM IDs into a matching list of probabilites.

    Arguments:
        observations: list of PFAM IDs ordered by position in genome

    Returns:
        a list of probabilities as floats, one float for each observation
    """
    return forward_backward(TRANSITION_PROBABILITIES, EMISSION_PROBABILITIES, observations, PFAM_INDICES)


def forward_backward(transitions: np.array, emissions: np.array,
                     observations: List[str], observation_indices: Dict[str, int]) -> List[float]:
    """ Forward-backward algorithm with precomputed emission and transition probabilities.

        Observation indices serve to find the correct index of an observation into
        the emission probabilities without having to reorder it.

        Arguments:
            transitions: a NxN matrix of floats representing the transition probabilities
            emissions: a NxM matrix of floats representing the emission probabilities
            observations: a list of observations in the order they were observed
            observation_indices: a dictionary mapping observations to the matching
                                 index (0 to M-1) in the emission probabilities

        Returns:
            the final probability for each observation as a list of floats
    """
    # all formulas in comments are in LaTeX format, if not familiar with the
    # syntax, there are online equation editors the formula can be pasted into
    # to render it

    # ensure the matrices are shaped correctly
    assert isinstance(observations, list) and observations
    assert len(transitions.shape) == 2 and transitions.shape[0] == transitions.shape[1]
    num_transitions = transitions.shape[0]
    assert len(emissions.shape) == 2 and emissions.shape[0] == num_transitions
    num_observations = len(observations)

    # extend the emission probabilities by one column to account for observations
    # not accounted for in the index set
    num_emissions = emissions.shape[1]
    emissions = np.hstack((emissions, [[1./num_emissions]] * num_transitions))
    # then point all those missing observations to that final position
    for observation in observations:
        if observation not in observation_indices:
            observation_indices[observation] = num_emissions

    # forward (F), a |O|+1 X |T| matrix, the extra observation adds the ground state
    # formula: F_o_j = \begin{cases}
    #                   E_i_j,&\text{if} i = 0 \\
    #                   F_{i-1}_j * T_i_j * E_i_j,&\text{otherwise}
    #                  \end{cases}
    forward = np.zeros((num_observations + 1, num_transitions))
    # ground state, start probabilities assumed to be 1.0
    emission_index = observation_indices[observations[0]]
    for i in range(num_transitions):
        forward[0][i] = emissions[i, emission_index]
    # following states
    for obs_index in range(1, num_observations):
        emission_index = observation_indices[observations[obs_index]]
        for j in range(num_transitions):
            total = 0.
            for i in range(num_transitions):
                total += forward[obs_index - 1, i] * transitions[i, j] * emissions[j, emission_index]
            forward[obs_index, j] = total
        # the original clusterfinder underflow protection
        # doesn't seem numerically valid, but no clusters are detected without it
        if sum(forward[obs_index]) < 1e20:
            for j in range(num_transitions):
                forward[obs_index, j] *= 1e10

    # backward (B), a |O| X |T| matrix
    # formula: B_i_j = \begin{cases}
    #                   1,& \text{if } i = |T| \\
    #                   B_{i+1}_j \times T_i_j \times E_i_j, &\text{otherwise}
    #                  \end{cases}
    backward = np.zeros((num_observations, num_transitions))
    # base state
    backward[-1] = np.ones((1, num_transitions))
    # previous states, building backwards from the end
    for obs_index in reversed(range(0, num_observations - 1)):
        emission_index = observation_indices[observations[obs_index + 1]]
        for i in range(num_transitions):
            total = 0.
            for j in range(num_transitions):
                total += backward[obs_index + 1, j] * transitions[i, j] * emissions[j, emission_index]
            backward[obs_index][i] = total

        # the original ClusterFinder underflow protection
        # doesn't seem numerically valid, but no clusters are detected without it
        if sum(backward[obs_index]) < 1e20:
            for i in range(num_transitions):
                backward[obs_index, i] *= 1e10

    # merge forward and backward into posterior marginals (P), a |T| X |o| matrix
    # formula: P_o_j = \frac{F_o_j \times B_o_j} {\sum_{i=0}^{|O|}{(F_i_j \times B_i_j)} }
    posterior = np.zeros((num_observations, num_transitions))
    for obs_index in range(num_observations):
        total = sum(forward[obs_index] * backward[obs_index])
        # skip division by zero and let the probablity stay at the 0 it was init as
        if not total:
            continue
        for j in range(num_transitions):
            posterior[obs_index, j] = (forward[obs_index, j] * backward[obs_index, j]) / total

    return list(posterior[:, 0])
