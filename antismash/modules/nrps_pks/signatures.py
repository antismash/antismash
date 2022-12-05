# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Handle extracting the 10 AA and 34 AA signatures of NRPS adenylation domains.
"""

from antismash.common import path, subprocessing, utils
from antismash.detection.nrps_pks_domains import ModularDomain

ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
# the NRPSPredictor2 HMM for aligning A domains
ADOMAINS_HMM = path.get_full_path(__file__, "data", "aa-activating.aroundLys.hmm")
# stachelhaus/10AA positions as per the NRPSPredictor2 HMM
HMM_SITE_POSITIONS_10 = [46, 47, 50, 92, 124, 126, 153, 161, 162]
# 34AA/8A positions as per the NRPSPredictor2 HMM
HMM_SITE_POSITIONS_34 = [
    12, 15, 16, 40, 45, 46, 47, 48, 49, 50, 51, 54, 92, 93, 124, 125, 126, 127,
    128, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165,
]


def verify_good_sequence(sequence: str) -> bool:
    """ Ensures a sequence is valid """
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True


def get_a_dom_signatures(domain: ModularDomain) -> tuple[str, str]:
    """ Extract 10 / 34 AA NRPS signatures from A domains """
    assert " " not in domain.get_name()
    assert verify_good_sequence(domain.translation)

    args = ["-T", "0", "-E", "0.1"]  # min score and max evalue
    results = subprocessing.hmmpfam.run_hmmpfam2(ADOMAINS_HMM, f">query\n{domain.translation}",
                                                 extra_args=args)
    if not (results and results[0].hsps):
        raise RuntimeError("no hits for query, is this a legitimate A-domain?")

    hit = results[0].hsps[0]
    profile = hit.aln[1].seq
    query = hit.aln[0].seq
    offset = hit.hit_start

    aa10 = utils.extract_by_reference_positions(query, profile, [p - offset for p in HMM_SITE_POSITIONS_10])
    aa10 += "K"  # since only 9 positions are included and the 10th is generally K
    aa34 = utils.extract_by_reference_positions(query, profile, [p - offset for p in HMM_SITE_POSITIONS_34])

    return aa10, aa34
