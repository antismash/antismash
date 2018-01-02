# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import re

from antismash.common import secmet, utils


def get_signature(query, hmm, positions, expected=None):
    ungapped = str(hmm).replace('.', '')
    if positions[-1] > len(ungapped):
        logging.warning("scaffold coordinate %s outside hsp!", positions[-1] - 1)
        return None, None
    ref_signature = "".join(ungapped[pos - 1] for pos in positions)
    if expected:
        assert ref_signature == "".join(expected)
    signature = utils.extract_by_reference_positions(query, hmm, [pos - 1 for pos in positions])
    return signature


def get_scaffold_annotation(result, positions, values, expected_emissions=None):
    "generate annotation from scaffold information"

    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq
    assert len(values) == len(positions)
    if not expected_emissions:
        emissions = ["n.d."] * len(positions)
    else:
        emissions = list(map(int, expected_emissions))  # TODO: should floats (0., 1.) really be converted to ints?
    assert len(emissions) == len(positions)

    # Check scaffold matches
    matches = []

    signature = get_signature(query_seq, hmm_seq, positions, values)
    assert len(signature) == len(positions), "%d != %d" % (len(signature), len(positions))
    for i, char in enumerate(signature):
        position = positions[i] - 1
        expected = values[i]

        # We have to use a RegEx here to allow negations and more complex queries; ignore case (?i)
        match = bool(re.match("(?i)" + expected, char))
        matches.append(str(match))
        logging.debug("Scaffold coordinate %s; query aa %s; "
                      "scaffold value %s; emission probability %s; match %s",
                      position, char, expected, emissions[i], match)

    overall_match = all(matches)

    logging.debug("Overall Scaffold Match: %s", str(overall_match).upper())

    # Generate Feature qualifiers
    return ("Scaffold coordinates: (%s); scaffold residues: (%s); expected: (%s); matchArray: (%s); "
            "emission probability array (%s); overall match: %s") % (
            ",".join(map(str, positions)), signature, ",".join(values),
            ",".join(map(str, matches)), ",".join(map(str, emissions)),
            str(overall_match).upper())


def get_prediction_annotation(result, positions, values, result_label, comment, expected_emissions=None):
    "gererate annotation from choices/prediction information, a single 'choice' at a time"

    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq
    if not expected_emissions:
        emissions = ["n.d."] * len(positions)
    else:
        emissions = list(map(int, expected_emissions))  # TODO: again, should floats (0., 1.) really be converted to ints?

    signature = get_signature(query_seq, hmm_seq, positions, values)

    logging.debug("testing %s (%s):", result_label, comment)

    matches = []
    for i, char in enumerate(signature):
        value = values[i]

        position = int(positions[i]) - 1

        matches.append(bool(re.match("(?i)" + value, char)))

        logging.debug("Offset %s; Expected: %s; observed in query %s; Emission %s; match %s",
                      position,
                      value,
                      char,
                      emissions[i],
                      matches[-1])

    overall_match = all(matches)

    logging.debug("Overall Match for prediction %s: %s", result_label, str(overall_match).upper())
    logging.debug("================================")

    description = (("Description: %s, choice result: %s, choice coordinates: (%s); residues: (%s); "
                   "expected for choice: (%s); matchArray: (%s); emission probability array (%s); overall match: %s") %
                   (comment,
                    result_label,
                    ",".join(map(str, positions)),
                    signature,
                    ",".join(values),
                    ",".join(map(str, matches)),
                    ",".join(map(str, emissions)),
                    str(overall_match).upper()))

    choice_string = ""
    if overall_match:
        choice_string = "Full match for prediction: %s" % result_label
    return description, choice_string
