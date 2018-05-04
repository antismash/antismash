# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
import sys
import tempfile
from typing import Any, Dict, List

from helperlibs.wrappers.io import TemporaryDirectory
from jinja2 import Markup

from antismash.common import fasta, path, subprocessing
from antismash.common.secmet import AntismashDomain
from antismash.common.secmet.qualifiers import DomainPrediction
from antismash.config import ConfigType

from .data_structures import Prediction

REF_SEQUENCE = "P0C062_A1"
A34_POSITIONS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A34positions.txt")
APOSITION_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "Apositions.txt")
KNOWN_CODES = path.get_full_path(__file__, "knowncodes.fasta")
ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
ADOMAINS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A_domains_muscle.fasta")
START_POSITION = 66


class PredictorSVMResult(Prediction):
    def __init__(self, angstrom_code: str, physicochemical_class: str, large_cluster_pred: List[str],
                 small_cluster_pred: List[str], single_amino_pred: str, stachelhaus_code: str, uncertain: bool) -> None:
        self.angstrom_code = str(angstrom_code)
        self.physicochemical_class = str(physicochemical_class)
        self.large_cluster_pred = list(large_cluster_pred)
        self.small_cluster_pred = list(small_cluster_pred)
        self.single_amino_pred = str(single_amino_pred)
        assert ',' not in self.single_amino_pred
        self.stachelhaus_code = str(stachelhaus_code)
        assert ',' not in self.stachelhaus_code
        self.uncertain = bool(uncertain)

    def get_classification(self) -> List[str]:
        classification = ["N/A"]
        if self.uncertain:
            return classification
        for prediction in [self.single_amino_pred, self.small_cluster_pred,  # TODO: include stach code?
                           self.large_cluster_pred, self.physicochemical_class]:
            if prediction not in ["N/A", ["N/A"]]:
                if isinstance(prediction, str):
                    classification = [prediction]
                else:
                    classification = list(prediction)
                break
        return classification

    def as_html(self) -> Markup:
        raw = ("\n"
               "<dl><dt>NRPSPredictor2 SVM prediction details:</dt>\n"
               " <dd>"
               "  %s"
               "  <dl>"
               "   <dt>Predicted physicochemical class:</dt>\n"
               "   <dd>%s</dd>\n"
               "   <dt>Large clusters prediction:</dt>\n"
               "   <dd>%s</dd>\n"
               "   <dt>Small clusters prediction:</dt>\n"
               "   <dd>%s</dd>\n"
               "   <dt>Single AA prediction:</dt>\n"
               "   <dd>%s</dd>\n"
               "   <dt>Nearest stachelhaus code:</dt>\n"
               "   <dd>%s</dd>\n"
               "  </dl>\n"
               "  </dd>\n"
               "</dl>\n" % ("NOTE: outside applicability domain<br>\n" if self.uncertain else "",
                            self.physicochemical_class, ", ".join(self.large_cluster_pred),
                            ", ".join(self.small_cluster_pred), self.single_amino_pred,
                            self.stachelhaus_code))
        return Markup(raw)

    def to_domain_prediction_qualifier(self) -> DomainPrediction:
        return DomainPrediction(",".join(self.get_classification()), "NRPSPredictor2")

    @classmethod
    def from_line(cls, line: str) -> "PredictorSVMResult":
        """ Generates a PredictorSVMResult from a line of NRPSPredictor2 output """
        parts = line.split("\t")
        # 0: sequence-id
        # 1: 8A-signature
        # 2: stachelhaus-code:
        # 3: 3class-pred
        # 4: large-class-pred
        # 5: small-class-pred
        # 6: single-class-pred
        # 7: nearest stachelhaus code
        # 8: NRPS1pred-large-class-pred
        # 9: NRPS2pred-large-class-pred
        # 10: outside applicability domain (1 or 0)
        # 11: coords
        # 12: pfam-score
        if not len(parts) == 13:
            raise ValueError("Invalid SVM result line: %s" % line)
        return cls(parts[1], parts[3], parts[4].split(","), parts[5].split(","),
                   parts[6], parts[7], parts[10] == "1")

    def __str__(self) -> str:
        return "PredictorSVMResult: " + str(vars(self))

    def to_json(self) -> Dict[str, Any]:
        return vars(self)

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "PredictorSVMResult":
        return PredictorSVMResult(json["angstrom_code"], json["physicochemical_class"],
                                  json["large_cluster_pred"], json["small_cluster_pred"],
                                  json["single_amino_pred"], json["stachelhaus_code"],
                                  json["uncertain"])


def read_positions(filename: str, start_position: int) -> List[int]:
    data = open(filename, "r")
    text = data.read().strip()
    results = []
    for i in text.split("\t"):
        results.append(int(i) - start_position)
    data.close()
    return results


def build_position_list(positions: List[int], reference_seq: str) -> List[int]:
    poslist = []
    position = 0
    for i, ref in enumerate(reference_seq):
        if ref != "-":
            if position in positions:
                poslist.append(i)
            position += 1
    return poslist


def verify_good_sequence(sequence: str) -> bool:
    for char in ILLEGAL_CHARS:
        if char in sequence:
            return False
    return True


def extract(query_seq: str, poslist: List[int]) -> str:
    seq = []
    for position in poslist:
        aa = query_seq[position]
        if aa == "-":
            if position - 1 not in poslist:
                aa = query_seq[position - 1]
            elif position + 1 not in poslist:
                aa = query_seq[position + 1]
        seq.append(aa)
    return "".join(seq)


def run_muscle(domain: AntismashDomain) -> Dict[str, str]:
    """ run muscle over the query and known A domains """
    # using mkstemp() here because NamedTemporaryFiles cause muscle
    # to think the file is infinite and use > 20Gb of RAM per invocation
    in_fd, in_filename = tempfile.mkstemp(prefix="antismash_musc")
    in_handle = os.fdopen(in_fd, "w")
    in_handle.write(">{}\n{}\n".format(domain.get_name(), domain.translation))
    in_handle.close()
    out_fd, out_filename = tempfile.mkstemp(prefix="antismash_musc")
    os.close(out_fd)
    try:
        res = subprocessing.execute(["muscle", "-profile", "-quiet",
                                     "-in1", ADOMAINS_FILENAME,
                                     "-in2", in_filename,
                                     "-out", out_filename])
        if not res.successful():
            logging.error(res.stderr)
            raise RuntimeError("Error while running muscle: %s" % res.stderr)
    except:
        raise
    else:
        results = fasta.read_fasta(out_filename)
    finally:
        os.unlink(in_filename)
        os.unlink(out_filename)
    return results


def get_34_aa_signature(domain: AntismashDomain) -> str:
    """ Extract 10 / 34 AA NRPS signatures from A domains """
    assert " " not in domain.get_name()
    assert verify_good_sequence(domain.translation)
    # Run muscle and collect sequence positions from file
    alignments = run_muscle(domain)

    domain_alignment = alignments[domain.get_name()]
    reference_alignment = alignments[REF_SEQUENCE]

    positions = read_positions(APOSITION_FILENAME, START_POSITION)
    # Count residues in ref sequence and put positions in list
    poslist = build_position_list(positions, reference_alignment)

    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    # repeat with 34 AA codes
    angpositions = read_positions(A34_POSITIONS_FILENAME, START_POSITION)
    poslist = build_position_list(angpositions, reference_alignment)

    return extract(domain_alignment, poslist)


def run_nrpspredictor(a_domains: List[AntismashDomain], options: ConfigType) -> List[PredictorSVMResult]:
    # NRPSPredictor: extract AMP-binding + 120 residues N-terminal of this domain,
    # extract 8 Angstrom residues and insert this into NRPSPredictor
    nrps_predictor_dir = path.get_full_path(__file__, "external", "NRPSPredictor2")
    data_dir = os.path.join(nrps_predictor_dir, 'data')
    lib_dir = os.path.join(nrps_predictor_dir, 'lib')
    jar_file = os.path.join(nrps_predictor_dir, 'build', 'NRPSpredictor2.jar')
    java_separator = ":"
    if sys.platform == "win32":
        java_separator = ";"
    classpath = java_separator.join([jar_file,
                                     os.path.join(lib_dir, 'java-getopt-1.0.13.jar'),
                                     os.path.join(lib_dir, 'Utilities.jar'),
                                     os.path.join(lib_dir, 'libsvm.jar')])
    input_filename = "signatures.fa"
    output_filename = "svm_output.txt"
    bacterial = "1" if options.taxon == "bacteria" else '0'

    signatures = [get_34_aa_signature(a_domain) for a_domain in a_domains]

    with TemporaryDirectory(change=True):
        # Get NRPSPredictor2 code predictions, output sig file for input for NRPSPredictor2 SVMs
        with open(input_filename, "w") as handle:
            for sig, domain in zip(signatures, a_domains):
                handle.write("%s\t%s\n" % (sig, domain.get_name()))
        # Run NRPSPredictor2 SVM
        commands = ['java',
                    '-Ddatadir=%s' % data_dir,
                    '-cp', classpath,
                    'org.roettig.NRPSpredictor2.NRPSpredictor2',
                    '-i', input_filename,
                    '-r', output_filename,
                    '-s', '1',
                    '-b', bacterial]
        result = subprocessing.execute(commands)
        if not result.successful():
            raise RuntimeError("NRPSPredictor2 failed: %s" % result.stderr)

        with open(output_filename) as handle:
            lines = handle.read().splitlines()[1:]

    return [PredictorSVMResult.from_line(line) for line in lines]
