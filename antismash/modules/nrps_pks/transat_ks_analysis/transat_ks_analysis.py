# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" transPACT analysis of KS domains """

from typing import Any, Dict, List, Tuple

from jinja2 import Markup

from antismash.common import path, subprocessing, utils, fasta
from antismash.modules.nrps_pks.data_structures import Prediction
from antismash.modules.nrps_pks.pks_names import get_long_form

_LEAF2CLADE_TBL = path.get_full_path(__file__, "data", "transPACT_leaf2clade.tsv")
_PPLACER_MASS_CUTOFF = 0.6
_PPLACER_REFERENCE_PKG = path.get_full_path(__file__, "data", "RAxML_bestTree.649KS_sequences_hmmalign_raxml_renamed.refpkg") ## Note: Reference package creation: taxit create --aln-fasta reference.fasta --tree-stats reference.log --tree-file reference.nwk -P reference.refpkg
_KS_REFERENCE_ALIGNMENT = '/'.join([_PPLACER_REFERENCE_PKG, '649KS_sequences_031218.fasta'])
_KS_REFERENCE_TREE = '/'.join([_PPLACER_REFERENCE_PKG, 'RAxML_bestTree.649KS_sequences_hmmalign_raxml_renamed.tre'])

class KSResult:
    """ A result for a specific KS domain """
    __slots__ = ["name", "clade", "specificity"]

    def __init__(self, name: str, clade: str, specificity: str) -> None:
        assert isinstance(name, str)
        assert isinstance(clade, str)
        assert isinstance(specificity, str)
        self.name = name
        self.clade = clade
        self.specificity = specificity

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "KSResult(name=%s, clade=%s, specificity=%s)" % (self.name, self.clade, self.specificity)

    def to_json(self) -> Tuple[str, str, str]:
        """ Serialises the instance """
        return (self.name, self.clade, self.specificity)

    @staticmethod
    def from_json(json: Tuple[str, str, str]) -> "KSResult":
        """ Deserialise an KSResult instance """
        assert len(json) == 3
        return KSResult(*json)


class KSPrediction(Prediction):
    """ Holds the transPACT predictions for a domain"""
    def __init__(self, prediction: str) -> None:
        super().__init__("transPACT_KS")
        self.prediction = prediction
        
    def as_html(self) -> Markup:
        if not self.prediction:
            return Markup("No matches")
        line = "<dd>%s</dd>\n" % (self.prediction)
        html = ((
            "<dl>\n"
            " <dt>transPACT assigned specificiy:</dt>\n"
            "%s"
            "</dl>\n"
        ) % line)
        return Markup(html)

    def to_json(self) -> Dict[str, Any]:
        return {
            "method": "transPACT_KS",
            "prediction": self.prediction,
        }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "KSPrediction":
        assert json["method"] == "transPACT_KS"
        return KSPrediction(json["prediction"])
    
def get_leaf2clade(leaf2cladetbl: str) -> [Dict[str, str], Dict[str, str]]:
    leaf2clade = {}
    clade2ann = {}
    with open(leaf2cladetbl) as c:
        for ln in c.read().splitlines():
            ksname, clade, ann = ln.split("\t")
            leaf2clade[ksname] = clade
            clade2ann[clade] = ann
    return(leaf2clade, clade2ann)


def run_transpact_pplacer(ks_name: str, alignment: Dict[str, str], reference_pkg: str, reference_aln: str, reference_tree: str, masscutoff: float, funClades: Dict[str, str]) -> Dict[str, str]:
    
    pplacer_tree = subprocessing.run_pplacer(ks_name, alignment, reference_pkg, reference_aln, reference_tree)
    print(pplacer_tree)
    import sys
    sys.exit("post pplacer exit")
    clade_assignment = parse_pplacer(ks_name, data_dir, masscutoff, funClades)
    for q in clade_assignment:
        out.write('%s: %s\n' % (q, clade_assignment[q]))
    out.close()
    return clade_assignment
    
    
def run_transpact_ks_analysis(domains: Dict[str, str]) -> Dict[str, Prediction]:
    """ Analyses PKS signature of KS domains

        Arguments:
            domains: a dictionary mapping domain identifier (e.g. 'locus_KS2')
                     to domain sequence

        Returns:
            a dictionary mapping domain identifier to
                a list of KSResults
    """
    ## Read clade to annotation maps from flat files
    funClades, clade2ann = get_leaf2clade(_LEAF2CLADE_TBL)
    
    for ks_name, ks_seq in domains.items():
        ## Align to reference
        alignment = subprocessing.run_muscle_single(ks_name, ks_seq, _KS_REFERENCE_ALIGNMENT)
        clade = run_transpact_pplacer(ks_name, alignment, _PPLACER_REFERENCE_PKG, _KS_REFERENCE_ALIGNMENT, _KS_REFERENCE_TREE, _PPLACER_MASS_CUTOFF, funClades)
    pass
