# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" transPACT analysis of KS domains """

from typing import Any, Dict, List, Tuple
from io import StringIO
import re
import copy
from jinja2 import Markup

from antismash.common import path, subprocessing
from antismash.modules.nrps_pks.data_structures import Prediction

from Bio import Phylo




_LEAF2CLADE_TBL = path.get_full_path(__file__, "data", "transPACT_leaf2clade.tsv")
_PPLACER_MASS_CUTOFF = 0.6  # transPACT default: 0.6; higher = more stringent

# Note: Reference package creation with the taxtastic package (install with pip):  taxit create --aln-fasta
# transAT_KS_ref.afa --tree-stats transAT_KS_ref.raxml.info --tree-file transAT_KS_ref.raxml.tre -P transAT_KS_refpkg
# -l transAT_KS
_PPLACER_REFERENCE_PKG = path.get_full_path(__file__, "data", "transAT_KS_refpkg")
_KS_REFERENCE_ALIGNMENT = path.get_full_path(__file__, "data", "transAT_KS_refpkg", 'transAT_KS_ref.afa')
_KS_REFERENCE_TREE = path.get_full_path(__file__, "data", "transAT_KS_refpkg", 'transAT_KS_ref.raxml.tre')


class KSResult:
    """ A result for a specific KS domain """
    __slots__ = ["clade", "specificity", "mass_score"]

    def __init__(self, clade: str, specificity: str, mass_score: float) -> None:
        assert isinstance(clade, str)
        assert isinstance(specificity, str)
        assert isinstance(mass_score, float)
        self.clade = clade
        self.specificity = specificity
        self.mass_score = mass_score

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return "KSResult(clade=%s, specificity=%s, mass_score=%s)" % (self.clade, self.specificity, self.mass_score)

    def to_json(self) -> Tuple[str, str, float]:
        """ Serialises the instance """
        return (self.clade, self.specificity, self.mass_score)

    @classmethod
    def from_json(cls, json: Tuple[str, str, float]) -> "KSResult":
        """ Deserialise a KSResult instance """
        assert len(json) == 3
        return cls(*json)


class KSPrediction(Prediction):
    """ Holds the transPACT predictions for a domain"""
    def __init__(self, predictions: Dict[str, KSResult]) -> None:
        super().__init__("transPACT")
        self.predictions = sorted(predictions.items(), key=lambda x: (-x[1].mass_score, x[1].clade))

    def get_classification(self) -> List[str]:
        results = []  # type: List[str]
        preds = []
        if not self.predictions:
            return results
        for spec, pred in self.predictions:
            results.append(spec)
            preds.append(pred)
        return results

    def as_html(self) -> Markup:
        if not self.predictions:
            return Markup("No matches")
        lines = []
        for spec, pred in self.predictions:
            desc = spec+' ('+pred.clade+')'
            lines.append("<dd>%s: %.1f%%</dd>\n" % (desc, pred.mass_score))
        html = ((
            "<dl>\n"
            " <dt>transPACT assigned specificiy:</dt>\n"
            "%s"
            "</dl>\n"
        ) % "".join(lines))
        return Markup(html)

    def to_json(self) -> Dict[str, Any]:
        return {
            "method": "transPACT",
            "predictions": {clade: pred.to_json() for clade, pred in self.predictions},
        }

    @classmethod
    def from_json(cls, json: Dict[str, Any]) -> "KSPrediction":
        assert json["method"] == "transPACT"
        return cls({specificity: KSResult.from_json(pred) for specificity, pred in json["predictions"].items()})


def get_leaf2clade(leaf2cladetbl: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """ Generates the dicts to map between reference
    tree leaf names, clade, and verbose annotation """
    leaf2clade = {}
    clade2ann = {}
    with open(leaf2cladetbl) as cfh:
        for ln in cfh.read().splitlines():
            ksname, clade, ann = ln.split("\t")
            leaf2clade[ksname] = clade
            clade2ann[clade] = ann.replace(' ', '_')
    return(leaf2clade, clade2ann)


def choose_elder(parent: Any, grandparent: Any) -> Any:
    """ Selects either parent or grandparent for monophyly test """
    ## Find the proper clade elder, if more than two siblings use parent, otherwise use grandparent
    if parent is None:
        return None
    else:
        if len(parent.get_terminals()) > 2:
            return parent
        else:
            if grandparent is None:
                return None
            else:
                return grandparent


def get_transpact_clade(query_name: str, tree: Any, funclades: Dict[str, str]) -> str:
    """
    tree: Bio.Phylo.Newick.Tree
    """
    ## Remove placements that aren't the query in question and note query's elders
    newtree = copy.deepcopy(tree) ## necessary to not screw up original tree for subsequent queries
    parent, grandparent = None, None
    for leaf in newtree.get_terminals():
        if leaf.name == query_name:
            node_path = newtree.get_path(leaf)
            parent, grandparent = node_path[-2], node_path[-3]
        else:
            ln = leaf.name.split("_")
            if re.match("^#\d+$", ln[-2]) is not None:
                newtree.prune(leaf)
    elder = choose_elder(parent, grandparent)
    if elder is None:
        raise ValueError("No leaf named %s in tree terminals." % query_name)
    ## Count number of occurances of each clade in elder descendants
    clade_count = {}  # type: Dict[str, int]
    for leaf in elder.get_terminals():
        if leaf.name != query_name:
            if funclades[leaf.name] in clade_count:
                clade_count[funclades[leaf.name]] += 1
            else:
                clade_count[funclades[leaf.name]] = 1
    clade_assignment = 'clade_not_conserved'
    if len(clade_count) == 1: ## clade consensus, monophyly
        clade_assignment = list(clade_count)[0]
    return clade_assignment


def get_best_transpact(totalmass: Dict[str, float],
                       masscutoff: float) -> Tuple[str, float]:
    """ Returns clade, spec and score for a transPACT query """
    best_clade: str = 'None'
    best_mass: float = float(0)
    for cld in totalmass:
        if totalmass[cld] > best_mass:
            best_clade, best_mass = cld, totalmass[cld]
    clade: str = 'clade_not_conserved'
    score: float = 0.0
    if best_clade != 'clade_not_conserved' and best_mass >= masscutoff:
        clade = best_clade
        score = round(best_mass, 2)
    return clade, score


def get_totalmass(tree: Any,
                  tree_hits: Dict[str, Any],
                  funclades: Dict[str, str]) -> Dict[str, float]:
    """ Returns the totalmass dict for each clade """
    totalmass: Dict[str, float] = {}
    query_prefix = None
    for placement_num in tree_hits:
        if query_prefix is None:
            query_prefix = re.sub(r"^(.+)_#\d+_M=\d+?\.?\d*$", "\g<1>", tree_hits[placement_num].name)
        mass = float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement_num].name))
        clade_assignment = get_transpact_clade(tree_hits[placement_num].name, tree, funclades)
        if clade_assignment in totalmass:
            totalmass[clade_assignment] += mass
        else:
            totalmass[clade_assignment] = mass
    return totalmass


def transpact_tree_prediction(pplacer_tree: str,
                              masscutoff: float,
                              funclades: Dict[str, str],
                              clade2ann: Dict[str, str]) -> KSPrediction:
    """ From a tree, calculate the max pplacer mass
    for each monophyletic specificity """
    tree = Phylo.read(StringIO(pplacer_tree), 'newick')
    tree_hits = {}
    for leaf in tree.get_terminals():
        if not leaf.name:
            continue
        ln = leaf.name.split("_") ## name ensures not truncated
        if re.match("^#\d+$", ln[-2]) is not None: ## Fits pplacer format
            placementn = re.sub(r"^#(\d+)$", "\g<1>", ln[-2]) ## placement number, zero indexed
            tree_hits[placementn] = leaf
            funclades[leaf.name] = 'query_seq'
    if not tree_hits: ## len = 0
        raise ValueError("There should be a leaf with name of minimal form #'int'_M='float' in the provided tree.")
    ## Look to see when threshold is met
    totalmass = get_totalmass(tree, tree_hits, funclades)
    clade, score = get_best_transpact(totalmass, masscutoff)
    spec = 'NA'
    if clade in clade2ann:
        spec = clade2ann[clade]
    return KSPrediction({spec: KSResult(clade, spec, score)})


def run_transpact_pplacer(ks_name: str,
                          alignment: Dict[str, str],
                          reference_tuple: Tuple[str, str, str],
                          masscutoff: float,
                          funclades: Dict[str, str],
                          clade2ann: Dict[str, str]) -> Prediction:
    """ Calls pplacer and runs prediction pipeline """
    pplacer_tree = subprocessing.run_pplacer(ks_name, alignment, reference_tuple[0], reference_tuple[1], reference_tuple[2])
    prediction = transpact_tree_prediction(pplacer_tree, masscutoff, funclades, clade2ann)
    return prediction


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
    funclades, clade2ann = get_leaf2clade(_LEAF2CLADE_TBL)
    reference_tuple = (_PPLACER_REFERENCE_PKG, _KS_REFERENCE_ALIGNMENT, _KS_REFERENCE_TREE)
    results = {}
    for ks_name, ks_seq in domains.items():
        ## Align to reference
        alignment = subprocessing.run_muscle_single(ks_name, ks_seq, _KS_REFERENCE_ALIGNMENT)
        results[ks_name] = run_transpact_pplacer(ks_name,
                                                 alignment,
                                                 reference_tuple,
                                                 _PPLACER_MASS_CUTOFF,
                                                 funclades,
                                                 clade2ann)
    return results
