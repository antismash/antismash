# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from dataclasses import dataclass
import unittest
from unittest.mock import patch

from markupsafe import Markup

import antismash
from antismash.config import build_config
from antismash.common.secmet.test.helpers import DummyAntismashDomain
from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name, SubstrateName
from antismash.modules.nrps_pks import nrpys, data_structures


class TestStachelhausMatch(unittest.TestCase):
    def test_str(self) -> None:
        substrates = [get_substrate_by_name(name) for name in ["Ala", "Ibu"]]
        signature = "FAKEDATAOK"
        match = nrpys.StachelhausMatch(substrates, signature, 1.0, 1.0)

        assert str(match) == "Ala or Ibu (FAKEDATAOK)"


class TestPredictorSVMResult(unittest.TestCase):
    def setUp(self) -> None:
        substrates = [get_substrate_by_name(name) for name in ["Ala", "Ibu"]]
        signature = "FAKEDATAOK"
        self.stach = nrpys.StachelhausMatch(substrates, signature, 1.0, 1.0)

        self.three = nrpys.SvmPrediction("hydrophobic-aromatic", 1.0,
                                         list(map(get_substrate_by_name, [
                                          "Phe", "Tyr", "2,3-dohBza", "Pgl", "R-ohTyr",
                                         ])))
        self.large = nrpys.SvmPrediction("cys", 1.0, [get_substrate_by_name("Cys")])

        self.small = nrpys.SvmPrediction("Unpolar aromatic ring", 1.0,
                                         [get_substrate_by_name("Phe"),
                                          get_substrate_by_name("Tyr")])

        self.single = nrpys.SvmPrediction("horn", 1.0, [get_substrate_by_name("ohOrn")])

        self.pred = nrpys.PredictorSVMResult(
            "ILIKEDATAEVENFAKEDATADIDISAYILIKED", "FAKEDATAOK", [self.stach], self.three,
            self.large, self.small, self.single)

    def test_valid(self) -> None:
        pred = self.pred

        assert isinstance(pred, data_structures.Prediction)
        assert pred.aa34 == "ILIKEDATAEVENFAKEDATADIDISAYILIKED"
        assert pred.aa10 == "FAKEDATAOK"
        assert pred.stachelhaus_matches == [self.stach]
        assert pred.stachelhaus_quality == 1.0
        assert not pred.uncertain

    def test_invalid(self) -> None:
        pred = nrpys.PredictorSVMResult(
            "ILIKEDATAEVENFAKE-----------------", "FAKEDATAOK", [self.stach], self.three,
            self.large, self.small, self.single)

        assert pred.uncertain

    def test_classification_certain(self) -> None:
        pred = self.pred
        pred.stachelhaus_quality = 0.6

        assert not pred.uncertain

        assert pred.get_classification() == [self.single.substrates[0].norine]

        pred.single_amino = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == list(
            map(nrpys._get_norine_if_not_x, self.small.substrates))

        pred.small_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == list(
            map(nrpys._get_norine_if_not_x, self.large.substrates))

        pred.large_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == list(
            map(nrpys._get_norine_if_not_x, self.three.substrates))

    def test_classification_uncertain(self) -> None:
        pred = self.pred
        pred.stachelhaus_quality = 0.6

        pred.uncertain = True

        pred.single_amino = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == []

        pred.small_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == []

        pred.large_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == []

        for i in range(10):
            pred.stachelhaus_quality = i / 10
            if i < 8:
                assert pred.get_classification() == []
            else:
                assert pred.get_classification() == list(
                    map(nrpys._get_norine_if_not_x, self.stach.substrates)), i

    def test_classification_stach_and_single(self) -> None:
        pred = self.pred
        assert not pred.uncertain

        # stach quality of 1.0 wins
        pred.stachelhaus_quality = 1.0
        assert sorted(pred.get_classification()) == sorted(list(
                    map(nrpys._get_norine_if_not_x, self.stach.substrates)))

        # stach quality of 0.8 loses
        pred.stachelhaus_quality = 0.8
        assert pred.get_classification() == [self.single.substrates[0].norine]

        pred.stachelhaus_quality = 0.9
        assert not set(pred.single_amino.substrates).intersection(set(self.stach.substrates))
        assert self.pred.get_classification() == []

        pred.single_amino = nrpys.SvmPrediction("ala", 1.0, [get_substrate_by_name("Ala")])
        assert self.pred.get_classification() == ["Ala"]

    def test_classification_stach_and_groups(self) -> None:
        pred = self.pred
        assert not pred.uncertain

        pred.single_amino = nrpys.SvmPrediction("N/A", 0.0, [])
        pred.physicochemical_class = nrpys.SvmPrediction("N/A", 0.0, [])

        pred.stachelhaus_matches = [
            nrpys.StachelhausMatch([get_substrate_by_name("Ala")], pred.aa10, 1.0, 1.0)]

        small = nrpys.SvmPrediction("Tiny, hydrophilic, transition to aliphatic",
                                    0.0, [get_substrate_by_name("Ala"),
                                          get_substrate_by_name("Gly")])
        pred.small_cluster = small

        large = nrpys.SvmPrediction("Hydroxy-benzoic acid derivates",
                                    0.0, [get_substrate_by_name("2,3-dohBza"),
                                          get_substrate_by_name("Sal")])
        pred.large_cluster = large

        for stach_matches in [8, 9]:
            pred.stachelhaus_quality = stach_matches / 10
            pred.small_cluster = small
            assert self.pred.get_classification() == ["Ala"]  # overlap with small

            pred.small_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
            pred.large_cluster = large
            assert self.pred.get_classification() == []  # no overlap with large, but predicted

            pred.large_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
            assert self.pred.get_classification() == ["Ala"]  # no pred from SVM at all

            # and check there's no classification when no SVM pred and stach < 8
            pred.stachelhaus_quality = 0.7
            assert self.pred.get_classification() == []

    def test_classification_stach_34aa_mismatch(self) -> None:
        best = nrpys.StachelhausMatch([get_substrate_by_name("Ser")], "FAKEDATAOK", 1.0, 1.0)
        second_best = nrpys.StachelhausMatch([get_substrate_by_name("Ala")], "FAKEDATAOK", 1.0, 0.9)
        pred = nrpys.PredictorSVMResult(
            "ILIKEDATAEVENFAKEDATADIDISAYILIKED", "FAKEDATAOK", [best, second_best], self.three,
            self.large, self.small, self.single)
        assert pred.get_classification() == ["Ser"]

    def test_classification_stach_34aa_partial(self) -> None:
        best = nrpys.StachelhausMatch([get_substrate_by_name("Ser")], "FAKEDATAOK", 1.0, 1.0)
        second_best = nrpys.StachelhausMatch([get_substrate_by_name("Asp")], "FAKEDATAOK", 1.0, 0.9)
        also_best = nrpys.StachelhausMatch([get_substrate_by_name("Ala")], "FAKEDATAOK", 1.0, 1.0)
        pred = nrpys.PredictorSVMResult(
            "ILIKEDATAEVENFAKEDATADIDISAYILIKED", "FAKEDATAOK", [second_best, best, also_best], self.three,
            self.large, self.small, self.single)
        assert pred.get_classification() == ["Ser", "Ala"]

    def test_html(self) -> None:
        self.pred.uncertain = False
        html = self.pred.as_html()
        assert isinstance(html, Markup)
        assert "uncertain match" not in str(html)

        self.pred.uncertain = True
        html = self.pred.as_html()
        assert "uncertain match" in str(html)

        self.pred.stachelhaus_quality = 0.8
        html = self.pred.as_html()
        assert "(moderate)" in str(html)

    def test_str(self) -> None:
        assert "PredictorSVMResult" in str(self.pred)


NAME_TO_SUBSTRATES: dict[str, tuple[str, list[SubstrateName]]] = {
    "hydrophilic": ("hydrophilic", list(map(get_substrate_by_name, [
        "Arg", "Asp", "Glu", "Asn", "Lys", "Gln",
        "Orn", "Aad",
    ]))),
    "hydrophobic-aliphatic": ("hydrophobic-aliphatic", list(map(get_substrate_by_name, [
        "Ala", "Gly", "Val", "Leu", "Ile", "Abu",
        "Iva", "Ser", "Thr", "Hpg", "dHpg", "Cys",
        "Pro", "Pip",
    ]))),
    "hydrophobic-aromatic": ("hydrophobic-aromatic", list(map(get_substrate_by_name, [
        "Phe", "Tyr", "2,3-dohBza", "Pgl", "R-ohTyr",
    ]))),
    "asp,asn,glu,gln,aad": ("Aliphatic chain with H-bond donor", list(map(get_substrate_by_name, [
        "Asp", "Asn", "Glu", "Gln", "Aad",
    ]))),
    "cys": ("Polar, uncharged (aliphatic with -SH)", [get_substrate_by_name("Cys")]),
    "dhb,sal": ("Hydroxy-benzoic acid derivates", [get_substrate_by_name("2,3-dohBza"),
                                                   get_substrate_by_name("Sal")]),
    "gly,ala,val,leu,ile,abu,iva": ("Apolar, aliphatic", list(map(get_substrate_by_name,
                                    "gly,ala,val,leu,ile,abu,iva".split(",")))),
    "orn,lys,arg": ("Long positively charged side chain", list(map(get_substrate_by_name,
                    "orn,lys,arg".split(",")))),
    "phe,trp,phg,tyr,bht": ("Aromatic side chain", list(map(get_substrate_by_name,
                            "phe,trp,Pgl,tyr,R-ohTyr".split(",")))),
    "pro,pip": ("Cyclic aliphatic chain (polar NH2 group)", list(map(get_substrate_by_name,
                "pro,pip".split(",")))),
    "ser,thr,dhpg,hpg": ("Aliphatic chain or phenyl group with -OH", list(map(get_substrate_by_name,
                         "ser,thr,dhpg,hpg".split(",")))),
    "dhpg,hpg": ("Polar, uncharged (hydroxy-phenyl)", list(map(get_substrate_by_name,
                 "dhpg,hpg".split(",")))),
    "gly,ala": ("Tiny, hydrophilic, transition to aliphatic", list(map(get_substrate_by_name,
                "gly,ala".split(",")))),
    "orn,horn": ("Orn and hydroxy-Orn specific", list(map(get_substrate_by_name,
                 "orn,ohOrn".split(",")))),
    "phe,trp": ("Unpolar aromatic ring", list(map(get_substrate_by_name,
                "phe,trp".split(",")))),
    "tyr,bht": ("Polar aromatic ring", list(map(get_substrate_by_name,
                "tyr,R-ohTyr".split(",")))),
    "val,leu,ile,abu,iva": ("Aliphatic, branched hydrophobic", list(map(get_substrate_by_name,
                            "val,leu,ile,abu,iva".split(",")))),

}


@dataclass
class FakePrediction:
    name: str
    score: float = 1.0


class TestSvmPrediction(unittest.TestCase):
    def test_resolve(self) -> None:
        for init, (name, substrates) in NAME_TO_SUBSTRATES.items():
            assert nrpys._resolve_name_substrates(FakePrediction(init)) == (name, substrates)


class TestMisc(unittest.TestCase):
    def setUp(self) -> None:
        self.config = build_config([], isolated=True, modules=antismash.get_all_modules())

    @patch("nrpys.run", return_value=[])
    @patch.object(nrpys, "get_a_dom_signatures", side_effect=[(None, None)])
    def test_run_nrpys(self, _mock_signature, _mock_run) -> None:
        ret = nrpys.run_nrpys([
            DummyAntismashDomain(locus_tag="A", domain_id="fakeA1")], self.config)
        assert ret == {}

    @patch.object(nrpys, "find_latest_database_version", side_effect=["1.0", "1.0"])
    def test_get_signature_path_fallback(self, _mock_version_lookup) -> None:
        config = build_config([
            "--databases", "/fake/path"], isolated=True, modules=antismash.get_all_modules())
        assert "1.0" in nrpys._get_signature_path(config)

        assert "1.1" in nrpys._get_signature_path(config, "1.1")


    @patch("os.path.exists", side_effect=[False, False])
    def test_check_prereqs_errors(self, _exists_mock) -> None:
        stach_path = nrpys._get_signature_path(self.config)
        db_path = nrpys._get_model_dir(self.config)
        ret = nrpys.check_prereqs(self.config)
        assert ret == [f"Failed to locate {stach_path}",
                       f"Failed to locate {db_path}"]
