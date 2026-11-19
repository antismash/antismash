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
from antismash.modules.nrps_pks.name_mappings import Category, get_category_by_name, get_substrate_by_name, lookup
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

        hydrophobic_aromatic = get_category_by_name("hydrophobic_aromatic")
        self.three = nrpys.SvmPrediction(hydrophobic_aromatic.name, 1.0, list(hydrophobic_aromatic.substrates))
        self.large = nrpys.SvmPrediction("cys", 1.0, [get_substrate_by_name("Cys")])

        unpolar_aromatic = get_category_by_name("Unpolar aromatic ring")
        self.small = nrpys.SvmPrediction(unpolar_aromatic.name, 1.0, list(unpolar_aromatic.substrates))

        self.single = nrpys.SvmPrediction("ohOrn", 1.0, [get_substrate_by_name("ohOrn")])

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

    def test_bad_signature_lengths(self) -> None:
        with self.assertRaises(ValueError):
            nrpys.PredictorSVMResult(
                "FAKEDATAOK", "FAKEDATAOK", [self.stach], self.three,
                self.large, self.small, self.single)

        with self.assertRaises(ValueError):
            nrpys.PredictorSVMResult(
                "ILIKEDATAEVENFAKEDATADIDISAYILIKED", "ILIKEDATAEVENFAKEDATADIDISAYILIKED",
                [self.stach], self.three, self.large, self.small, self.single)

    def test_classification_certain(self) -> None:
        pred = self.pred
        pred.stachelhaus_quality = 0.6

        assert not pred.uncertain

        assert pred.get_classification() == [self.single.substrates[0].short]

        pred.single_amino = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == [substrate.short for substrate in self.small.substrates]

        pred.small_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == [substrate.short for substrate in self.large.substrates]

        pred.large_cluster = nrpys.SvmPrediction("N/A", 0.0, [])
        assert pred.get_classification() == [substrate.short for substrate in self.three.substrates]

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
                assert pred.get_classification() == [substrate.short for substrate in self.stach.substrates], i

    def test_classification_stach_and_single(self) -> None:
        pred = self.pred
        assert not pred.uncertain

        # stach quality of 1.0 wins
        pred.stachelhaus_quality = 1.0
        assert sorted(pred.get_classification()) == sorted([substrate.short for substrate in self.stach.substrates])

        # stach quality of 0.8 loses
        pred.stachelhaus_quality = 0.8
        assert pred.get_classification() == [self.single.substrates[0].short]

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


@dataclass
class FakePrediction:
    name: str
    score: float = 1.0


class TestSvmPrediction(unittest.TestCase):
    def test_str(self) -> None:
        pred = nrpys.SvmPrediction("hydrophobic-aromatic", 1.0,
                                   list(map(get_substrate_by_name, [
                                    "Phe", "Tyr", "2,3-dohBza", "Pgl", "R-ohTyr",
                                   ])))
        assert str(pred) == "hydrophobic-aromatic: Phe, Tyr, 2,3-dohBza, Pgl, R-ohTyr"


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
        assert ret == [f"Signature file not found '{stach_path}'",
                       f"Model directory not found '{db_path}'"]
