# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from jinja2 import Markup
from minimock import mock, restore

from antismash.common import path, fasta, subprocessing  # mocked # pylint: disable=unused-import
from antismash.common.secmet import AntismashDomain, FeatureLocation
from antismash.modules.nrps_pks import nrps_predictor, data_structures


def read_file():
    output = path.get_full_path(__file__, 'data', 'nrps_pred2_Y16952.txt')
    with open(output) as handle:
        content = handle.read()
    assert content.startswith('#')
    assert content.count('\n') == 9
    return content.splitlines()[1:]


class TestResultsReading(unittest.TestCase):
    def setUp(self):
        self.lines = read_file()

    def test_single(self):
        results = nrps_predictor.read_output(self.lines[0:1])
        assert list(results) == ["nrpspksdomains_bpsA_A1"]
        pred = results["nrpspksdomains_bpsA_A1"]
        assert isinstance(pred, data_structures.Prediction)
        assert pred.angstrom_code == "L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW"
        assert pred.physicochemical_class == "hydrophobic-aliphatic"
        assert pred.large_cluster_pred == ["N/A"]
        assert pred.small_cluster_pred == ['val', 'leu', 'ile', 'abu', 'iva']
        assert pred.single_amino_pred == "leu"
        assert pred.stachelhaus_predictions == ["leu"]
        assert pred.uncertain is False

        assert pred.get_classification() == [pred.single_amino_pred]

    def test_multiple(self):
        results = nrps_predictor.read_output(self.lines)
        assert len(results) == 8
        for name, single in [("nrpspksdomains_bpsA_A1", "leu"),
                             ("nrpspksdomains_bpsA_A2", "bht"),
                             ("nrpspksdomains_bpsA_A3", "asn"),
                             ("nrpspksdomains_bpsB_A1", "hpg"),
                             ("nrpspksdomains_bpsB_A2", "hpg"),
                             ("nrpspksdomains_bpsB_A3", "bht"),
                             ("nrpspksdomains_bpsC_A1", "dhpg"),
                             ("nrpspksdomains_bpsD_A1", "tyr")]:
            assert results[name].single_amino_pred == single, name


class TestPrediction(unittest.TestCase):
    def setUp(self):
        self.pred = nrps_predictor.PredictorSVMResult("L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW",
                                                      "hydrophobic-aliphatic",
                                                      ['asp', 'asn', 'glu', 'gln', 'aad'],  # arbitrary
                                                      ['val', 'leu', 'ile', 'abu', 'iva'],
                                                      "leu", ["ile"], False, 5)

    def test_valid(self):
        pred = self.pred

        assert isinstance(pred, data_structures.Prediction)
        assert pred.angstrom_code == "L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW"
        assert pred.physicochemical_class == "hydrophobic-aliphatic"
        assert pred.large_cluster_pred == ['asp', 'asn', 'glu', 'gln', 'aad']
        assert pred.small_cluster_pred == ['val', 'leu', 'ile', 'abu', 'iva']
        assert pred.single_amino_pred == "leu"
        assert pred.stachelhaus_predictions == ["ile"]
        assert pred.uncertain is False
        assert pred.stachelhaus_match_count == 5

    def test_bad_line(self):
        line = (
            "nrpspksdomains_bpsD_A1	DAATLAAVAK	"
            "hydrophobic-aromatic	phe,trp,phg,tyr,bht	tyr,bht	tyr	tyr	"
            "phe,trp,phg,tyr,bht	tyr,bht	0	0:0	0.000000e+00"
        )
        with self.assertRaisesRegex(ValueError, "Invalid SVM result line:"):
            nrps_predictor.PredictorSVMResult.from_line(line)

    def test_stach_scoring(self):
        for i in range(10):
            stach = "X" * i + "DAFYLGMMCK"[i:]
            line = (
                "nrpspksdomains_bpsA_AMP-binding.1	L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW	"
                "%s	hydrophobic-aliphatic	N/A	val,leu,ile,abu,iva	leu	leu	"
                "orn,lys,arg	gly,ala	0	0:0	0.000000e+00"
            ) % stach
            pred = nrps_predictor.PredictorSVMResult.from_line(line)
            assert pred.stachelhaus_match_count == 10 - i

    def test_classification_certain(self):
        pred = self.pred
        assert not pred.uncertain

        assert pred.get_classification() == [pred.single_amino_pred]

        pred.single_amino_pred = "N/A"
        assert pred.get_classification() == pred.small_cluster_pred

        pred.small_cluster_pred = ["N/A"]
        assert pred.get_classification() == pred.large_cluster_pred

        pred.large_cluster_pred = ["N/A"]
        assert pred.get_classification() == [pred.physicochemical_class]

    def test_classification_uncertain(self):
        pred = self.pred
        pred.uncertain = True

        assert pred.single_amino_pred != "N/A"
        assert pred.get_classification() == []

        pred.single_amino_pred = "N/A"
        assert pred.small_cluster_pred != ["N/A"]
        assert pred.get_classification() == []

        pred.small_cluster_pred = ["N/A"]
        assert pred.large_cluster_pred != ["N/A"]
        assert pred.get_classification() == []

        pred.large_cluster_pred = ["N/A"]
        assert pred.get_classification() == []

        for i in range(10):
            pred.stachelhaus_match_count = i
            if i < 8:
                assert pred.get_classification() == []
            else:
                assert pred.get_classification() == pred.stachelhaus_predictions

    def test_classification_stach_and_single(self):
        assert not self.pred.uncertain

        # stach 10 wins
        self.pred.stachelhaus_match_count = 10
        assert self.pred.get_classification() == ["ile"]

        # stach 8 loses
        self.pred.stachelhaus_match_count = 8
        assert self.pred.get_classification() == ["leu"]

        self.pred.stachelhaus_match_count = 9
        # mismatch, so none
        assert self.pred.single_amino_pred not in self.pred.stachelhaus_predictions
        assert self.pred.get_classification() == []
        # single in stach, result is intersection
        self.pred.single_amino_pred = "ile"
        self.pred.stachelhaus_predictions = ["ile", "leu"]
        assert self.pred.get_classification() == ["ile"]

    def test_classification_stach_and_groups(self):
        assert not self.pred.uncertain
        self.pred.single_amino_pred = "N/A"
        self.pred.physicochemical_class = "N/A"

        old_small = self.pred.small_cluster_pred
        old_large = self.pred.large_cluster_pred

        assert self.pred.stachelhaus_predictions[0] in old_small
        assert old_large != ["N/A"] and self.pred.stachelhaus_predictions[0] not in old_large

        for stach_count in [8, 9]:
            self.pred.stachelhaus_match_count = stach_count
            self.pred.small_cluster_pred = old_small
            assert self.pred.get_classification() == ["ile"]  # overlap with small

            self.pred.small_cluster_pred = ["N/A"]
            self.pred.large_cluster_pred = old_large
            assert self.pred.get_classification() == []  # no overlap with large, but predicted

            self.pred.large_cluster_pred = ["N/A"]
            assert self.pred.get_classification() == ["ile"]  # no pred from SVM at all

            # and check there's no classification when no SVM pred and stach < 8
            self.pred.stachelhaus_match_count = 7
            assert self.pred.get_classification() == []

    def test_json(self):
        json = self.pred.to_json()
        new = nrps_predictor.PredictorSVMResult.from_json(json)
        assert vars(self.pred) == vars(new)

    def test_html(self):
        self.pred.uncertain = False
        html = self.pred.as_html()
        assert isinstance(html, Markup)
        assert "outside applicability" not in str(html)

        self.pred.uncertain = True
        html = self.pred.as_html()
        assert "outside applicability" in str(html)


class TestAngstromGeneration(unittest.TestCase):
    def setUp(self):
        self.aligns = fasta.read_fasta(path.get_full_path(__file__, 'data', 'nrpspred_aligns.fasta'))
        mock("subprocessing.run_muscle_single", returns=self.aligns)

    def tearDown(self):
        restore()

    def test_angstrom(self):
        domain = AntismashDomain(FeatureLocation(1, 2), "test")
        domain.domain_id = "query"
        domain.translation = self.aligns[domain.domain_id].replace("-", "")

        sig = nrps_predictor.get_34_aa_signature(domain)
        assert sig == "L--SFDASLFEMYLLTGGDRNMYGPTEATMCATW"


class TestHelpers(unittest.TestCase):
    def test_good_sequence(self):
        check = nrps_predictor.verify_good_sequence

        assert check("LSFDASLFEMYLLTGGDRNMYGPTEATMCATW")
        assert not check("!LSFD")
        assert not check("LSFD-LL")
        assert not check("{LS")
