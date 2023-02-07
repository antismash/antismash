# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest
from unittest.mock import patch

from helperlibs.wrappers.io import TemporaryDirectory

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import nrps_pks


class IntegrationWithoutNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal"], isolated=True, modules=antismash.get_all_modules())
        assert not nrps_pks.is_enabled(self.options)

    def tearDown(self):
        destroy_config()

    def test_minimal(self):
        with TemporaryDirectory(change=True) as tempdir:
            self.options = build_config(["--minimal", "--output-dir", tempdir],
                                        isolated=True, modules=antismash.get_all_modules())
            with patch.object(nrps_pks, "run_on_record", side_effect=RuntimeError("shouldn't run")):
                antismash.main.run_antismash(helpers.get_path_to_balhymicin_genbank(),
                                             self.options)


def ensure_module_lids_shown(output_dir: str):
    starter = '<input class="show-module-domains" type="checkbox"'
    shown_snippet = f"{starter}>"
    hidden_snippet = f"{starter} checked>"
    with open(os.path.join(output_dir, "index.html"), encoding="utf-8") as handle:
        content = handle.read()
        assert starter in content, f"all variants missing of {starter}"
        shown = shown_snippet in content
        hidden = hidden_snippet in content
        # one of the above has to be true, if both are missing that's a big (and different) problem
        assert not hidden, f"'{hidden_snippet}' found in HTML, lids would be hidden"
        assert shown, f"'{shown_snippet}' missing from HTML, lids wouldn't be shown"


class IntegrationNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-nrps-pks", "--enable-html"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_balhymicin(self):
        filename = helpers.get_path_to_balhymicin_genbank()
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options,
                                                                callback=ensure_module_lids_shown)
        assert len(results.domain_predictions) == 9
        a_domains = [("bpsA", 1), ("bpsA", 2), ("bpsA", 3),
                     ("bpsB", 1), ("bpsB", 2), ("bpsB", 3),
                     ("bpsC", 1),
                     ("bpsD", 1)]
        nrps_names = [f"nrpspksdomains_{name}_AMP-binding.{index}" for name, index in a_domains]
        feature_names = nrps_names + ["nrpspksdomains_pks_CAL_domain.1"]
        assert set(results.domain_predictions) == set(feature_names)

        assert set(results.domain_predictions[feature_names[0]]) == {"NRPSPredictor2", "Stachelhaus"}
        nrpspred2_results = {}
        nrpspred2_norine_results = {}
        stachelhaus_results = {}
        stachelhaus_norine_results = {}
        for domain, methods in results.domain_predictions.items():
            if "CAL" in domain:
                continue
            nrpspred2_results[domain] = methods["NRPSPredictor2"].get_classification()
            nrpspred2_norine_results[domain] = methods["NRPSPredictor2"].get_classification(as_norine=True)
            stachelhaus_results[domain] = methods["Stachelhaus"].get_classification()
            stachelhaus_norine_results[domain] = methods["Stachelhaus"].get_classification(as_norine=True)

        expected_preds = [["leu"], ["bht"], ["asn"], ["hpg"], ["hpg"], ["bht"], ["dhpg"], ["tyr"], ["pk"]]
        expected_stachelhaus_preds = [["Leu"], ["R-ohTyr"], ["Asn"], ["Hpg"], ["Hpg"], ["R-ohTyr"], ["dHpg"], ["Tyr"]]
        expected_norine_preds = [["Leu"], ["bOH-Tyr"], ["Asn"], ["Hpg"], ["Hpg"], ["bOH-Tyr"], ["Dhpg"], ["Tyr"]]

        expected_nrps2 = dict(zip(nrps_names, expected_preds))
        assert nrpspred2_results == expected_nrps2
        expected_nrps2_norine = dict(zip(nrps_names, expected_norine_preds))
        assert nrpspred2_norine_results == expected_nrps2_norine
        expected_stachelhaus = dict(zip(nrps_names, expected_stachelhaus_preds))
        assert stachelhaus_results == expected_stachelhaus
        expected_stachelhaus_norine = dict(zip(nrps_names, expected_norine_preds))
        assert stachelhaus_norine_results == expected_stachelhaus_norine

        cal = results.domain_predictions["nrpspksdomains_pks_CAL_domain.1"]["minowa_cal"]
        assert len(cal.predictions) == 5
        assert cal.predictions[0] == ["AHBA", 167.0]

        assert len(results.region_predictions[1]) == 4
        # as does this, though it still won't use domain docking
        pred = results.region_predictions[1][1]
        monomers = '(Leu - D-bOH-Tyr - Asn) + (D-Hpg - D-Hpg - bOH-Tyr) + (Dhpg) + (Tyr)'
        assert pred.polymer == monomers
        assert not pred.domain_docking_used
        assert pred.ordering == ['bpsA', 'bpsB', 'bpsC', 'bpsD']

        pred = results.region_predictions[1][2]
        assert pred.polymer == "(Tyr) + (Dhpg)"
        assert not pred.domain_docking_used
        assert pred.ordering == ['bpsD', 'bpsC']

    def test_cp002271_c19(self):
        filename = path.get_full_path(__file__, 'data', 'CP002271.1.cluster019.gbk')
        results = helpers.run_and_regenerate_results_for_module(filename, nrps_pks, self.options,
                                                                callback=ensure_module_lids_shown)
        # catch ordering changes along with ensuring ATResults are there
        pred = results.domain_predictions["nrpspksdomains_STAUR_3982_PKS_AT.1"]
        assert pred["signature"].predictions[0][1].score == 87.5
        # ensure all genes are present and have the right consensus
        assert results.consensus == {'nrpspksdomains_STAUR_3982_PKS_AT.1': 'mmal',
                                     'nrpspksdomains_STAUR_3983_PKS_AT.1': 'mmal',
                                     'nrpspksdomains_STAUR_3984_PKS_AT.1': 'mmal',
                                     'nrpspksdomains_STAUR_3985_PKS_AT.1': 'pk',
                                     'nrpspksdomains_STAUR_3985_PKS_AT.2': 'mmal'}
        assert len(results.region_predictions) == 1
        assert list(results.region_predictions) == [1]
        assert len(results.region_predictions[1]) == 1
        # check the gene ordering and, in this case, that it used domain docking
        sc_pred = results.region_predictions[1][0]
        assert sc_pred.polymer == '(Me-ccmal) + (Me-ccmal) + (Me-ccmal)'
        assert sc_pred.domain_docking_used
        assert sc_pred.ordering == ['STAUR_3984', 'STAUR_3983', 'STAUR_3982']
        assert len(results.domain_predictions) == 10
        expected_domains = {'nrpspksdomains_STAUR_3982_PKS_AT.1',
                            'nrpspksdomains_STAUR_3983_PKS_AT.1',
                            'nrpspksdomains_STAUR_3984_PKS_AT.1',
                            'nrpspksdomains_STAUR_3985_PKS_AT.1',
                            'nrpspksdomains_STAUR_3985_PKS_AT.2',
                            'nrpspksdomains_STAUR_3972_PKS_KR.1',
                            'nrpspksdomains_STAUR_3984_PKS_KR.1',
                            'nrpspksdomains_STAUR_3985_PKS_KR.1',
                            'nrpspksdomains_STAUR_3983_PKS_KR.1',
                            'nrpspksdomains_STAUR_3983_PKS_KR.1',
                            'nrpspksdomains_STAUR_3982_PKS_KR.1'}
        assert set(results.domain_predictions) == expected_domains

    def test_get_a_dom_signatures(self):
        filename = path.get_full_path(__file__, 'data', 'dom_signatures.txt')
        data = []
        for line in open(filename, 'r', encoding="utf-8"):
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            data.append(line.split("\t"))

        for name, seq, aa34, aa10 in data:
            domain = helpers.DummyAntismashDomain(locus_tag="A", domain_id=name)
            domain._translation = seq
            assert (aa10, aa34) == nrps_pks.signatures.get_a_dom_signatures(domain), name
