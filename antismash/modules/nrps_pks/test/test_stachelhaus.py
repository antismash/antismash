# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import os
import tempfile
import unittest
from unittest.mock import patch

import antismash
from antismash.config import build_config
from antismash.common.secmet.test.helpers import DummyAntismashDomain
from antismash.modules.nrps_pks import stachelhaus
from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name


class TestStachSignature(unittest.TestCase):
    def test_json_roundtrip(self):
        short_name = "Tyr"
        substrate = get_substrate_by_name(short_name)
        sig = stachelhaus.StachSignature(
            "FAKEDATAOK", "ILIKEDATAEVENFAKEDATADIDISAYILIKED", {substrate}, short_name)
        data = json.dumps(sig.to_json())
        print(data)
        restored = stachelhaus.StachSignature.from_json(json.loads(data))
        assert sig == restored


class TestStachelhausPrediction(unittest.TestCase):
    def setUp(self) -> None:
        short_name = "Tyr"
        self.substrate = get_substrate_by_name(short_name)
        sig = stachelhaus.StachSignature(
            "FAKEDATAOK", "ILIKEDATAEVENFAKEDATADIDISAYILIKED", {self.substrate}, short_name)
        self.prediction = stachelhaus.StachelhausPrediction(sig)

    def test_json_roundtrip(self):
        nothing = stachelhaus.StachelhausPrediction(original_signature="FAKEDATAOK")
        data = json.dumps(nothing.to_json())
        restored = stachelhaus.StachelhausPrediction.from_json(json.loads(data))
        assert restored.get_classification() == []
        assert restored.aa10 == "FAKEDATAOK"


        data = json.dumps(self.prediction.to_json())
        restored = stachelhaus.StachelhausPrediction.from_json(json.loads(data))
        assert restored.get_classification() == self.prediction.get_classification()

        with self.assertRaisesRegex(ValueError, "provide either a StachSignature or an original 10 AA signature"):
            stachelhaus.StachelhausPrediction()

    def test_get_classification(self):
        assert [] == stachelhaus.StachelhausPrediction(original_signature="FAKEDATAOK").get_classification()

        assert [self.substrate.short] == self.prediction.get_classification()

    def test_html(self):
        nothing = stachelhaus.StachelhausPrediction(original_signature="FAKEDATAOK")
        assert "No exact match" in str(nothing.as_html())

        assert self.substrate.long in str(self.prediction.as_html())


FAKE_A_DOMAIN_SIGNATURES = (
    ("FAKEDATAOK", "ILIKEDATAEVENFAKEDATADIDISAYILIKED"),  # shouldn't hit
    ("DVSAIGCVTK", "LDLIFDVSVSEMAMIVGGEINCYGPTETTVTATL"),  # should hit once, Trp
    ("DAACIGNVVK", "MWHVFDASIAEPCLITGGDYNNYGPTENTVVATW"),  # should hit two 10 AA sigs but only one 34 AA sig, Trp
    ("DAACIGNVVK", "LWHVFDASIAEPCLITGGDYNNYGPTENTVVATW"),  # shouldn't hit any 34 AA sig, but consensus Trp anyway
    ("DFWNIGMVHK", "LALAFDFSVWEGNQIFGGEVNMYGITETTVHVSY"),  # no aa34 hits, consensus Thr with a more complicated consensus
    ("DVMFYTALVK", "LRQMFDVSFMECFLYTTGEINAYGPSEAHLVSAR"),  # multiple substrates in consensus results, Trp|Tyr
    (None, None),                                          # Should not break on this
)

FAKE_KNOWN_STACH_CODES = {
    "FAKEDATAOK": [
        stachelhaus.StachSignature(
            "FAKEDATAOK", "FAKEDATAOKFAKEDATAOKFAKEDATAOKFAKE", set([stachelhaus.get_substrate_by_name("Val")]), "Val"),
        stachelhaus.StachSignature(
            "FAKEDATAOK", "FAKEFAKEDATAOKFAKEDATAOKFAKEDATAOK", set([stachelhaus.get_substrate_by_name("Ile")]), "Ile"),
    ],
}


class TestStachelhaus(unittest.TestCase):
    def setUp(self) -> None:
        self.config = build_config([], isolated=True, modules=antismash.get_all_modules())
        self.domains = []
        for domain_id in [
                "no_hit", "one_hit", "aa34_hit", "simple_consensus", "tricky_consensus", "multi",
                "none"]:
            self.domains.append(DummyAntismashDomain(locus_tag="A", domain_id=domain_id))


    @patch.object(stachelhaus, "get_a_dom_signatures", side_effect=FAKE_A_DOMAIN_SIGNATURES)
    def test_run_stachelhaus(self, _patched_get_a_dom_signatures):
        predictions = stachelhaus.run_stachelhaus(self.domains, self.config)
        expected_results = {
            "no_hit": [],
            "one_hit": ["Trp"],
            "aa34_hit": ["Trp"],
            "simple_consensus": ["Trp"],
            "tricky_consensus": ["Thr"],
            "multi": ["Trp", "Tyr"],
        }
        for name, expected in expected_results.items():
            assert expected == predictions[name].get_classification()

        assert "none" not in predictions

    @patch.object(stachelhaus, "init_data", return_value=FAKE_KNOWN_STACH_CODES)
    @patch.object(stachelhaus, "get_a_dom_signatures", side_effect=[("FAKEDATAOK", "ILIKEDATAEVENFAKEDATADIDISAYILIKED")])
    def test_run_stachelhaus_combine_consensus(self, _patched_sigs, _patched_codes):
        domains = [
            DummyAntismashDomain(locus_tag="A", domain_id="different_multi"),
        ]
        predictions = stachelhaus.run_stachelhaus(domains, self.config)

        assert "different_multi" in predictions
        assert predictions["different_multi"].get_classification() == ["Val", "Ile"]

    def test_init_data(self):
        stachelhaus.KNOWN_STACH_CODES.clear()
        assert not stachelhaus.KNOWN_STACH_CODES
        mappings = stachelhaus.init_data(self.config)
        assert mappings == stachelhaus.KNOWN_STACH_CODES

    def test_check_prereqs(self):
        with tempfile.TemporaryDirectory(prefix="aS.stachelhaustest") as temp_dbdir:
            dir_name = os.path.join(temp_dbdir, "nrps_pks", "stachelhaus", "1.0")
            os.makedirs(dir_name)
            with open(os.path.join(dir_name, "dummy.file"), "w", encoding="utf-8") as handle:
                handle.write("dummy content")
            config = build_config(["--databases", temp_dbdir], isolated=True, modules=antismash.get_all_modules())
            assert stachelhaus.check_prereqs(config)
