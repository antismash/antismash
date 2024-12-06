# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from unittest import TestCase
from unittest.mock import Mock, patch

from antismash.common import json
from antismash.common.layers import OptionsLayer, RecordLayer, RegionLayer
from antismash.common.test.helpers import DummyHMMResult, get_simple_options
from antismash.common.secmet.test.helpers import DummyRecord, DummyRegion
from antismash.modules.terpene.data_loader import (
    CompoundGroup,
    MissingCompoundError,
    MissingHmmError,
    Reaction,
    TerpeneHMM,
    load_hmm_properties,
)
from antismash.modules.terpene.results import (
    DomainPrediction,
    ProtoclusterPrediction,
    TerpeneResults,
)
from antismash.modules import terpene
from antismash.modules.terpene import (
    terpene_analysis,
    html_output,
)

_hmm_properties = load_hmm_properties()


class DummyCompoundGroup(CompoundGroup):
    def __init__(self,
                 name="Dummy",
                 extended_name="Long dummy",
                 single_compound=False,
                 biosynthetic_class="diterpene",
                 biosynthetic_subclass=None,
                 chain_length=10,
                 initial_cyclisations=("unknown",),
                 functional_groups=("PP",)):
        super().__init__(name, extended_name, single_compound, biosynthetic_class,
                         biosynthetic_subclass, chain_length, initial_cyclisations, functional_groups)


class DummyDomainPrediction(DomainPrediction):
    def __init__(self, domain_type="T1TS", subtypes=None, start=1, end=200, reactions=None):
        super().__init__(domain_type, subtypes or tuple(), start, end, reactions or tuple())


class DummyTerpeneHMM(TerpeneHMM):
    def __init__(self,
                 name="PT_FPP_bact",
                 description="Prenyltransferase; Farnesyl diphosphate synthase, bacterial",
                 domain_type="PT_FPPS_like",
                 length=265,
                 cutoff=250,
                 subtypes=tuple(),
                 reactions=tuple()):
        super().__init__(name, description, type, length, cutoff, subtypes, reactions)


EXISTING_COMPOUND_GROUPS = {
    "GFPP": CompoundGroup(
        name="GFPP",
        extended_name="geranylfarnesyl diphosphate",
        single_compound=True,
        biosynthetic_class="sesterterpene precursor",
        biosynthetic_subclass="",
        chain_length=25,
        initial_cyclisations=tuple(),
        functional_groups=("PP",),
    ),
    "C40_carotenoid": CompoundGroup(
        name="C40_carotenoid",
        extended_name=None,
        single_compound=False,
        biosynthetic_class="tetraterpene",
        biosynthetic_subclass="C40 carotenoid",
        chain_length=40,
        initial_cyclisations=("unknown",),
        functional_groups=("unknown",),
    )
}


def build_dummy_reaction(substrates=(EXISTING_COMPOUND_GROUPS["GFPP"],),
                         products=(EXISTING_COMPOUND_GROUPS["C40_carotenoid"],)):
    return Reaction(substrates=substrates, products=products)


def build_dummy_domain():
    return DomainPrediction(domain_type="T1TS", subtypes=("T1TS_Bas_a",),
                            start=1, end=200, reactions=(build_dummy_reaction(),))


def build_dummy_cds_preds():
    return {
        "cds1": [build_dummy_domain()],
        "cds2": [],
    }


def build_dummy_cluster_pred():
    return ProtoclusterPrediction(cds_predictions=build_dummy_cds_preds(),)


def build_dummy_terpene_hmms():
    return {
        "PT_FPP_bact": TerpeneHMM(
            name="PT_FPP_bact",
            description="Prenyltransferase; Farnesyl diphosphate synthase, bacterial",
            domain_type="PT_FPPS_like",
            length=265,
            cutoff=250,
            subtypes=tuple(),
            reactions=(build_dummy_reaction(),),
        )
    }


FAKE_REACTION_DATA = {"substrates": ["Fake1", "Fake2"], "products": ["Fake3"]}

FAKE_HMM_DATA = {
    "name": "PT_FPP_bact",
    "description": "Prenyltransferase; Farnesyl diphosphate synthase, bacterial",
    "length": 265,
    "cutoff": 250,
    "subtypes": ["Fake_subtype"],
    "reactions": [FAKE_REACTION_DATA],
}

FAKE_CLUSTER_DATA = {"cds_predictions": {},
                     "products": ["Fake1", "Fake2"]}


class TestJSONConversion(TestCase):
    def test_compound_regeneration(self):
        compound = EXISTING_COMPOUND_GROUPS["GFPP"]
        regenerated = CompoundGroup.from_json(json.loads(json.dumps(compound.to_json())))
        assert regenerated == compound

    def test_reaction_regeneration(self):
        reaction = build_dummy_reaction()
        regenerated = Reaction.from_json(json.loads(json.dumps(reaction.to_json())),
                                         EXISTING_COMPOUND_GROUPS)
        assert regenerated == reaction

    def test_domain_pred_regeneration(self):
        pred = build_dummy_domain()
        regenerated = DomainPrediction.from_json(json.loads(json.dumps(pred.to_json())),
                                                 EXISTING_COMPOUND_GROUPS)
        assert regenerated == pred

    def test_cluster_pred_regeneration(self):
        pred = build_dummy_cluster_pred()
        pred.add_product(EXISTING_COMPOUND_GROUPS["GFPP"])
        raw = json.loads(json.dumps(pred.to_json()))
        regenerated = ProtoclusterPrediction.from_json(raw, EXISTING_COMPOUND_GROUPS)
        assert regenerated.to_json() == pred.to_json()
        assert regenerated == pred

    def test_regeneration(self):
        record = DummyRecord()
        results = TerpeneResults(record.id)
        results.cluster_predictions = {1: build_dummy_cluster_pred()}
        regenerated = TerpeneResults.from_json(json.loads(json.dumps(results.to_json())), record)
        assert regenerated.cluster_predictions == results.cluster_predictions

    def test_missing_compounds(self):
        with self.assertRaises(MissingCompoundError):
            Reaction.from_json(FAKE_REACTION_DATA, EXISTING_COMPOUND_GROUPS)
        with self.assertRaises(MissingCompoundError):
            ProtoclusterPrediction.from_json(FAKE_CLUSTER_DATA, EXISTING_COMPOUND_GROUPS)
        with self.assertRaises(MissingHmmError):
            TerpeneHMM.from_json(FAKE_HMM_DATA, build_dummy_terpene_hmms(), EXISTING_COMPOUND_GROUPS)

    def test_missing_compounds_top_level(self):
        record = DummyRecord()
        data = {
            "record_id": record.id,
            "schema_version": TerpeneResults._schema_version,
            "protocluster_predictions": {1: build_dummy_cluster_pred().to_json()}
        }
        # changed compounds should discard existing results
        with patch.object(terpene.results, "load_compounds", return_value={}):
            assert TerpeneResults.from_json(data, record) is None


class TestDatatypes(TestCase):
    def test_hmm_substrates(self):
        reaction1 = build_dummy_reaction(substrates=(DummyCompoundGroup(name="compound1"),
                                                     DummyCompoundGroup(name="compound2"),))
        reaction2 = build_dummy_reaction(substrates=(DummyCompoundGroup(name="compound1"),))
        with self.assertRaises(ValueError):
            DummyTerpeneHMM(reactions=(reaction1, reaction2,))


class TestAnalysis(TestCase):
    def test_filter_by_score(self):
        hmm_results = [DummyHMMResult(label="T1TS_III-IV_b", bitscore=300),
                       DummyHMMResult(label="T2TS_squalene", bitscore=300)]
        hmm_results_by_id = {"cds1": [hmm_results[0]],
                             "cds2": [hmm_results[1]]}
        filtered = terpene_analysis.filter_by_score(hmm_results_by_id, _hmm_properties)
        expected_filtered = {"cds1": [hmm_results[0]]}
        assert filtered == expected_filtered

    def test_group_hmms(self):
        hmm_results = [DummyHMMResult(label="independent_dom", start=250, end=400),
                       DummyHMMResult(label="PT_GGPP", start=20, end=200),
                       DummyHMMResult(label="PT_noFPP_bact", start=1, end=100)]
        groups = terpene_analysis.group_hmm_results(hmm_results)
        expected_groups = [[hmm_results[2], hmm_results[1]], [hmm_results[0]]]
        assert groups == expected_groups

    def test_merge_reactions_by_substrate(self):
        reaction1 = build_dummy_reaction(
            substrates=(DummyCompoundGroup(name="compound1"), DummyCompoundGroup(name="compound2"),),
            products=(DummyCompoundGroup(name="compound3"), DummyCompoundGroup(name="compound4"),)
        )
        reaction2 = build_dummy_reaction(
            substrates=(DummyCompoundGroup(name="compound1"), DummyCompoundGroup(name="compound2"),),
            products=(DummyCompoundGroup(name="compound3"),)
        )
        profiles = [
            DummyTerpeneHMM(reactions=(reaction1,)),
            DummyTerpeneHMM(reactions=(reaction2,)),
            DummyTerpeneHMM(reactions=tuple()),
        ]
        merged = terpene_analysis.merge_reactions_by_substrate(profiles)[0]

        def names(reaction):
            return set(s.name for s in reaction.substrates)

        assert names(merged) == names(reaction1) == names(reaction2)
        assert merged.products == reaction2.products

    def test_get_domain_prediction(self):
        hmm_results = [DummyHMMResult(label="T1TS", start=1, end=50),
                       DummyHMMResult(label="T1TS_ARIS", start=30, end=100)]
        domain_pred = terpene_analysis.get_domain_prediction(hmm_results, _hmm_properties)

        assert domain_pred.domain_type == "T1TS"
        assert domain_pred.subtypes == ("T1TS_ARIS",)
        assert domain_pred.start == 1
        assert domain_pred.end == 100
        assert domain_pred.reactions[0].products[0].name == "aristolochene"

        hmm_results2 = [DummyHMMResult(label="T1TS", start=1, end=50),
                        DummyHMMResult(label="T2TS", start=30, end=100)]
        domain_pred = terpene_analysis.get_domain_prediction(hmm_results2, _hmm_properties)
        assert domain_pred.domain_type == "ambiguous"
        assert domain_pred.subtypes == tuple()

    def test_get_cds_predictions(self):
        hmm_results_per_cds = {"cds1": [DummyHMMResult(label="T1TS", bitscore=300)]}
        cds_preds = terpene_analysis.get_cds_predictions(hmm_results_per_cds, _hmm_properties)
        assert cds_preds == {"cds1": [DomainPrediction(domain_type="T1TS", subtypes=tuple(),
                                                       start=1, end=10, reactions=tuple())]}

    def test_get_cluster_prediction(self):
        bad_product = DummyCompoundGroup(name="BadProduct")
        good_product1 = DummyCompoundGroup(name="GoodProduct1")
        good_product2 = DummyCompoundGroup(name="GoodProduct2")
        good_product3 = DummyCompoundGroup(name="GoodProduct3")
        cds_preds = {
            "cds1": [
                DummyDomainPrediction(domain_type="ambiguous",
                                      reactions=tuple()),
                DummyDomainPrediction(domain_type="T1TS",
                                      reactions=(Mock(products=(good_product1, good_product2,)),)),
                DummyDomainPrediction(domain_type="T2TS",
                                      reactions=(Mock(products=(bad_product,)),))
            ],
            "cds2": [
                DummyDomainPrediction(domain_type="T1TS",
                                      reactions=(Mock(products=(good_product1, good_product3,)),)),
                DummyDomainPrediction(domain_type="Lycopene_cycl",
                                      reactions=tuple())
            ],
        }
        cluster_pred = terpene_analysis.get_cluster_prediction(cds_preds)
        assert cluster_pred.products == [good_product1, good_product2, good_product3]


class TestHTML(TestCase):
    def test_get_domain_description(self):
        domain_pred = DummyDomainPrediction(domain_type="T1TS")
        description = html_output.get_domain_description(domain_pred)
        expected_description = "Type 1 terpene synthase"
        assert description == expected_description

    def test_format_subtype(self):
        domain_pred = DummyDomainPrediction(subtypes=("T1TS_Bas_a", "T1TS_Bas_b"))
        assert html_output.format_subtype(domain_pred) == "Basidiomycete, group a or Basidiomycete, group b"
        domain_pred = DummyDomainPrediction(domain_type="T1TS", subtypes=tuple())
        assert html_output.format_subtype(domain_pred) == "unknown"
        domain_pred = DummyDomainPrediction(domain_type="TS_UbiA", subtypes=tuple())
        assert html_output.format_subtype(domain_pred) == "none"

    def test_format_reactions(self):
        reaction1 = build_dummy_reaction(substrates=(DummyCompoundGroup(name="substrate"),),
                                         products=(DummyCompoundGroup(name="product1"),
                                                   DummyCompoundGroup(name="product2"),))
        domain_pred = DummyDomainPrediction(reactions=(reaction1,))
        markup = html_output.format_reactions(domain_pred)
        expected_markup = "<dd>substrate &rarr; product1</dd><dd>substrate &rarr; product2</dd>"
        assert markup == expected_markup

    def test_format_domain_types(self):
        domain_preds = [DummyDomainPrediction(domain_type="T1TS"),
                        DummyDomainPrediction(domain_type="T2TS")]
        assert html_output.format_domain_types(domain_preds) == "T1TS + T2TS"

    def test_get_preds_by_cluster(self):
        region = DummyRegion()
        record = DummyRecord(features=[region])
        results = TerpeneResults(record.id)
        options_layer = OptionsLayer(get_simple_options(terpene, []), modules=[terpene])
        region_layer = RegionLayer(RecordLayer(record=record, results=results, options=options_layer), region)
        clusters = []
        for i, category in enumerate(["terpene", "terpene", "NRPS"]):
            clusters.append(Mock(get_protocluster_number=lambda: i, product_category=category))
        results.cluster_predictions = {i: build_dummy_cluster_pred() for i in range(1, 3)}
        with patch.object(region_layer, "get_unique_protoclusters", return_value=clusters):
            preds_by_cluster = html_output.get_preds_by_cluster(region_layer, results)
        assert len(preds_by_cluster) == 2
        assert clusters[0] in preds_by_cluster
        assert clusters[1] in preds_by_cluster

    def test_get_glossary_data(self):
        predictions = [build_dummy_cluster_pred()]
        name_mappings = html_output.get_glossary_data(predictions)
        assert name_mappings == {'GFPP': 'Geranylfarnesyl diphosphate'}
