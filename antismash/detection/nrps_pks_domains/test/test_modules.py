# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,too-many-public-methods

import json
import unittest

from antismash.common import path
from antismash.common.test.helpers import DummyHMMResult
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.detection import nrps_pks_domains
from antismash.detection.nrps_pks_domains.module_identification import (
    CDSModuleInfo,
    CLASSIFICATIONS,
    Component as Component_actual,
    DOUBLE_TRANSPORTER_CASES,
    Module,
    build_modules_for_cds as build_modules_for_cds_actual,
    classify,
    combine_modules,
)

# chosen arbitrarily, these exist to make future profile renames easier
NRPS_START = "Cglyc"
NRPS_LOAD = "AMP-binding"
PKS_START = "PKS_KS"
PKS_LOAD = "PKS_AT"
TRANS_AT_SUBTYPE = "Trans-AT-KS"
CP = "ACP"


class Component(Component_actual):
    """ a tiny wrapper to avoid always supplying a dummy CDS name """
    def __init__(self, domain, cds_name="test_name", subtype=None):
        if subtype:
            sub = DummyHMMResult(subtype)
            domain.add_internal_hits([sub])
        super().__init__(domain, cds_name)


def build_modules_for_cds(domains, subtypes, cds_name="test_name"):
    subs = iter(subtypes)
    for domain in domains:
        if domain.hit_id == PKS_START:
            sub = DummyHMMResult(next(subs), domain.query_start, domain.query_end)
            domain.add_internal_hits([sub])
    return build_modules_for_cds_actual(domains, cds_name)


def add_component(module, name, sub=None, start=1, end=10, cds_name="test_name", lookaheads=None):
    assert cds_name
    if lookaheads is None:
        lookaheads = []
    module.add_component(Component(DummyHMMResult(name, start, end), cds_name, sub), lookaheads)


def build_module(names, subtypes=None, first_in_cds=True, cds_name="test_name"):
    module = Module(first_in_cds=first_in_cds)
    subs = iter(subtypes or [])
    for domain in names:
        sub = next(subs) if domain == PKS_START and subtypes else None
        add_component(module, domain, sub, cds_name=cds_name)
    return module


class TestClassify(unittest.TestCase):
    def test_existing(self):
        assert CLASSIFICATIONS
        for classification, group in CLASSIFICATIONS.items():
            for label in group:
                assert classify(label) == classification

    def test_all_domains_classified(self):
        domain_names = []
        with open(path.get_full_path(nrps_pks_domains.__file__, "data", "nrpspksdomains.hmm"),
                  encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("NAME"):
                    domain_names.append(line.strip().split()[1])
        missing = [name for name in domain_names if not classify(name)]
        assert not missing, missing

    def test_unclassifiable(self):
        with self.assertRaisesRegex(ValueError, "could not classify"):
            classify("bad-domain-name")


class TestComponent(unittest.TestCase):
    def test_construction(self):
        domain = DummyHMMResult("ACP")
        component = Component(domain)
        assert component.domain == domain
        assert component.label == domain.hit_id
        assert component.subtype is None
        assert component.classification == "CP"

        domain._hit_id = PKS_START
        component = Component(domain, subtype="some-subtype")
        assert component.subtype == "some-subtype"
        assert component.subtypes == ["some-subtype"]
        assert component.classification == "KS"

    def test_json(self):
        component = Component(DummyHMMResult(PKS_START), subtype=TRANS_AT_SUBTYPE)
        intermediate = component.to_json()
        new = Component.from_json(json.loads(json.dumps(intermediate)))
        assert new.to_json() == intermediate
        assert new.domain == component.domain
        assert new.subtypes == component.subtypes
        assert new.classification == component.classification

    def test_condensation(self):
        for cond in CLASSIFICATIONS["C"]:
            assert Component(DummyHMMResult(cond)).is_condensation()
        assert not Component(DummyHMMResult(NRPS_LOAD)).is_condensation()

    def test_loader(self):
        for cond in [PKS_LOAD, "AMP-binding", "A-OX"]:
            assert Component(DummyHMMResult(cond)).is_loader()
        assert not Component(DummyHMMResult(NRPS_START)).is_loader()


class TestModule(unittest.TestCase):
    def setUp(self):
        self.pks = Module()
        for i, domain in enumerate([PKS_START, PKS_LOAD, CP]):
            start = (i + 1) * 10
            add_component(self.pks, domain, start=start, end=start + 5)
        assert self.pks.is_pks()

        self.nrps = Module()
        for i, domain in enumerate([NRPS_START, NRPS_LOAD, CP]):
            start = (i + 1) * 10
            add_component(self.nrps, domain, start=start, end=start + 5)
        assert self.nrps.is_nrps()

    def test_json(self):
        intermediate = self.pks.to_json()
        assert intermediate
        new = Module.from_json(json.loads(json.dumps(intermediate)))
        assert new.to_json() == intermediate
        assert new._starter.label == PKS_START
        assert not new._end

    def test_methylations(self):
        pks = build_module([PKS_START, PKS_LOAD])
        # mmal is a methylated malonyl-CoA, so it should look the same
        assert pks.get_monomer("mmal") == "Me-mal"
        add_component(pks, "PKS_KR")
        assert pks.get_monomer("mmal") == "Me-ohmal"
        add_component(pks, "nMT")
        assert pks.get_monomer("mmal") == "NMe-Me-ohmal"
        assert pks.get_monomer("mal") == "NMe-ohmal"
        add_component(pks, "oMT")
        assert pks.get_monomer("mal") == "NMe-OMe-ohmal"
        assert pks.get_monomer("pk") == "NMe-OMe-?"

        nrps = build_module([NRPS_START, NRPS_LOAD])
        assert nrps.get_monomer("ala") == "ala"
        add_component(nrps, "cMT")
        assert nrps.get_monomer("ala") == "Me-ala"

        # even with an unknown, the methylation should be indicated
        assert nrps.get_monomer("X") == "Me-?"

    def test_methylation_before_load(self):
        # methylation before an A domain, currently impossible
        with self.assertRaisesRegex(ValueError, "loader after other non-starter"):
            build_module([NRPS_START, "cMT", NRPS_LOAD])

    def test_epimerase(self):
        # should always add a D- at the beginning
        nrps = build_module([NRPS_START, NRPS_LOAD, CP])
        assert nrps.get_monomer("ile") == "ile"
        add_component(nrps, "Epimerization")
        assert nrps.get_monomer("ile") == "D-ile"

        nrps = build_module([NRPS_START, NRPS_LOAD, "cMT", CP])
        assert nrps.get_monomer("gly") == "Me-gly"
        nrps = build_module([NRPS_START, NRPS_LOAD, "cMT", CP, "Epimerization"])
        assert nrps.get_monomer("gly") == "D-Me-gly"

    def test_start(self):
        assert self.pks.start == 10
        self.pks._starter.domain._query_start = 9
        assert self.pks.start == 9
        self.pks._components.pop(0)
        assert self.pks.start == self.pks._loader.domain.query_start

    def test_end(self):
        assert self.nrps.end == 35
        add_component(self.nrps, "Epimerization", start=40, end=60)
        assert self.nrps.end == 60
        self.nrps._end = None
        self.nrps._components.pop(-1)
        assert self.nrps.end == 35
        # TE/TD not included in the length, but are included in the module
        add_component(self.nrps, "Thioesterase", start=70, end=80)
        assert self.nrps.end == 35

    def test_termination(self):
        domains = [NRPS_LOAD, CP]
        module = build_module(domains)
        assert not module.is_terminated()

        add_component(module, "Epimerization")
        assert module.is_terminated()
        assert not module.is_termination_module()

        module = build_module(domains + ["Thioesterase"])
        assert module.is_terminated()
        assert module.is_termination_module()

    def test_starter_module(self):
        assert not self.pks.is_starter_module()
        for starter in ["Condensation_Starter", "CAL_domain", "SAT", NRPS_LOAD]:
            assert build_module([starter]).is_starter_module(), starter

        for domains in [[NRPS_START, NRPS_LOAD],
                        [PKS_START, PKS_LOAD],
                        [CP],
                        ]:
            assert not build_module(domains).is_starter_module(), domains

    def test_no_monomer(self):
        assert self.pks.get_monomer() == ""
        assert self.pks.get_monomer("") == ""

    def test_pk_nrp_indicated(self):
        assert self.pks.get_monomer("pk") == "?"
        assert self.pks.get_monomer("mal") != "?"
        assert self.nrps.get_monomer("X") == "?"
        assert self.nrps.get_monomer("ala") != "?"

    def test_monomer_trans_at_default(self):
        trans_at = build_module([PKS_START, CP], [TRANS_AT_SUBTYPE])
        assert trans_at._starter.subtype == TRANS_AT_SUBTYPE
        assert trans_at._loader is None
        assert trans_at.is_trans_at()
        assert trans_at.get_monomer("") == "mal"

    def test_trans_at(self):
        assert not self.pks.is_trans_at()
        assert self.pks.is_complete()

        module = build_module([PKS_START, PKS_LOAD, CP], [TRANS_AT_SUBTYPE])
        assert not module.is_trans_at()

        module = build_module([PKS_START, CP])
        assert not module.is_trans_at()
        assert not module.is_complete()

        module = build_module([PKS_START, CP], [TRANS_AT_SUBTYPE])
        assert module.is_trans_at()
        assert module.is_complete()

        # without matching subtype is fine if the matching docking domain is present
        module = build_module([PKS_START, "Trans-AT_docking", CP])
        assert module.is_trans_at()
        assert module.is_complete()

        # but the loader must still be missing
        module = build_module([PKS_START, PKS_LOAD, "Trans-AT_docking", CP])
        assert not module.is_trans_at()
        assert module.is_complete()  # still complete, just a little strange

    def test_trailing_modifiers(self):
        error = "modification domain after carrier protein"
        # not allowed for NRPSs and PKSs
        with self.assertRaisesRegex(ValueError, error):
            build_module([NRPS_START, NRPS_LOAD, CP, "nMT"])
        with self.assertRaisesRegex(ValueError, error):
            build_module([PKS_START, PKS_LOAD, CP, "nMT"])
        # except for trans-at-pks, and only KR
        assert build_module([PKS_START, CP, "PKS_KR"], [TRANS_AT_SUBTYPE]).is_trans_at()
        with self.assertRaisesRegex(ValueError, error):
            assert build_module([PKS_START, CP, "nMT"], [TRANS_AT_SUBTYPE])

    def test_pks_chaining(self):
        module = Module()
        for comp in list(self.pks)[:-1]:
            module.add_component(comp, [])
        assert module.get_monomer("mal") == "mal"
        add_component(module, "PKS_KR")
        assert module.get_monomer("mal") == "ohmal"
        add_component(module, "PKS_DH")
        assert module.get_monomer("mal") == "ccmal"
        add_component(module, "PKS_ER")
        assert module.get_monomer("mal") == "redmal"

        # and that the order doesn't matter
        module = build_module([PKS_START, PKS_LOAD, "PKS_ER", "PKS_DH", "PKS_KR"])
        assert module.get_monomer("mal") == "redmal"

    def test_completness(self):
        for complete in [[PKS_START, PKS_LOAD, CP],
                         [NRPS_START, NRPS_LOAD, CP],
                         [PKS_LOAD, CP],
                         [NRPS_LOAD, CP],
                         ]:
            assert build_module(complete).is_complete()

        for incomplete in [
                           [PKS_START, CP],
                           [NRPS_START, CP],
                           [NRPS_START, NRPS_LOAD],
                           [NRPS_START, CP],
                           [NRPS_LOAD],
                           ["PKS_KR"],
                           ]:
            assert not build_module(incomplete).is_complete(), incomplete

    def test_is_iterative(self):
        module = build_module(["PKS_KS"], ["Iterative-KS"])
        assert module.is_iterative()

        for other in ["Trans-AT-KS", "Modular-KS"]:
            module = build_module(["PKS_KS"], [other])
            assert not module.is_iterative()

    def test_component_after_end(self):
        add_component(self.nrps, "Epimerization")
        assert self.nrps.is_terminated()
        for other in [NRPS_LOAD, NRPS_START, "nMT"]:
            with self.assertRaisesRegex(ValueError, "adding extra component after end"):
                add_component(self.nrps, other)

    def test_starter_after_others(self):
        with self.assertRaisesRegex(ValueError, "starter after other components"):
            add_component(self.nrps, NRPS_START)

    def test_duplicate_loader(self):
        assert self.nrps._loader
        with self.assertRaisesRegex(ValueError, "duplicate loader"):
            add_component(self.nrps, NRPS_LOAD)

    def test_incompatible_loader(self):
        nrps = build_module([NRPS_START])
        with self.assertRaisesRegex(ValueError, "adding a PKS loader to a NRPS starter"):
            add_component(nrps, PKS_LOAD)

        pks = build_module([PKS_START])
        with self.assertRaisesRegex(ValueError, "adding a NRPS loader to a PKS starter"):
            add_component(pks, NRPS_LOAD)

    def test_duplicate_carrier(self):
        assert self.nrps._carrier_protein
        with self.assertRaisesRegex(ValueError, "duplicate carrier protein"):
            add_component(self.nrps, CP)

    def test_adding_ignored(self):
        start = len(self.nrps._components)
        add_component(self.nrps, "NRPS-COM_Nterm")
        assert len(self.nrps._components) == start

    def test_unknown_component_type(self):
        component = Component(DummyHMMResult(NRPS_LOAD))
        component._domain._hit_id = "unclassifiable"
        component.classification = "unclassifiable"
        with self.assertRaisesRegex(AssertionError, "invalid classification"):
            Module().add_component(component, [])


class TestBuildModules(unittest.TestCase):
    def test_mismatching_ks_subtypes(self):
        with self.assertRaises(StopIteration):
            build_modules_for_cds([DummyHMMResult(PKS_START)], [])

    def test_no_empties(self):
        assert build_modules_for_cds([], []) == []
        assert len(build_modules_for_cds([DummyHMMResult("ACPS")], [])) == 1

    def test_unclassifiable(self):
        with self.assertRaisesRegex(ValueError, "could not classify domain"):
            build_modules_for_cds([DummyHMMResult("UNCLASSIFIABLE")], [])

    def test_module_for_each_starter(self):
        modules = build_modules_for_cds([DummyHMMResult("Condensation_DCL"), DummyHMMResult("Condensation_LCL")], [])
        assert len(modules) == 2

    def test_module_for_orphans(self):
        for name in [NRPS_START, NRPS_LOAD, PKS_START, PKS_LOAD, "cMT", "ACP", "Trans-AT_docking"]:
            modules = build_modules_for_cds([DummyHMMResult(name)], [DummyHMMResult("")])
            assert len(modules) == 1
            assert not modules[0].is_complete()

    def test_bad_add_makes_new_module(self):
        modules = build_modules_for_cds([DummyHMMResult(NRPS_LOAD)] * 2, [])
        assert len(modules) == 2
        assert not modules[0].is_complete()

    def test_starters(self):
        for domain_type in [NRPS_LOAD, PKS_LOAD]:
            domains = [DummyHMMResult(i) for i in [domain_type, "ACP", domain_type, "ACP"]]
            modules = build_modules_for_cds(domains, [])
            print(modules)
            assert len(modules) == 2
            assert modules[0]._first_in_cds
            assert modules[0].is_complete()
            assert not modules[1]._first_in_cds
            print(modules[1], modules[1].is_complete())
            print(modules[1]._starter, modules[1]._loader, modules[1]._carrier_protein,
                  modules[1]._starter is modules[1]._loader)
            assert not modules[1].is_complete()

    def test_double_transporters(self):
        initial_domains = [DummyHMMResult(i) for i in [PKS_START, CP, CP]]
        for special_case in DOUBLE_TRANSPORTER_CASES:
            tail = [DummyHMMResult(i) for i in list(special_case) + ["cMT"]]  # dummy domain to test extent
            # when the exact domains are present, they should extend
            domains = initial_domains + tail
            modules = build_modules_for_cds(domains, ["Trans-AT-PKS"])
            assert len(modules) == 2
            assert len(modules[0].components) == len(domains) - 1
            assert len(modules[1].components) == 1 and modules[1].components[0].domain.hit_id == "cMT"
            # but if there's something else mixed in, they should be treated separately
            domains = initial_domains + [DummyHMMResult("cMT")] + tail
            assert len(build_modules_for_cds(domains, ["Trans-AT-PKS"])) == 3

    def test_double_transporters_miss(self):
        domains = [DummyHMMResult(i) for i in [PKS_START, CP, CP, NRPS_START]]
        modules = build_modules_for_cds(domains, ["Trans-AT-PKS"])
        # the two CP domains should be in separate modules, with the trailing NRPS a third
        assert len(modules) == 3
        assert len(modules[1].components) == 1
        assert modules[1].components[0].domain.hit_id == CP


class TestMerging(unittest.TestCase):
    def setUp(self):
        self.nrps_head = build_module([NRPS_START, NRPS_LOAD])
        self.nrps_tail = build_module([NRPS_LOAD, CP])
        self.nrps_tail._first_in_cds = False
        assert not self.nrps_tail.is_complete()
        self.pks_head = build_module([PKS_START, PKS_LOAD])
        self.pks_tail = build_module([PKS_LOAD, CP])
        self.pks_tail._first_in_cds = False
        assert not self.pks_tail.is_complete()
        self.generic_tail = build_module([CP])
        self.trans_at_head = build_module([PKS_START], subtypes=[TRANS_AT_SUBTYPE])
        self.complete = build_module([PKS_START, PKS_LOAD, CP])

    def build(self, early, late, strand=1, tail_strand_multiplier=1):
        if strand == -1:
            head = late
            tail = early
            second = CDSModuleInfo(DummyCDS(start=50, end=110, strand=strand), [tail])
            first = CDSModuleInfo(DummyCDS(start=500, end=560, strand=strand * tail_strand_multiplier), [head])
        else:
            head = early
            tail = late
            first = CDSModuleInfo(DummyCDS(start=50, end=110, strand=strand), [head])
            second = CDSModuleInfo(DummyCDS(start=500, end=560, strand=strand * tail_strand_multiplier), [tail])

        first_modules = list(first.modules)
        second_modules = list(second.modules)

        if strand == -1:
            module = combine_modules(second, first)
        else:
            module = combine_modules(first, second)

        if not module:
            # nothing should be changed
            assert first_modules == first.modules
            assert second_modules == second.modules
        else:
            # head is replaced
            assert head not in first.modules
            assert len(first_modules) == len(first.modules), (first, second)
            assert module in first.modules
            # tail removed
            assert tail not in second.modules
            # and not replaced
            assert len(second_modules) - 1 == len(second.modules)

        return module

    def test_pks_nrps_combinations(self):
        assert self.nrps_head.is_nrps() and self.pks_tail.is_pks()
        assert not self.nrps_head.is_complete() and not self.pks_tail.is_complete()
        assert not self.build(self.nrps_head, self.pks_tail)
        assert not self.build(self.pks_head, self.nrps_tail)

    def test_invalid_components(self):
        assert not self.build(self.nrps_head, self.nrps_head)
        assert not self.build(self.complete, self.nrps_head)
        assert not self.build(self.nrps_head, self.complete)

    def test_still_incomplete(self):
        # no loaders
        assert not self.build(build_module([NRPS_START]), self.generic_tail)
        assert not self.build(build_module([PKS_START]), self.generic_tail)

    def test_valid_merges(self):
        for head in [self.nrps_head, self.pks_head, self.trans_at_head]:
            module = self.build(head, self.generic_tail)
            assert module.is_complete()

    def test_strand_ordering(self):
        # builds correctly when strands forward
        assert self.build(self.nrps_head, self.generic_tail)
        # but not when strands are reverse
        assert not self.build(self.nrps_head, self.generic_tail, strand=-1)

        # flipping the fragments to match should work in reverse
        assert self.build(self.generic_tail, self.nrps_head, strand=-1)
        # but not in forward
        assert not self.build(self.generic_tail, self.nrps_head)

    def test_strands_must_match(self):
        # works with same strand
        assert self.build(self.nrps_head, self.generic_tail)
        # doesn't work with mismatching strands
        assert not self.build(self.nrps_head, self.generic_tail, tail_strand_multiplier=-1)

    def test_missing_modules(self):
        missing_modules = CDSModuleInfo(DummyCDS(start=50, end=110), [])
        has_modules = CDSModuleInfo(DummyCDS(start=150, end=210), [self.generic_tail])
        assert not combine_modules(missing_modules, has_modules)
        assert not combine_modules(has_modules, missing_modules)

    def test_merge_trailing_kr(self):  # kirromycin/AF484556.1
        head = CDSModuleInfo(DummyCDS(start=50, end=110), [self.trans_at_head])
        tail = CDSModuleInfo(DummyCDS(start=150, end=210), [self.generic_tail, build_module(["PKS_KR"])])
        new = combine_modules(head, tail)
        assert len(tail.modules) == 0
        assert len(head.modules) == 1
        assert head.modules[0] is new
        assert len(new.components) == len(self.trans_at_head.components) + len(self.generic_tail.components) + 1
        assert new.components[-1].domain.hit_id == "PKS_KR"

        # the exceptions to the "rule"
        exceptions = [
            [self.pks_tail, build_module(["PKS_KR"])],  # not a trans-AT
            [self.generic_tail, build_module(["PKS_KR", CP])],  # the KR should stay with its CP
            [self.generic_tail, self.nrps_head],  # completely unrelated next module
        ]
        for modules in exceptions:
            head = CDSModuleInfo(DummyCDS(start=50, end=110), [self.trans_at_head])
            tail = CDSModuleInfo(DummyCDS(start=150, end=210), modules)
            last = [comp.domain.hit_id for comp in modules[-1].components]
            new = combine_modules(head, tail)
            assert len(tail.modules) == 1
            assert len(head.modules) == 1
            # the module after the merged split should be unchanged
            assert [comp.domain.hit_id for comp in tail.modules[0].components] == last
