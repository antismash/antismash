# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet import Module
from antismash.common.test.helpers import DummyCDS, DummyHMMResult, DummyRecord
from antismash.detection.nrps_pks_domains import domain_drawing
from antismash.detection.nrps_pks_domains.domain_identification import (
    CDSModuleInfo,
    CDSResult,
    build_modules_for_cds,
    combine_modules,
)


class TestModuleJSON(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord()
        self.head = None
        self.tail = None
        self.head_hits = []
        self.tail_hits = []

    def add_module(self):
        self.record.add_cds_feature(self.tail)
        self.record.add_cds_feature(self.head)

        head_modules = build_modules_for_cds(self.head_hits, self.head.get_name())
        tail_modules = build_modules_for_cds(self.tail_hits, self.tail.get_name())
        assert str(tail_modules[0]) == "[CP]"
        assert str(head_modules[0]) == "[C,A]"
        tail_info = CDSModuleInfo(self.tail, tail_modules)
        head_info = CDSModuleInfo(self.head, head_modules)
        if self.head.location.strand == -1:
            combine_modules(tail_info, head_info)
        else:
            combine_modules(head_info, tail_info)
        assert not tail_modules
        assert str(head_modules[0]) == "[C,A,CP]"

        CDSResult(self.head_hits, [], head_modules).annotate_domains(self.record, self.head)
        CDSResult(self.tail_hits, [], []).annotate_domains(self.record, self.tail)

        domains = []
        for cds in [self.tail, self.head][::self.head.location.strand]:  # keep strand in mind
            domains.extend(self.record.get_antismash_domains_in_cds(cds))

        module = Module(domains, Module.types.NRPS, complete=True)
        self.record.add_module(module)

        # ensure it's constructed as expected
        assert module.is_complete()
        # head should always be first, reverse strand or not
        assert module.parent_cds_names == ("head", "tail")
        assert len(module.domains) == 3
        assert module.domains[0].locus_tag != module.domains[-1].locus_tag
        assert module.domains[0].locus_tag == "head"

        # and is linked as expected
        assert module in self.head.modules
        assert module in self.tail.modules

        return module

    def check_conversion(self, module):
        match_gen = iter(["A", "B"])
        head_json = domain_drawing._build_module_js(module, self.head, {}, match_gen)

        assert head_json.start == 31
        assert head_json.end == 674
        # these, and the matching tail section below, are the vital parts
        assert head_json.match_id == "A"
        assert head_json.multi_cds == "head"

        tail_json = domain_drawing._build_module_js(module, self.tail, {}, match_gen)
        assert tail_json.start == 45
        assert tail_json.end == 114
        assert tail_json.match_id == "B"
        assert tail_json.multi_cds == "tail"

    def test_full_path_reverse_strand(self):
        self.tail = DummyCDS(locus_tag="tail", start=19_580, end=19_997, strand=-1)
        self.tail_hits = [DummyHMMResult("PCP", start=45, end=114)]

        self.head = DummyCDS(locus_tag="head", start=20_000, end=22_547, strand=-1)
        self.head_hits = [
            DummyHMMResult("Condensation_DCL", start=31, end=170),
            DummyHMMResult("AMP-binding", start=354, end=674),
        ]

        module = self.add_module()

        assert module.domains[0].location.strand == -1

        self.check_conversion(module)

    def test_full_path_forward_strand(self):
        self.head = DummyCDS(locus_tag="head", start=20_000, end=22_547, strand=1)
        self.head_hits = [
            DummyHMMResult("Condensation_DCL", start=31, end=170),
            DummyHMMResult("AMP-binding", start=354, end=674),
        ]
        self.tail = DummyCDS(locus_tag="tail", start=22_580, end=22_997, strand=1)
        self.tail_hits = [DummyHMMResult("PCP", start=45, end=114)]

        module = self.add_module()

        assert module.domains[0].location.strand == 1
        self.check_conversion(module)
