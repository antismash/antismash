# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,assigning-non-slot

import unittest

from antismash.common.secmet.locations import FeatureLocation
from antismash.common.secmet.features.module import Module, ModuleType
from antismash.common.secmet.test.helpers import DummyCDS, DummyAntismashDomain
from antismash.common.test.helpers import DummyRecord


def create_module(domains=None, **kwargs):
    if domains is None:
        domains = [DummyAntismashDomain()]
    return Module(domains, **kwargs)


def add_module_references_to_record(module, record):
    for i, domain in enumerate(module.domains):
        record.add_antismash_domain(domain)
        try:
            record.get_cds_by_name(domain.locus_tag)
        except KeyError:
            record.add_cds_feature(DummyCDS(start=max(i, module.location.start - 10),
                                            end=module.location.end + 10 + i,
                                            locus_tag=domain.locus_tag))


class TestModuleType(unittest.TestCase):
    def test_conversion(self):
        for mtype in ModuleType:
            assert ModuleType.from_string(str(mtype)) is mtype

        with self.assertRaisesRegex(ValueError, "unknown module type"):
            ModuleType.from_string("other")


class TestModule(unittest.TestCase):
    def test_minimal_biopython_conversion(self):
        original = create_module()

        bio = original.to_biopython()
        assert isinstance(bio, list) and len(bio) == 1
        assert bio[0].location == original.location
        expected_qualifiers = {"domains", "locus_tags", "incomplete"}
        assert set(bio[0].qualifiers.keys()).issuperset(expected_qualifiers)
        assert bio[0].qualifiers["incomplete"] is None
        assert "complete" not in bio[0].qualifiers

        record = DummyRecord()
        add_module_references_to_record(original, record)

        final = Module.from_biopython(bio[0], record=record)

        assert final.location == original.location
        assert not final.is_starter_module()
        assert not final.is_final_module()
        assert not final.is_iterative()
        assert not final.monomers
        assert not final.is_complete()
        assert final.location == original.location

    def test_detailed_biopython_conversion(self):
        expected_qualifiers = {"domains", "locus_tags", "complete", "starter_module",
                               "final_module", "monomer_pairings"}
        original = create_module(starter=True, final=True, complete=True, iterative=True)
        original.add_monomer("from", "to")
        bio = original.to_biopython()
        assert set(bio[0].qualifiers.keys()).issuperset(expected_qualifiers)
        assert bio[0].qualifiers["complete"] is None
        assert "incomplete" not in bio[0].qualifiers
        assert all(domain and " " not in domain for domain in bio[0].qualifiers["domains"])
        # fake an inserted space due to length
        first = bio[0].qualifiers["domains"][0]
        half = len(first)//2
        bio[0].qualifiers["domains"][0] = first[:half] + " " + first[half:]

        record = DummyRecord()
        add_module_references_to_record(original, record)

        final = Module.from_biopython(bio[0], record=record)

        assert final.location == original.location
        assert final.is_starter_module()
        assert final.is_final_module()
        assert final.is_iterative()
        assert final.monomers == original.monomers
        assert final.is_complete()
        assert final.location == original.location

    def test_bad_biopython_conversion(self):
        record = DummyRecord()
        module = create_module()
        add_module_references_to_record(module, record)

        for removed_key in ["domains", "incomplete"]:
            bio = module.to_biopython()[0]
            bio.qualifiers.pop(removed_key)
            with self.assertRaisesRegex(ValueError, f"missing .* '?{removed_key}"):
                Module.from_biopython(bio, record=record)

        bio = module.to_biopython()[0]
        bio.qualifiers["type"] = ["not a valid type"]
        with self.assertRaisesRegex(ValueError, "unknown module type"):
            Module.from_biopython(bio, record=record)

        bio = module.to_biopython()[0]
        for bad in ["missing", "has -> some -> extra", "missing ->", "-> missing", "->"]:
            bio.qualifiers["monomer_pairings"] = [bad]
            with self.assertRaisesRegex(ValueError, "invalid monomer pairing"):
                Module.from_biopython(bio, record=record)

    def test_conversion_bad_record(self):
        bio = create_module().to_biopython()[0]

        with self.assertRaisesRegex(ValueError, "record instance required"):
            Module.from_biopython(bio)

        with self.assertRaisesRegex(ValueError, "does not contain domain referenced by module"):
            Module.from_biopython(bio, record=DummyRecord())

    def test_parents(self):
        cds = DummyCDS(0, 6, locus_tag="testCDS")
        domain = DummyAntismashDomain(2, 5)
        domain.locus_tag = "testCDS"
        module = create_module([domain])
        assert module.parent_cds_names[0] == cds.get_name()

    def test_starter(self):
        assert not create_module().is_starter_module()
        assert create_module(starter=True).is_starter_module()
        assert not create_module(starter=False).is_starter_module()

    def test_final(self):
        assert not create_module().is_final_module()
        assert create_module(final=True).is_final_module()
        assert not create_module(final=False).is_final_module()

    def test_complete(self):
        assert not create_module().is_complete()
        assert create_module(complete=True).is_complete()
        assert not create_module(complete=False).is_complete()

    def test_iterative(self):
        assert not create_module().is_iterative()
        assert create_module(iterative=True).is_iterative()
        assert not create_module(iterative=False).is_iterative()

    def test_cross_gene(self):
        module = create_module()
        assert len(module.parent_cds_names) == 1 and not module.is_multigene_module()
        module._parent_cds_names.append("more")
        assert len(module.parent_cds_names) == 2 and module.is_multigene_module()

    def test_cross_gene_ordering(self):
        domains = [DummyAntismashDomain(locus_tag="A"), DummyAntismashDomain(locus_tag="B")]
        for doms in [domains, domains[::-1]]:
            assert Module(doms).parent_cds_names == tuple([dom.locus_tag for dom in doms])

    def test_type(self):
        assert create_module().module_type == ModuleType.UNKNOWN

        for kind in ModuleType:
            assert create_module(module_type=kind).module_type == kind

        for bad in [None, "unknown"]:
            with self.assertRaises(TypeError):
                create_module(module_type=bad)

    def test_no_domains(self):
        with self.assertRaisesRegex(ValueError, "at least one domain required"):
            create_module(domains=[])

    def test_no_random_attributes(self):
        with self.assertRaises(AttributeError):
            create_module().something = "value"

    def test_monomers_and_substrates(self):
        module = create_module()
        module.add_monomer("from", "to")
        module.add_monomer("other", "thing")
        assert module.monomers == (("from", "to"), ("other", "thing"))

        with self.assertRaisesRegex(ValueError, "substrate is required"):
            module.add_monomer("", "to")
        with self.assertRaisesRegex(ValueError, "monomer is required"):
            module.add_monomer("from", "")

    def test_multi_cds_sorting(self):
        # regardless of strand, the domains should be ordered as per a module
        forward = [DummyAntismashDomain(start=20, end=50, protein_start=7, protein_end=17, locus_tag="A"),
                   DummyAntismashDomain(start=70, end=91, protein_start=3, protein_end=10, locus_tag="B")]
        reverse = [DummyAntismashDomain(start=70, end=91, strand=-1, protein_start=3, protein_end=10, locus_tag="A"),
                   DummyAntismashDomain(start=20, end=50, strand=-1, protein_start=7, protein_end=17, locus_tag="B")]

        for domains in [forward, reverse]:
            module = create_module(domains=domains)
            assert module.is_multigene_module()
            alternate = create_module(domains=domains[::-1])
            assert module.domains == alternate.domains
            assert list(module.domains) == domains

    def test_multi_cds_mismatching_strands(self):
        domains = [DummyAntismashDomain(start=20, end=50, strand=1, protein_start=3, protein_end=10, locus_tag="A"),
                   DummyAntismashDomain(start=70, end=91, strand=-1, protein_start=7, protein_end=17, locus_tag="B")]
        with self.assertRaisesRegex(ValueError, "cannot be on different strands"):
            create_module(domains)

    def test_multi_cds_tracking(self):
        domains = [DummyAntismashDomain(locus_tag=i) for i in "AB"]
        module = create_module(domains=domains)
        assert module.is_multigene_module()
        record = DummyRecord()
        add_module_references_to_record(module, record)
        record.add_cds_feature(DummyCDS(locus_tag="C"))
        for cds in record.get_cds_features():
            assert not cds.modules
        assert not record.get_modules()
        record.add_module(module)
        # make sure it's not added to every CDS
        assert not record.get_cds_by_name("C").modules
        # but that it is added to all CDSes with a domain included
        for i in "AB":
            assert record.get_cds_by_name(i).modules == (module,)

    def test_multi_cds_protein_location(self):
        domains = [DummyAntismashDomain(locus_tag=i, protein_start=n, protein_end=n+5) for n, i in enumerate("AB")]
        module = create_module(domains=domains)
        assert module.is_multigene_module()
        with self.assertRaisesRegex(ValueError, "cannot generate protein location for multi"):
            _ = module.protein_location

        assert module.get_parent_protein_location("A") == FeatureLocation(0, 5)
        assert module.get_parent_protein_location("B") == FeatureLocation(1, 6)
        with self.assertRaisesRegex(ValueError, "has no parent"):
            module.get_parent_protein_location("C")
