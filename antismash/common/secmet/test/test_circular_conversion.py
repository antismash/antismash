# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from antismash.common.secmet.locations import FeatureLocation, CompoundLocation
from antismash.common.secmet.record import Record


class TestBridgeConversion(unittest.TestCase):
    def setUp(self):
        self.seqrec = SeqRecord(UnknownSeq(21))
        loc = CompoundLocation([FeatureLocation(12, 15, strand=1),
                                FeatureLocation(18, 21, strand=1),
                                FeatureLocation(0, 3, strand=1),
                                FeatureLocation(6, 9, strand=1)],
                                operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"

    def test_bridge_in_linear_record(self):
        self.seqrec.annotations["topology"] = "linear"
        self.seqrec.features.append(self.seqcds)
        with self.assertRaisesRegex(ValueError, "Features that bridge"):
            Record.from_biopython(self.seqrec, taxon='bacteria')
        self.seqrec.features[0] = self.seqgene
        with self.assertRaisesRegex(ValueError, "Features that bridge"):
            Record.from_biopython(self.seqrec, taxon='bacteria')

    def test_cds_split(self):
        self.seqrec.features.append(self.seqcds)
        print(self.seqcds)
        for id_name in ["locus_tag", "gene"]:
            self.seqcds.qualifiers[id_name] = ["test"]
            rec = Record.from_biopython(self.seqrec, taxon="bacteria")
            cdses = rec.get_cds_features()
            assert len(cdses) == 2

            assert cdses[0].location.start == 0
            assert cdses[0].location.end == 9
            assert getattr(cdses[0], id_name) == "test_LOWER"
            assert cdses[0].get_name() == "test_LOWER"

            assert cdses[1].location.start == 12
            assert cdses[1].location.end == 21
            assert getattr(cdses[1], id_name) == "test_UPPER"
            assert cdses[1].get_name() == "test_UPPER"

            self.seqcds.qualifiers.pop(id_name)

    def test_gene_split(self):
        self.seqrec.features.append(self.seqgene)
        for id_name in ["locus_tag", "gene"]:
            self.seqgene.qualifiers[id_name] = [id_name + "_test"]
            expected = id_name + "_test"
            rec = Record.from_biopython(self.seqrec, taxon="bacteria")
            self.seqgene.qualifiers.pop(id_name)
            genes = rec.get_genes()
            assert len(genes) == 2

            if id_name == "gene":
                id_name = "gene_name"  # since a Gene doesn't have a gene member

            assert genes[0].location.start == 12
            assert genes[0].location.end == 21
            assert getattr(genes[0], id_name) == expected + "_UPPER"
            assert genes[0].get_name() == expected + "_UPPER"

            assert genes[1].location.start == 0
            assert genes[1].location.end == 9
            assert getattr(genes[1], id_name) == expected + "_LOWER"
            assert genes[1].get_name() == expected + "_LOWER"


    def test_cds_with_no_id(self):
        self.seqrec.features.append(self.seqcds)
        rec = Record.from_biopython(self.seqrec, taxon="bacteria")
        cdses = rec.get_cds_features()
        assert len(cdses) == 2
        assert cdses[0].location.start == 0
        assert cdses[0].location.end == 9
        assert cdses[0].get_name() == "bridge_LOWER"

        assert cdses[1].location.start == 12
        assert cdses[1].location.end == 21
        assert cdses[1].get_name() == "bridge_UPPER"

    def test_gene_with_no_id(self):
        self.seqrec.features.append(self.seqgene)
        rec = Record.from_biopython(self.seqrec, taxon="bacteria")
        genes = rec.get_genes()
        assert len(genes) == 2

        assert genes[0].location.start == 12
        assert genes[0].location.end == 21
        assert genes[0].get_name() == "bridge_UPPER"

        assert genes[1].location.start == 0
        assert genes[1].location.end == 9
        assert genes[1].get_name() == "bridge_LOWER"


class TestSingleLower(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(UnknownSeq(21))
        loc = CompoundLocation([FeatureLocation(12, 15, strand=1),
                                FeatureLocation(18, 21, strand=1),
                                FeatureLocation(0, 9, strand=1)],
                                operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"


class TestSingleUpper(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(UnknownSeq(21))
        loc = CompoundLocation([FeatureLocation(12, 21, strand=1),
                                FeatureLocation(0, 3, strand=1),
                                FeatureLocation(6, 9, strand=1)],
                                operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"

class TestSingleBoth(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(UnknownSeq(21))
        loc = CompoundLocation([FeatureLocation(12, 21, strand=1),
                                FeatureLocation(0, 9, strand=1)],
                                operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"
