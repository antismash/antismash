# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from antismash.common.secmet.locations import FeatureLocation, CompoundLocation
from antismash.common.secmet.record import Record, SecmetInvalidInputError


class TestBridgeConversion(unittest.TestCase):
    def setUp(self):
        self.seqrec = SeqRecord(Seq("A"*21))
        loc = CompoundLocation([FeatureLocation(12, 15, strand=1),
                                FeatureLocation(18, 21, strand=1),
                                FeatureLocation(0, 3, strand=1),
                                FeatureLocation(6, 9, strand=1)],
                               operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"
        self.seqrec.annotations["molecule_type"] = "DNA"

    def test_bridge_in_linear_record(self):
        self.seqrec.annotations["topology"] = "linear"
        self.seqrec.features.append(self.seqcds)
        with self.assertRaisesRegex(SecmetInvalidInputError, "cannot determine correct exon ordering"):
            Record.from_biopython(self.seqrec, taxon='fungi')
        self.seqrec.features[0] = self.seqgene
        with self.assertRaisesRegex(SecmetInvalidInputError, "cannot determine correct exon ordering"):
            Record.from_biopython(self.seqrec, taxon='fungi')


class TestSingleLower(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(Seq("A"*21))
        loc = CompoundLocation([FeatureLocation(12, 15, strand=1),
                                FeatureLocation(18, 21, strand=1),
                                FeatureLocation(0, 9, strand=1)],
                               operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"


class TestSingleUpper(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(Seq("A"*21))
        loc = CompoundLocation([FeatureLocation(12, 21, strand=1),
                                FeatureLocation(0, 3, strand=1),
                                FeatureLocation(6, 9, strand=1)],
                               operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"


class TestSingleBoth(TestBridgeConversion):
    def setUp(self):
        self.seqrec = SeqRecord(Seq("A"*21))
        loc = CompoundLocation([FeatureLocation(12, 21, strand=1),
                                FeatureLocation(0, 9, strand=1)],
                               operator="join")
        self.seqcds = SeqFeature(loc, type="CDS")
        self.seqgene = SeqFeature(loc, type="gene")
        self.seqrec.annotations["topology"] = "circular"
