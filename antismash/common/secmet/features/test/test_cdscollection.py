# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import CDSCollection, FeatureLocation
from antismash.common.secmet.locations import CompoundLocation
from antismash.common.secmet.test.helpers import DummyCDS, DummyRecord


class TestCDSCollection(unittest.TestCase):
    def test_cds_caching(self):
        collection = CDSCollection(FeatureLocation(0, 200, 1), "dummy")
        assert not collection.cds_children
        collection.add_cds(DummyCDS(10, 40, strand=1))
        original = collection.cds_children
        assert len(original) == 1
        assert collection.cds_children is original  # same object, since no change
        collection.add_cds(DummyCDS(110, 140, strand=-1))
        updated = collection.cds_children
        assert len(updated) == 2
        assert updated is not original

    def test_parent_linkage(self):
        child = CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[])
        assert child.parent is None
        parent = CDSCollection(FeatureLocation(10, 50), feature_type="test", child_collections=[child])
        assert child.parent is parent

    def test_bad_child(self):
        with self.assertRaises(AssertionError):
            child = CDSCollection(FeatureLocation(10, 50), feature_type="test", child_collections=[])
            CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[child])

        with self.assertRaises(AssertionError):
            cds = DummyCDS(25, 35)
            CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[cds])

    def test_root(self):
        child = CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[])
        assert child.get_root() is child
        parent = CDSCollection(FeatureLocation(10, 50), feature_type="test", child_collections=[child])
        assert child.get_root() is parent
        grandparent = CDSCollection(FeatureLocation(0, 60), feature_type="test", child_collections=[parent])
        for col in [child, parent, grandparent]:
            assert col.get_root() is grandparent

    def test_add_cds(self):
        collection = CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[])
        cds = DummyCDS(20, 40)
        collection.add_cds(cds)
        assert cds in collection.cds_children

        cds = DummyCDS(120, 140)
        with self.assertRaisesRegex(ValueError, "not contained by"):
            collection.add_cds(cds)

    def test_contig_edge_transitivity(self):
        inner = CDSCollection(FeatureLocation(30, 40), feature_type="test")
        mid = CDSCollection(FeatureLocation(20, 50), feature_type="test", child_collections=[inner])
        outer = CDSCollection(FeatureLocation(10, 60), feature_type="test", child_collections=[mid])
        collections = [inner, mid, outer]

        # without a record, this breaks
        with self.assertRaisesRegex(ValueError, "Cannot determine"):
            outer.contig_edge  # pylint: disable=pointless-statement

        # add the record, ensure for testing that everything starts not on the edge
        record = DummyRecord(seq="A"*70)
        for collection in collections:
            collection.parent_record = record
            assert not collection.contig_edge

        # set the inner to be on the edge, ensure that all parents now respect that
        inner._contig_edge = True
        for collection in collections:
            assert collection.contig_edge

        # reset
        for collection in collections:
            collection._contig_edge = False

        # set the mid only to be on the edge, to be sure it doesn't apply to children
        mid._contig_edge = True
        assert not inner.contig_edge
        assert mid.contig_edge and outer.contig_edge

    def test_ordering(self):
        shorter = CDSCollection(FeatureLocation(20, 40), feature_type="test", child_collections=[])
        longer = CDSCollection(FeatureLocation(20, 50), feature_type="test", child_collections=[])
        assert shorter > longer

        cross_origin = CDSCollection(CompoundLocation([FeatureLocation(60, 80, 1), FeatureLocation(0, 40, 1)]),
                                     feature_type="test", child_collections=[])
        assert sorted([shorter, longer, cross_origin]) == [cross_origin, longer, shorter]

    def test_without_cross_origin(self):
        # simple locations are fine and obviously can't cross the origin
        area = CDSCollection(FeatureLocation(50, 70), feature_type="dummy_type")
        assert not area.crosses_origin()

        # compound locations that don't cross the origin aren't allowed
        with self.assertRaisesRegex(ValueError, "without crossing the origin"):
            area = CDSCollection(CompoundLocation([
                FeatureLocation(0, 20),
                FeatureLocation(50, 70),
            ]), feature_type="dummy_type")

    def test_cross_origin(self):
        area = CDSCollection(CompoundLocation([
            FeatureLocation(150, 200, 1),
            FeatureLocation(0, 20, 1),
        ]), feature_type="dummy_type")
        assert area.crosses_origin()

        inner = CDSCollection(CompoundLocation([
            FeatureLocation(150, 200, 1),
            FeatureLocation(0, 20, 1),
        ]), feature_type="dummy_type")
        outer = CDSCollection(FeatureLocation(0, 200), feature_type="test", child_collections=[inner])
        assert not outer.crosses_origin()
        assert inner.crosses_origin()
