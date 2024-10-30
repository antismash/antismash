# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.features import CDSCollection
from antismash.common.secmet.locations import CompoundLocation, FeatureLocation
from antismash.common.secmet.test.helpers import (
    DummyCandidateCluster,
    DummyProtocluster,
    DummyRegion,
    DummySubRegion,
)
from antismash.outputs.html import area_packing
from antismash.outputs.html.area_packing import (
    Area,
    Row,
)


class TestArea(unittest.TestCase):
    def test_default_neighbouring(self):
        area = Area(start=5, end=8, kind="dummykind")
        assert area.neighbouring_start == area.start == 5
        assert area.neighbouring_end == area.end == 8

        area = Area(start=5, end=8, kind="dummykind", neighbouring_start=3, neighbouring_end=9)
        assert area.neighbouring_start == 3
        assert area.start == 5
        assert area.end == 8
        assert area.neighbouring_end == 9

    def test_minimal(self):
        area = Area(start=5, end=8, kind="dummykind")
        # neighbouring coordinates that are the same shouldn't be included
        converted = area.to_minimal_json()
        # including height, since it will always be present for drawing
        assert set(converted) == {"start", "end", "kind", "height"}

        # but when they differ, they should be included
        area.neighbouring_start = 2
        converted = area.to_minimal_json()
        assert set(converted) == {"start", "end", "kind", "height", "neighbouring_start"}

    def test_cloning(self):
        original = Area(start=5, end=8, kind="dummykind", product="test product")
        # no grouping should exist
        assert not original.group

        clone = original.clone()
        # and now grouping should exist on both and match
        assert original.group and original.group == clone.group

        # and if the clone is then cloned once more, the groups should all match
        cloned_clone = clone.clone()
        assert original.group
        assert len({original.group, clone.group, cloned_clone.group}) == 1

    def test_protocluster_conversion(self):
        proto = DummyProtocluster()
        area = Area.from_feature(proto)
        assert area.start == proto.core_start
        assert area.end == proto.core_end
        assert area.neighbouring_start == proto.start
        assert area.neighbouring_end == proto.end
        assert area.kind == proto.FEATURE_TYPE
        assert area.product == proto.product
        assert area.category == proto.product_category
        assert area.tool == proto.tool

    def test_candidate_cluster_conversion(self):
        candidate = DummyCandidateCluster()
        area = Area.from_feature(candidate)
        assert area.start == candidate.start
        assert area.end == candidate.end
        assert area.kind == "candidatecluster"  # matches the expected types of the javascript
        assert str(candidate.kind) in area.product

    def test_subregion_conversion(self):
        sub = DummySubRegion()
        area = Area.from_feature(sub)
        assert area.start == sub.start
        assert area.end == sub.end
        assert area.kind == sub.FEATURE_TYPE
        assert area.product == sub.label

    def test_generic_collection(self):
        feature = CDSCollection(FeatureLocation(1, 5, 1), feature_type="some_type")
        area = Area.from_feature(feature)
        assert area.start == feature.start
        assert area.end == feature.end

    def test_extra_args(self):
        # not provided by the feature
        proto = DummyProtocluster()
        area = Area.from_feature(proto, height=7)
        assert area.height == 7
        # overriding the feature
        area = Area.from_feature(DummyProtocluster(), product="override")
        assert area.product == "override"
        assert area.product != proto.product

    def test_offset(self):
        start = 5
        end = 7
        neighbouring_start = 3
        neighbouring_end = 8
        for offset in [2, -1]:
            area = Area(start=start, end=end, neighbouring_start=neighbouring_start,
                        neighbouring_end=neighbouring_end, kind="dummykind")
            # no cheating, the initial start must match the provided start
            assert area.start == start
            area.offset(offset)
            assert area.start == start + offset
            assert area.end == end + offset
            assert area.neighbouring_start == neighbouring_start + offset
            assert area.neighbouring_end == neighbouring_end + offset


class TestAreaAdjustment(unittest.TestCase):
    def check(self, feature, length, region_crosses_origin):
        assert feature.crosses_origin()
        area = Area.from_feature(feature)
        assert area.crosses_origin()
        extra = area_packing.adjust_cross_origin_area(area, feature, region_crosses_origin, length)
        # regardless of the region crossing the origin, the original area must no longer cross the origin
        assert not area.crosses_origin()
        # but it does affect whether or not the area is split
        if region_crosses_origin:
            assert extra is None
        else:
            assert extra
            # if a split occurs, the new portion must not cross the origin
            assert not extra.crosses_origin()
        return area, extra

    def test_invalid(self):
        area = Area(start=1, end=5, kind="test")
        for region_crosses_origin in [False, True]:
            with self.assertRaisesRegex(ValueError, "feature must cross the origin"):
                area_packing.adjust_cross_origin_area(area, DummySubRegion(), region_crosses_origin, 100)

    def test_protocluster_left_neighbourhood(self):
        length = 100
        proto = DummyProtocluster(start=80, end=50, core_start=10, core_end=20, record_length=length)
        # when the region crosses the origin, it's only coordinate extensions
        area, extra = self.check(proto, length, region_crosses_origin=True)
        assert extra is None
        assert area.neighbouring_start == proto.start
        assert area.start == proto.core_start + length
        assert area.end == proto.core_end + length
        assert area.neighbouring_end == proto.end + length

        # but when the region doesn't cross the origin, then the area must split
        area, extra = self.check(proto, length, region_crosses_origin=False)
        assert extra

        assert area.neighbouring_start == proto.start
        assert area.start == length
        assert area.end == length
        assert area.neighbouring_end == length

        assert extra.neighbouring_start == 0
        assert extra.start == proto.core_start
        assert extra.end == proto.core_end
        assert extra.neighbouring_end == proto.end

    def test_protocluster_right_neighbourhood(self):
        length = 100
        proto = DummyProtocluster(start=50, end=20, core_start=70, core_end=80, record_length=length)
        # cross-origin regions mean extensions of the area coordinates
        area, extra = self.check(proto, length, region_crosses_origin=True)
        assert extra is None
        assert area.neighbouring_start == proto.start
        assert area.start == proto.core_start
        assert area.end == proto.core_end
        assert area.neighbouring_end == proto.end + length

        # but when the region doesn't cross the origin, then the area must split
        area, extra = self.check(proto, length, region_crosses_origin=False)
        assert extra

        assert area.neighbouring_start == proto.start
        assert area.start == proto.core_start
        assert area.end == proto.core_end
        assert area.neighbouring_end == length

        assert extra.neighbouring_start == 0
        assert extra.start == 0
        assert extra.end == 0
        assert extra.neighbouring_end == proto.end

    def test_protocluster_core(self):
        length = 100
        proto = DummyProtocluster(start=80, end=20, core_start=90, core_end=10, record_length=length)

        area, extra = self.check(proto, length, region_crosses_origin=True)
        assert extra is None

        assert area.neighbouring_start == proto.start
        assert area.start == proto.core_start
        assert area.end == proto.core_end + length
        assert area.neighbouring_end == proto.end + length

        area, extra = self.check(proto, length, region_crosses_origin=False)
        assert extra

        assert area.neighbouring_start == proto.start
        assert area.start == proto.core_start
        assert area.end == length
        assert area.neighbouring_end == length

        assert extra.neighbouring_start == 0
        assert extra.start == 0
        assert extra.end == proto.core_end
        assert extra.neighbouring_end == proto.end

    def test_subregion(self):
        length = 100
        subregion = DummySubRegion(start=50, end=30, record_length=length)
        area, extra = self.check(subregion, length, region_crosses_origin=True)
        assert extra is None
        assert area.start == area.neighbouring_start == subregion.start
        assert area.end == area.neighbouring_end == subregion.end + length

        area, extra = self.check(subregion, length, region_crosses_origin=False)
        assert extra
        assert area.start == area.neighbouring_start == subregion.start
        assert area.end == area.neighbouring_end == length

        assert extra.start == extra.neighbouring_start == 0
        assert extra.end == extra.neighbouring_end == subregion.end

    def test_candidate_cluster(self):
        length = 100
        candidate = DummyCandidateCluster(start=50, end=30, circular_wrap_point=length)
        area, extra = self.check(candidate, length, region_crosses_origin=True)
        assert extra is None
        # ignoring cores, since these are dependent on the protoclusters anyway
        assert area.neighbouring_start == candidate.start
        assert area.neighbouring_end == candidate.end + length

        area, extra = self.check(candidate, length, region_crosses_origin=False)
        assert extra
        assert area.neighbouring_start == candidate.start
        assert area.neighbouring_end == length

        assert extra.neighbouring_start == 0
        assert extra.neighbouring_end == candidate.end


class TestRow(unittest.TestCase):
    def test_init(self):
        row = Row()
        assert row.start == 0
        assert not row.contents

        row = Row(start=20, initial_contents=[DummySubRegion(start=21, end=25)])
        assert row.start == 26
        assert len(row.contents) == 1

    def test_init_contents(self):
        # all provided contents must fit, in the order they're given in
        with self.assertRaisesRegex(ValueError, "cannot fit area"):
            Row(initial_contents=[DummySubRegion(start=1, end=5),
                                  DummySubRegion(start=3, end=7)])

        # the provided start is enforced for the contents
        with self.assertRaisesRegex(ValueError, "cannot fit area"):
            Row(start=10, initial_contents=[DummySubRegion(start=5, end=15)])

        areas = [DummySubRegion(start=1, end=5), DummySubRegion(start=10, end=15)]
        row = Row(initial_contents=areas)
        assert row.start == 16
        assert row.contents == tuple(areas)

    def test_adding(self):
        row = Row()
        assert row.start == 0
        assert not row.contents

        row.add(DummySubRegion(start=5, end=20))
        assert len(row.contents) == 1
        assert row.start == 21

        row.add(DummySubRegion(start=30, end=40))
        assert len(row.contents) == 2
        assert row.start == 41

    def test_adding_invalid(self):
        row = Row()
        row.add(DummySubRegion(start=5, end=10))
        with self.assertRaisesRegex(ValueError, "cannot fit area"):
            row.add(DummySubRegion(start=8, end=20))

    def test_can_fit(self):
        row = Row()
        areas = [
            DummySubRegion(start=5, end=10),
            DummySubRegion(start=8, end=15),
            DummySubRegion(start=20, end=30),
        ]
        for area in areas:
            assert row.can_fit(area)

        row.add(areas[0])
        assert row.start == 11
        assert areas[1].start < row.start
        assert areas[2].start > row.start
        assert not row.can_fit(areas[1])
        assert row.can_fit(areas[2])

    def test_cross_origin_insert_before(self):
        # a cross-origin area can pack into the area before the first existing area
        row = Row()
        row.add(DummySubRegion(start=10, end=20))
        assert row.start == 21
        cross_origin = DummySubRegion(location=CompoundLocation([
            FeatureLocation(25, 30, 1),
            FeatureLocation(0, 5, 1),
        ]))
        row.add(cross_origin)
        assert row.start == 30
        assert len(row.contents) == 2
        # and now that the end is bounded, attempting to add anything after it should fail
        with self.assertRaisesRegex(ValueError, "cannot fit area"):
            row.add(DummySubRegion(start=35, end=45))


class TestHeights(unittest.TestCase):
    def test_only_protoclusters(self):
        # since protoclusters have labels above the bars, the heights need to allow for that
        # and this will happen in regions where there's a single protoclusters and thus a
        # single, matching candidate cluster
        protocluster = DummyProtocluster()
        candidate = DummyCandidateCluster(clusters=[protocluster])
        region = DummyRegion(candidate_clusters=[candidate])

        results = area_packing.build_area_rows(region, 100)
        assert len(results) == 1
        assert results[0]["height"] == 1

    def test_only_subregions(self):
        region = DummyRegion(subregions=[DummySubRegion()])
        results = area_packing.build_area_rows(region, 100)
        assert len(results) == 1
        assert results[0]["height"] == 0

    def test_complicated(self):
        # a circular record, with multiple protoclusters, some crossing the origin,
        # some not, with a parent candidate cluster that also crosses the origin
        # and a single subregion that spans the gap between while also crossing
        # the origin

        record_length = 98206
        protoclusters = [
            # core crosses the origin
            DummyProtocluster(start=65398, end=42440, core_start=85398, core_end=22440,
                              product="core_over_origin", record_length=record_length),
            # neighbouring crosses the origin
            DummyProtocluster(start=7053, end=52159, core_start=27053, core_end=32159,
                              product="mid_contig", record_length=record_length),
            # simple, non-crossing cluster
            DummyProtocluster(start=93497, end=52159, core_start=15291, core_end=32159,
                              product="neighbourhood_over_origin", record_length=record_length),
        ]

        # and the bridging subregion, covering the previously uncovered area of the record
        subregion = DummySubRegion(start=90000, end=70000, record_length=record_length)

        # finally, the region doesn't cross the origin, but spans the whole record
        region = DummyRegion(candidate_clusters=[DummyCandidateCluster(clusters=protoclusters)],
                             subregions=[subregion])

        areas = area_packing.build_area_rows(region, record_length, circular=True)
        # all cross-origin features will have been split into two
        # so no area should cross the origin
        for area in areas:
            start = area.get("neighbouring_start", area["start"])
            end = area.get("neighbouring_end", area["end"])
            assert start <= end, area
            assert area["start"] <= area["end"], area
        # and those features should count double
        proto_areas = [area for area in areas if area["kind"] == "protocluster"]
        assert len(proto_areas) == 5  # two cross-origin, one not

        candidate_areas = [area for area in areas if area["kind"] == "candidatecluster"]
        assert len(candidate_areas) == 2

        subregion_areas = [area for area in areas if area["kind"] == "subregion"]
        assert len(subregion_areas) == 2
        assert len(areas) == 9


        step_size = 2
        # candidates are on top, then subregions, then protoclusters
        post_origin, pre_origin = sorted(candidate_areas, key=lambda x: x["start"])
        assert post_origin["height"] == pre_origin["height"] == 0

        post_origin, pre_origin = sorted(subregion_areas, key=lambda x: x["start"])
        assert post_origin["height"] == step_size  # only one candidate above

        for area in proto_areas:
            assert area["height"] > step_size  # below the candidate
            assert area["height"] % step_size == 0  # each should be on the step, not part way
