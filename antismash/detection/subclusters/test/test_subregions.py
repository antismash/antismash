# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.test.helpers import (
    DummyProtocluster,
    DummyRecord,
    DummySubRegion,
)
from antismash.detection.subclusters.subregions import (
    LABEL,
    SubRegionMode,
    build_subregions,
    gather_foreign_areas,
)

RECORD_LENGTH = 1000
TOOL = "test_tool"


def create_protocluster(start, end, neighbourhood=0, product="test-product"):
    """ Builds a protocluster with a core covering exactly the given coordinates,
        extended on each side by the given neighbourhood

        A start greater than the end is treated as a core crossing the origin,
        otherwise the neighbourhood is clipped to the record rather than wrapping
    """
    if start > end:  # crossing the origin, so let the record length handle the wrapping
        return DummyProtocluster(core_start=start, core_end=end, product=product,
                                 neighbourhood_range=neighbourhood, record_length=RECORD_LENGTH)
    return DummyProtocluster(start=max(0, start - neighbourhood),
                             end=min(RECORD_LENGTH, end + neighbourhood),
                             core_start=start, core_end=end, product=product,
                             neighbourhood_range=neighbourhood)


def create_record(protoclusters=None, subregions=None, circular=False):
    """ Builds a record containing the given areas of other detection modules """
    record = DummyRecord(seq="A" * RECORD_LENGTH, circular=circular)
    for protocluster in protoclusters or []:
        record.add_protocluster(protocluster)
    for subregion in subregions or []:
        record.add_subregion(subregion)
    return record


def build_for_mode(mode, areas, foreign=None, foreign_subregions=None, *, circular=False):
    """ Builds subregions for the given subcluster areas in the given mode, from a
        record containing the given areas of other detection modules
    """
    record = create_record(protoclusters=foreign, subregions=foreign_subregions, circular=circular)
    return build_subregions(record, areas, tool=TOOL, mode=mode)


def get_coordinates(subregions):
    """ Converts subregions into a list of coordinate tuples for simpler comparisons """
    return [(int(sub.start), int(sub.end)) for sub in subregions]


class TestGatherForeignAreas(unittest.TestCase):
    def test_empty(self):
        assert gather_foreign_areas(create_record()) == []

    def test_protoclusters_only(self):
        protoclusters = [create_protocluster(50, 100), create_protocluster(200, 250)]
        record = create_record(protoclusters=protoclusters)
        assert gather_foreign_areas(record) == protoclusters

    def test_subregions_only(self):
        subregions = [DummySubRegion(10, 50), DummySubRegion(100, 150)]
        record = create_record(subregions=subregions)
        assert gather_foreign_areas(record) == subregions

    def test_both(self):
        protoclusters = [create_protocluster(50, 100, 10)]
        subregions = [DummySubRegion(10, 50)]
        record = create_record(protoclusters=protoclusters, subregions=subregions)
        assert gather_foreign_areas(record) == protoclusters + subregions


class TestBuildSubregionsGeneral(unittest.TestCase):
    def test_bad_modes(self):
        record = create_record()
        for bad in ["invalid", "", None, 5]:
            with self.assertRaisesRegex(ValueError, "Unknown subcluster region mode"):
                build_subregions(record, [create_protocluster(10, 50)], tool=TOOL, mode=bad)

    def test_no_areas(self):
        record = create_record(protoclusters=[create_protocluster(50, 100, 10)])
        areas = []
        for mode in SubRegionMode:
            assert build_subregions(record, areas, tool=TOOL, mode=mode) == []

    def test_tool_and_label(self):
        record = create_record(protoclusters=[create_protocluster(50, 100, 10)])
        areas = [create_protocluster(20, 80)]
        for mode in SubRegionMode:
            subregions = build_subregions(record, areas, tool=TOOL, mode=mode)
            assert subregions, mode
            for subregion in subregions:
                assert subregion.tool == TOOL
                assert subregion.label == LABEL


class TestCreateMode(unittest.TestCase):
    mode = SubRegionMode.CREATE

    def test_non_overlapping_areas(self):
        areas = [create_protocluster(10, 50), create_protocluster(100, 150)]
        assert get_coordinates(build_for_mode(self.mode, areas)) == [(10, 50), (100, 150)]

    def test_overlapping_areas_merge(self):
        areas = [create_protocluster(10, 50), create_protocluster(40, 80),
                 create_protocluster(150, 170)]
        assert get_coordinates(build_for_mode(self.mode, areas)) == [(10, 80), (150, 170)]


class TestExtendMode(unittest.TestCase):
    mode = SubRegionMode.EXTEND

    def test_without_foreign_areas(self):
        assert build_for_mode(self.mode, [create_protocluster(10, 50)]) == []

    def test_non_overlapping_discarded(self):
        areas = [create_protocluster(60, 100)]
        assert build_for_mode(self.mode, areas, foreign=[create_protocluster(100, 200)]) == []
        assert build_for_mode(self.mode, areas, foreign_subregions=[DummySubRegion(100, 200)]) == []

    def test_overlapping_kept_in_full(self):
        areas = [create_protocluster(60, 101)]
        subregions = build_for_mode(self.mode, areas, foreign=[create_protocluster(100, 200)])
        assert get_coordinates(subregions) == [(60, 101)]
        subregions = build_for_mode(self.mode, areas, foreign_subregions=[DummySubRegion(100, 200)])
        assert get_coordinates(subregions) == [(60, 101)]

    def test_only_overlapping_areas_used(self):
        # the two subclusters overlap each other, but only the second reaches
        # the foreign area, so the first must not be dragged in with it
        areas = [create_protocluster(10, 50), create_protocluster(40, 120)]
        subregions = build_for_mode(self.mode, areas, foreign=[create_protocluster(100, 200)])
        assert get_coordinates(subregions) == [(40, 120)]

    def test_neighbourhood_overlap_is_enough(self):
        # only the neighbourhoods overlap
        foreign = [create_protocluster(300, 350, neighbourhood=60)]
        subregions = build_for_mode(self.mode, [create_protocluster(150, 250)], foreign=foreign)
        assert get_coordinates(subregions) == [(150, 250)]


class TestClipMode(unittest.TestCase):
    mode = SubRegionMode.CLIP

    def test_without_foreign_areas(self):
        assert build_for_mode(self.mode, [create_protocluster(10, 50)]) == []

    def test_non_overlapping_discarded(self):
        foreign = [create_protocluster(100, 200)]
        assert build_for_mode(self.mode, [create_protocluster(10, 50)], foreign=foreign) == []

    def test_truncated_at_start(self):
        foreign = [create_protocluster(100, 200)]
        subregions = build_for_mode(self.mode, [create_protocluster(50, 150)], foreign=foreign)
        assert get_coordinates(subregions) == [(100, 150)]

    def test_truncated_at_end(self):
        foreign = [create_protocluster(100, 200)]
        subregions = build_for_mode(self.mode, [create_protocluster(150, 250)], foreign=foreign)
        assert get_coordinates(subregions) == [(150, 200)]

    def test_truncated_both_ends(self):
        # only the section shared with the foreign area survives
        foreign = [create_protocluster(120, 150)]
        subregions = build_for_mode(self.mode, [create_protocluster(100, 200)], foreign=foreign)
        assert get_coordinates(subregions) == [(120, 150)]

    def test_contained_subcluster_kept_in_full(self):
        foreign = [create_protocluster(100, 200)]
        subregions = build_for_mode(self.mode, [create_protocluster(120, 180)], foreign=foreign)
        assert get_coordinates(subregions) == [(120, 180)]

    def test_multiple_foreign_areas(self):
        # one subcluster spanning two separate foreign areas results in one
        # subregion for each of the shared sections
        subregions = build_for_mode(self.mode, [create_protocluster(100, 200)],
                                    foreign=[create_protocluster(50, 120)],
                                    foreign_subregions=[DummySubRegion(180, 260)])
        assert get_coordinates(subregions) == [(100, 120), (180, 200)]

    def test_single_base_overlap_discarded(self):
        # too small to be a useful area, so it's dropped
        foreign = [create_protocluster(99, 200)]
        assert build_for_mode(self.mode, [create_protocluster(10, 100)], foreign=foreign) == []

    def test_clipped_to_foreign_neighbourhood(self):
        # the limit is the full extent of the foreign area, including its neighbourhood
        foreign = [create_protocluster(280, 350, neighbourhood=100)]
        subregions = build_for_mode(self.mode, [create_protocluster(100, 200)], foreign=foreign)
        assert get_coordinates(subregions) == [(180, 200)]

    def test_multiple_areas(self):
        # make sure it works with multiple subcluster areas
        areas = [create_protocluster(10, 50), create_protocluster(40, 80)]
        subregions = build_for_mode(self.mode, areas, foreign=[create_protocluster(70, 120)])
        assert get_coordinates(subregions) == [(70, 80)]


class TestCrossOrigin(unittest.TestCase):
    # the mode varies with each test here, rather than being shared by the class

    def test_create(self):
        areas = [create_protocluster(950, 50)]
        subregions = build_for_mode(SubRegionMode.CREATE, areas, circular=True)
        assert len(subregions) == 1
        assert subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(950, 50)]

    def test_extend(self):
        areas = [create_protocluster(950, 50)]
        foreign = [create_protocluster(900, 970)]
        subregions = build_for_mode(SubRegionMode.EXTEND, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(950, 50)]

    def test_clip_both_crossing(self):
        # the foreign area ends before the subcluster
        areas = [create_protocluster(950, 50)]
        foreign = [create_protocluster(900, 20)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(950, 20)]

        # the foreign area starts after the subcluster
        areas = [create_protocluster(950, 50)]
        foreign = [create_protocluster(980, 100)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(980, 50)]

    def test_clip_subcluster_crossing(self):
        # only the part of the subcluster after the origin overlaps the foreign area
        areas = [create_protocluster(950, 50)]
        foreign = [create_protocluster(10, 100)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert not subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(10, 50)]

        # only the part of the subcluster before the origin overlaps the foreign area
        foreign = [create_protocluster(900, 970)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert not subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(950, 970)]

    def test_clip_foreign_crossing(self):
        # the subcluster starts after the origin, inside the end of the foreign area
        foreign = [create_protocluster(900, 50)]
        areas = [create_protocluster(10, 100)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert not subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(10, 50)]

        # the subcluster ends before the origin, entirely within the foreign area
        areas = [create_protocluster(920, 990)]
        subregions = build_for_mode(SubRegionMode.CLIP, areas, foreign=foreign, circular=True)
        assert len(subregions) == 1
        assert not subregions[0].crosses_origin()
        assert get_coordinates(subregions) == [(920, 990)]
