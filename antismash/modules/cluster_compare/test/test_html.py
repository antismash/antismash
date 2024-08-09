# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from antismash.common.secmet.locations import FeatureLocation
from antismash.common.test.helpers import (
    DummyCandidateCluster,
    DummyCDS,
    DummyProtocluster,
    DummyRecord,
    DummyRegion,
)
from antismash.modules.cluster_compare import html_output as html
from antismash.modules.cluster_compare.data_structures import (
    Hit,
    Mode,
    ReferenceCDS,
    ReferenceProtocluster,
    ReferenceRegion,
    ReferenceScorer,
)
from antismash.modules.cluster_compare.results import (
    ClusterCompareResults,
    DatabaseResults,
    ProtoToProtoScores,
    VariantResults,
)

from .test_analysis import DummyHit


def dummy_results(region, ref_accession):
    query_cdses = region.cds_children
    ref_cdses = [ReferenceCDS(name, function="biosynthetic", components={}, location=loc)
                 for name, loc in [
                    ("A", FeatureLocation(1, 15, 1)),
                    ("B", FeatureLocation(35, 40, 1)),
                 ]]
    cds_mapping = {str(i): ref.name for i, ref in enumerate(ref_cdses)}
    ref_by_name = {ref.name: ref for ref in ref_cdses}
    ref_proto = ReferenceProtocluster(accession=ref_accession, start=0, end=50, cds_mapping=cds_mapping,
                                      cdses=ref_by_name, cores=ref_cdses, product="prod")
    ref_region = ReferenceRegion(ref_accession, 0, 50, protoclusters=[ref_proto], cdses=ref_by_name,
                                 products=["test_product"], cds_mapping={"?": "?"}, description="desc", organism="org")
    hits: dict[ReferenceRegion, dict[str, Hit]] = {
        ref_region: {
            "B": DummyHit(ref_name="B", cds=query_cdses[1], ref_acc=ref_accession),
            "A": DummyHit(ref_name="A", cds=query_cdses[0], ref_acc=ref_accession),
        }
    }
    ref_scorer = ReferenceScorer(best_hits=hits[ref_region], reference=ref_region, identity=0.7,
                                 order=0.5, component=1.0, mode=Mode.BEST)
    by_region = [(ref_region, 0.5)]
    assert ref_region.end
    details = ProtoToProtoScores({region.get_region_number(): {ref_region: {ref_proto: ref_scorer}}})
    variant = VariantResults("dummy_variant", by_region, details, hits)
    db = DatabaseResults("dummy_db", url="https://some.url/{accession}/index.html",
                         description="desc", by_region={region.get_region_number(): {ref_accession: variant}})
    return ClusterCompareResults(ref_accession, {"dummy_db": db})


class TestJavascriptData(unittest.TestCase):
    def setUp(self):
        self.record = DummyRecord(seq="A" * 100)
        self.cdses = [
            DummyCDS(start=10, end=20, locus_tag="cds1"),
            DummyCDS(start=30, end=40, locus_tag="cds2"),
        ]
        self.protocluster = DummyProtocluster(start=10, end=40)
        self.candidate = DummyCandidateCluster(clusters=[self.protocluster])
        self.region = DummyRegion(candidate_clusters=[self.candidate])
        self.record.add_region(self.region)
        for cds in self.cdses:
            self.record.add_cds_feature(cds)
            self.region.add_cds(cds)
        self.func = html.generate_javascript_data
        self.ref_acc = "ref_accession5"

    def check_base_conversion(self, results, js_data):
        # make sure the data essentially matches what's expected
        assert list(js_data) == list(results.by_database)
        json_ref_acc = list(results.by_database["dummy_db"].by_region[self.region.get_region_number()])[0]
        assert json_ref_acc == self.ref_acc
        db = js_data["dummy_db"]
        expected_for_region = results.by_database["dummy_db"].by_region[self.region.get_region_number()]
        json_ref_clusters = db[self.ref_acc]["reference_clusters"]
        assert len(json_ref_clusters) == len(expected_for_region)
        expected_ref_region = expected_for_region[self.ref_acc].scores_by_region[0][0]
        ref_region_id = expected_ref_region.get_identifier()
        assert list(json_ref_clusters) == [ref_region_id]
        json_ref_cluster = json_ref_clusters[ref_region_id]
        assert json_ref_cluster["start"] == expected_ref_region.start
        assert json_ref_cluster["end"] == expected_ref_region.end
        expected_ref_hits = expected_for_region[self.ref_acc].hits_by_region
        assert len(json_ref_cluster["links"]) == len(expected_ref_hits[expected_ref_region])
        return json_ref_cluster["links"]

    def test_linear(self):
        results = dummy_results(self.region, self.ref_acc)
        # ensure some results exist
        assert results.by_database.get("dummy_db").by_region

        # convert to the JS-ready JSON
        js_data = self.func(self.record, self.region, results)
        # make sure the data essentially matches what's expected for the common case
        links = self.check_base_conversion(results, js_data)
        # ensure the calculated midpoint of each CDS is within the CDS and the record
        for link in links:
            cds = self.record.get_cds_by_name(link["query"])
            assert cds.location.start < link["query_loc"] < cds.location.end
            assert 0 < link["query_loc"] < len(self.record)

    def test_circular(self):
        assert list(self.record.all_features)
        self.record.rotate(cut_point=(self.cdses[0].location.end + self.cdses[1].location.start) // 2, padding=50)
        # before getting more complicated, make sure the region is cross-origin
        assert len(self.record.get_regions()[0].location.parts) == 2

        results = dummy_results(self.region, self.ref_acc)

        # convert to the JS-ready JSON
        js_data = self.func(self.record, self.region, results)
        # make sure the data essentially matches what's expected for the common case
        links = self.check_base_conversion(results, js_data)

        def check_js_cds(link, expected_name, extended=False):
            assert link["query"] == expected_name
            cds = self.record.get_cds_by_name(link["query"])
            mid_point = (cds.location.start + cds.location.end) // 2
            if extended:
                mid_point += len(self.record)
            js_midpoint = link["query_loc"]
            assert js_midpoint == mid_point

        # now that cds1 is after cds2, around the origin, make sure that it's extended
        check_js_cds(links[0], "cds2", extended=True)
        check_js_cds(links[1], "cds1", extended=False)
