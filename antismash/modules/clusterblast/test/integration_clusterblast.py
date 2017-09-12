# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import unittest

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.main import get_all_modules, detect_signature_genes
from antismash.common import deprecated
from antismash.common.module_results import ModuleResults
import antismash.common.test.helpers as helpers
from antismash.config import args
from antismash.modules import clusterblast

class Base(unittest.TestCase):
    def setUp(self):
        options = args.build_parser(modules=get_all_modules()).parse_args(self.get_args())
        self.old_config = args.Config().__dict__
        self.options = args.Config(options)

        assert clusterblast.check_prereqs() == []
        assert clusterblast.check_options(self.options) == []
        assert clusterblast.is_enabled(self.options)

        self.record = None # set it or build it with build_record(genbank)

    def get_args(self):
        """ override with the args you'll need to use in setUp(),
            format is as on the commandline, e.g. ["--tta", "--minimal"]
        """
        self.fail() # not overridden

    def build_record(self, genbank):
        # construct a working record
        self.record = deprecated.parse_input_sequence(genbank, self.options)[0]
        detect_signature_genes(self.record, self.options)
        clusters = self.record.get_clusters()
        # make sure it's worth using
        assert clusters
        for cluster in clusters:
            assert cluster.cds_children
        assert deprecated.get_cds_features_within_clusters(self.record)

    def tearDown(self):
        args.Config({})

    def get_results(self):
        """ override with a function that runs *blast, verifies basics and
            returns results instance for further testing
        """
        self.fail() # wasn't overriden

    def check_nisin(self, expected):
        with TemporaryDirectory(change=True):
            self.build_record(helpers.get_path_to_nisin_genbank())
            results = self.get_results()
            assert len(results.cluster_results) == 1
            cluster = results.cluster_results[0]
            if expected > clusterblast.get_result_limit():
                assert cluster.total_hits == expected
                expected = clusterblast.get_result_limit()
            assert len(cluster.ranking) == expected # will change if database does
        return results

    def check_balhymicin(self, expected):
        with TemporaryDirectory(change=True):
            self.build_record(helpers.get_path_to_balhymicin_genbank())
            results = self.get_results()
            assert len(results.cluster_results) == 1
            cluster = results.cluster_results[0]
            if expected > clusterblast.get_result_limit():
                assert cluster.total_hits == expected
                expected = clusterblast.get_result_limit()
            assert len(cluster.ranking) == expected # will change if database does
        return results

# TODO: test with a small sequence instead (grab a CDS that hit and take it's translation)

class GeneralIntegrationTest(Base):
    def get_args(self):
        return ["--cb-general"]

    def get_results(self):
        results = clusterblast.run_on_record(self.record, self.options)
        assert isinstance(results, ModuleResults)
        assert results.general
        assert results.general.search_type == "clusterblast"
        assert results.knowncluster is None
        assert results.subcluster is None
        json = results.to_json()
        assert json
        new = clusterblast.ClusterBlastResults.from_json(json, self.record)
        assert json == new.to_json()
        return results.general

    def test_nisin(self):
        self.check_nisin(2452)

class KnownIntegrationTest(Base):
    def get_args(self):
        return ["--cb-knowncluster"]

    def get_results(self):
        results = clusterblast.run_on_record(self.record, self.options)
        assert isinstance(results, ModuleResults)
        assert results.general is None
        assert results.knowncluster
        assert results.knowncluster.search_type == "knownclusterblast"
        assert results.subcluster is None
        json = results.to_json()
        assert json
        new = clusterblast.ClusterBlastResults.from_json(json, self.record)
        assert json == new.to_json()
        return results.knowncluster

    def test_nisin(self):
        self.check_nisin(7)

    def test_balhymicin(self):
        self.check_balhymicin(102)

class SubIntegrationTest(Base):
    def get_args(self):
        return ["--cb-subcluster"]

    def get_results(self):
        results = clusterblast.run_on_record(self.record, self.options)
        assert isinstance(results, ModuleResults)
        assert results.general is None
        assert results.knowncluster is None
        assert results.subcluster
        assert results.subcluster.search_type == "subclusterblast"
        json = results.to_json()
        assert json
        new = clusterblast.ClusterBlastResults.from_json(json, self.record)
        assert json == new.to_json()
        return results.subcluster

    def test_nisin(self):
        self.check_nisin(0)

    def test_balhymicin(self):
        self.check_balhymicin(21)
