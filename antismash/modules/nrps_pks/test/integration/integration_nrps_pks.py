# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
import unittest

from helperlibs.wrappers.io import TemporaryDirectory
import minimock

import antismash
from antismash.common import path
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config
from antismash.modules import nrps_pks

class IntegrationWithoutNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal"],
                                    isolated=True, modules=antismash.get_all_modules())
        assert not nrps_pks.is_enabled(self.options)
        self.tracker = minimock.TraceTracker()
        minimock.mock("nrps_pks.run_on_record", tracker=self.tracker)

    def tearDown(self):
        destroy_config()
        minimock.restore()

    def test_minimal(self):
        with TemporaryDirectory(change=True):
            antismash.main.run_antismash(helpers.get_path_to_balhymicin_genbank(), self.options)
        # make sure it didn't run
        minimock.assert_same_trace(self.tracker, "")

class IntegrationNRPSPKS(unittest.TestCase):
    def setUp(self):
        self.options = build_config(["--minimal", "--enable-nrps-pks"],
                                    isolated=True, modules=antismash.get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_balhymicin(self):
        with TemporaryDirectory(change=True):
            assert not os.path.exists("structures/genecluster1.png")
            antismash.main.run_antismash(helpers.get_path_to_balhymicin_genbank(), self.options)
            assert os.path.exists("structures/genecluster1.png")

    def test_CP002271_c19(self):
        with TemporaryDirectory(change=True):
            assert not os.path.exists("structures/genecluster1.png")
            antismash.main.run_antismash(path.get_full_path(__file__, 'data', 'CP002271.1.cluster019.gbk'), self.options)
            assert os.path.exists("structures/genecluster1.png")
