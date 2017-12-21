# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from argparse import Namespace
import os.path  # mocked, pylint: disable=unused-import
import unittest

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from Bio import SearchIO  # mocked, pylint: disable=unused-import
from minimock import mock, restore, TraceTracker, assert_same_trace, Printer
from helperlibs.bio import seqio

from antismash.common import path  # mocked, pylint: disable=unused-import
from antismash.common import record_processing, subprocessing
from antismash.common.test import helpers
from antismash.modules import active_site_finder


class TestASF(unittest.TestCase):
    "Test Active Site Finder class"

    def setUp(self):
        "set-up required variables, etc; load demo sequence and parse XML"
        self.asf_path = os.path.dirname(os.path.dirname(__file__))

        self.options = Namespace()
        qualifier_tags = Namespace()
        qualifier_tags.ASF_scaffold = 'aSASF_scaffold'
        qualifier_tags.ASF_choice = 'aSASF_choice'
        qualifier_tags.ASF_prediction = 'aSASF_prediction'

        result = Namespace()
        result.id = "fullhmmer_oxyB_0001"
        result.hsps = [Namespace]
        result.hsps[0] = Namespace()
        result.hsps[0].query_start = 1
        result.hsps[0].query_end = 99
        result.hsps[0].hit_start = 322
        result.hsps[0].hit_end = 428
        result.hsps[0].aln = [Namespace(), Namespace()]
        result.hsps[0].aln[0].seq = "RAVDELIRYLTVPYGPTPRIAKQDVTVGDQVIKAGESVICSLPAANRDPALVPDADRLDVTR--------DPVPHVAFGHGIHHCLGAALARLELRTVFTALWRRF"
        result.hsps[0].aln[1].seq = "avikEtLRlhpvvplllpRevtkdvvirgylipkGtevivnlyalhrdpevfpnPeeFdpeRFldekgsrksfaflPFGaGpRnCiGerlArmelklflatlLqnF"

        self.result = result

        self.options.qualifier_tags = qualifier_tags

        self.record = record_processing.parse_input_sequence(path.get_full_path(__file__, 'Y16952.3.final.gbk'))[0]

        self.assertEqual(360, self.record.get_feature_count())
        assert len(self.record.get_pfam_domains()) == 130
        assert len(self.record.get_antismash_domains()) == 34

        xmltext = """<?xml version="1.0" encoding="UTF-8"?>

<resource xmlns:xsd="http://www.w3.org/2001/XMLSchema"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <analysis name='ASP_P450Oxy' type="active_site">
        <Prerequisite>
            <primary_tag_type>PFAM_domain</primary_tag_type>
            <tag>domain</tag>
            <tag_value>p450</tag_value>
        </Prerequisite>
        <Execute program="hmmscan" CaptureConsole="TRUE">
            <!-- Currently, the location of the hmmpfam2 binary and the database location
                is inferred from the antismash configuration file! -->

            <parameters>
                <!-- If prefixes for parameters are required they can be added as attribute
                    prefix -->
                <parameter name="evalue" prefix="--domE">0.1</parameter>
                <parameter name="cpus" prefix="--cpu">1</parameter>
            </parameters>
            <database>p450.hmm3</database>
            <db_source>PFAM24</db_source>
            <BioPythonParser>hmmer3-text</BioPythonParser>
        </Execute>
        <Alignment>
            <scaffold>
                <scaffoldOffset>327,330,400,403,409</scaffoldOffset>
                <scaffoldValue>E,R,F,G,G</scaffoldValue>
            </scaffold>
            <choice result="active site cystein present">
                <offset>407</offset>
                <value>C</value>
                <comment>Cytochrome P450 oxygenase active site cystein; coordinates heme Fe ligand</comment>
            </choice>
        </Alignment>
        <description>Pediction of cytochrome P450 active site cystein</description>
        <referenceList>
            <reference>Del Vecchio, F., H. Petkovic, S. G. Kendrew, L. Low, B.
                Wilkinson, R. Lill, J. Cortes, B. A. Rudd, J. Staunton, and P. F.
                Leadlay. 2003. Active-site residue, domain and module swaps in
                modular polyketide synthases. J Ind. Microbiol Biotechnol
                30:489-494.</reference>
        </referenceList>
    </analysis>
</resource>
"""

        etree_obj = ET.fromstring(xmltext)
        etree_obj.find('./analysis')
        self.etree_obj = etree_obj

        # now acutally generate ASF object
        asf = active_site_finder.active_site_finder(self.record, self.options)

        assert isinstance(asf, active_site_finder.active_site_finder)

        self.my_ASF = asf

    def test__init(self):
        "Test active_site_finder.__init__ method"

        self.tt = TraceTracker()

        class dummyET:
            def getroot(self):
                return dummyET()

            def findall(self, _dummy):
                return "ET-Subtree"

        mock('ET.parse', returns=dummyET(), tracker=self.tt)
        mock('active_site_finder.active_site_finder.check_prereqs', returns=[], tracker=self.tt)

        expected = """    Called ET.parse(
        '{0}/config/SignatureResources.xml')""".format(self.asf_path)

        active_site_finder.active_site_finder(self.record, self.options)

        assert_same_trace(self.tt, expected)

        # restore mocked objects that other tests can run with the original code
        restore()
        self.tt = TraceTracker()

        # now acutally generate ASF object
        asf = active_site_finder.active_site_finder(self.record, self.options)

        assert isinstance(asf, active_site_finder.active_site_finder)

        self.my_ASF = asf

    def test__get_scaffold_annotation(self):
        "Test active_site_finder._get_scaffold_annotation method"

        my_ScaffXML = self.etree_obj.find('./analysis/Alignment/scaffold')

        resultline = self.my_ASF._get_scaffold_annotation(self.result, my_ScaffXML)

        expectation = "Scaffold coordinates: (327,330,400,403,409); scaffold residues: (E,R,F,G,G); expected: (E,R,F,G,G); matchArray: (True,True,True,True,True); emission probability array (n.d.,n.d.,n.d.,n.d.,n.d.); overall match: TRUE"

        self.assertEqual(resultline, expectation, "scaffold line mismatch")

    def test_get_prediction_annotation(self):
        "Test active_site_finder._get_prediction_annotation method"
        my_ChoicesXML = self.etree_obj.findall('./analysis/Alignment/choice')

        (choiceList, predictionList) = active_site_finder.get_prediction_annotation(self.result, my_ChoicesXML)

        expected = ["Description: Cytochrome P450 oxygenase active site cystein; coordinates heme Fe ligand, choice result: active site cystein present, choice coordinates: (407); residues: (C); expected for choice: (C); matchArray: (True); emission probability array (n.d.); overall match: TRUE"]
        self.assertListEqual(choiceList, expected, "prediction choices mismatch")

        expected = ["Full match for prediction: active site cystein present"]
        self.assertListEqual(predictionList, expected, "prediction string mismatch")

    def test__execute_tool(self):
        "Test active_site_finder._execute_tool method"
        self.tt = TraceTracker()


        class NamePrinter(TraceTracker):
            def call(self, func_name, *args, **kw):
                msg = 'Called %s()\n' % func_name
                self.file.write(msg)

        self.tt2 = NamePrinter()
        dummy_run_result = subprocessing.RunResult('', "shell output".encode("utf-8"), "", 0, piped_out=True, piped_err=False)
        mock('subprocessing.execute', returns=dummy_run_result, tracker=self.tt)

        mock('SearchIO.parse', returns=["SearchIO object"], tracker=self.tt2)

        result = self.my_ASF._execute_tool(self.etree_obj.find('./analysis'), fileName="testTempfile")

        self.assertListEqual(result, ["SearchIO object"])

        expected = """        Called subprocessing.execute(
            ['hmmscan', '--domE', '0.1', '--cpu', '1', '{0}/data/p450.hmm3', 'testTempfile'])
      """.format(self.asf_path)
        assert_same_trace(self.tt, expected)

        expected = "   Called SearchIO.parse()"
        assert_same_trace(self.tt2, expected)

        restore()

        self.tt.clear()
        self.tt2.clear()
        mock('subprocessing.execute', returns=dummy_run_result, tracker=self.tt)

        mock('SearchIO.parse', returns=["SearchIO object"], tracker=self.tt2)

        result = self.my_ASF._execute_tool(self.etree_obj.find('./analysis'), stdin_data="fasta sequence from stdin")

        self.assertListEqual(result, ["SearchIO object"])

        expected = """Called subprocessing.execute(
    ['hmmscan', '--domE', '0.1', '--cpu', '1', '{0}/data/p450.hmm3'],
    stdin='fasta sequence from stdin')
      """.format(self.asf_path)
        assert_same_trace(self.tt, expected)

        expected = "Called SearchIO.parse()"
        assert_same_trace(self.tt2, expected)

        restore()
        self.tt = TraceTracker()
        del self.tt2

    def test__run_external_tool(self):
        "Test active_site_finder._run_external_tool method"

        self.tt = TraceTracker()

        targets = [pfam for pfam in self.record.get_pfam_domains() if pfam.domain == "p450"]
        assert len(targets) == 6

        mock('active_site_finder.active_site_finder._execute_tool', returns=["external program call successful"], tracker=self.tt)

        result = self.my_ASF._run_external_tool(self.etree_obj.find('./analysis'), targets)
        expected = ["external program call successful"]
        self.assertListEqual(result, expected)

    def test_fix_coordinates(self):
        "Test active_site_finder._fix_coordinates method"

        seq = "lsgpeavkevlikkgeefs..grgdeallatsrkafkgkgvlfangekwkklRrfltptltsf.klsleelveeeaedlv"

        coord = active_site_finder.fix_coordinates(10, seq)
        self.assertEqual(coord, 10, "assignment without gap")

        coord = active_site_finder.fix_coordinates(25, seq)
        self.assertEqual(coord, 27, "assignment with 2 gaps (at pos 20,21 of original string)")

        coord = active_site_finder.fix_coordinates(20, seq)
        self.assertEqual(coord, 22, "assignment with 2 gaps (at end position 20 and 21)")

        coord = active_site_finder.fix_coordinates(76, seq)
        self.assertEqual(coord, 79, "assignment with 3 gaps (at end position 20,21 and 64)")

        with self.assertRaises(ValueError):
            coord = active_site_finder.fix_coordinates(77, seq)

    def test_check_prereqs(self):
        "Test active_site_finder.check_prereqs method"
        # THIS TEST HAS TO BE UPDATED WHEN NEW PROFILES ARE ADDED TO THE MODULE!

        self.tt = TraceTracker()
        discard_tracker = TraceTracker()

        mock('path.locate_executable', returns="/my/path/to/executable", tracker=self.tt)
        mock('path.locate_file', returns="/my/path/to/file", tracker=self.tt)
        mock('os.path.getmtime', returns_iter=[2, 1, 4, 3, 6, 5, 8, 7], tracker=discard_tracker)

        result = self.my_ASF.check_prereqs()

        self.assertListEqual(result, [], "return empty list if executables/files are found")

        expected = """Called path.locate_executable('blastp')
Called path.locate_executable('hmmpfam2')
Called path.locate_executable('hmmscan')
Called path.locate_executable('hmmpress')
Called path.locate_file(
    '{0}/data/PKSI-KR.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-KS_N.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-KS_C.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-AT.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-ACP.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-DH.hmm2')
Called path.locate_file(
    '{0}/data/Thioesterase.hmm2')
Called path.locate_file(
    '{0}/data/PKSI-ER.hmm2')
Called path.locate_file(
    '{0}/data/aa-activating.aroundLys.hmm2')
Called path.locate_file(
    '{0}/data/p450.hmm3')
Called path.locate_file(
    '{0}/data/p450.hmm3.h3f')
Called path.locate_file(
    '{0}/data/p450.hmm3.h3i')
Called path.locate_file(
    '{0}/data/p450.hmm3.h3m')
Called path.locate_file(
    '{0}/data/p450.hmm3.h3p')""".format(self.asf_path)
        assert_same_trace(self.tt, expected)

        restore()
        self.tt = TraceTracker()

        mock('path.locate_executable', returns=None, tracker=self.tt)
        mock('path.locate_file', returns="/my/path/to/file", tracker=self.tt)
        mock('os.path.getmtime', returns_iter=[2, 1, 4, 3, 6, 5, 8, 7], tracker=discard_tracker)

        result = self.my_ASF.check_prereqs()
        expected = [
            "Failed to locate file: 'blastp'",
            "Failed to locate file: 'hmmpfam2'",
            "Failed to locate file: 'hmmscan'",
            "Failed to locate file: 'hmmpress'",
        ]
        self.assertListEqual(result, expected, "test result if file not found. %r vs %r" % (result, expected))
        restore()
