# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import json
import unittest

from helperlibs.bio import seqio

import antismash
from antismash.common import secmet
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.test import helpers
from antismash.config import build_config, destroy_config, get_config, update_config
from antismash.detection.genefunctions import core, smcogs, prepare_data


class TestSMCOGs(unittest.TestCase):
    def setUp(self):
        options = build_config(self.get_args(), isolated=True,
                               modules=antismash.get_all_modules())
        self.old_config = get_config().__dict__
        self.options = update_config(options)

        self.record = self.build_record(helpers.get_path_to_nisin_with_detection())

        prepare_data()

    def tearDown(self):
        destroy_config()
        update_config(self.old_config)

    def get_args(self):
        return ["--minimal", "--enable-genefunctions"]

    def build_record(self, genbank):
        # construct a working record
        with open(genbank, encoding="utf-8") as handle:
            seq_record = seqio.read(handle, "genbank")
        record = secmet.Record.from_biopython(seq_record, taxon="bacteria")
        assert record.get_protoclusters()
        assert record.get_protocluster(0).cds_children
        return record

    def test_classification_with_colon(self):
        # since SMCOG id and description are stored in a string separated by :,
        # ensure that descriptions containing : are properly handled
        # test gene is AQF52_5530 from CP013129.1
        translation = ("MDTHQREEDPVAARRDRTHYLYLAVIGAVLLGIAVGFLAPGVAVELKPLGTGFVN"
                       "LIKMMISPIIFCTIVLGVGSVRKAAKVGAVGGLALGYFLVMSTVALAIGLLVGNL"
                       "LEPGSGLHLTKEIAEAGAKQAEGGGESTPDFLLGIIPTTFVSAFTEGEVLQTLLV"
                       "ALLAGFALQAMGAAGEPVLRGIGHIQRLVFRILGMIMWVAPVGAFGAIAAVVGAT"
                       "GAAALKSLAVIMIGFYLTCGLFVFVVLGAVLRLVAGINIWTLLRYLGREFLLILS"
                       "TSSSESALPRLIAKMEHLGVSKPVVGITVPTGYSFNLDGTAIYLTMASLFVAEAM"
                       "GDPLSIGEQISLLVFMIIASKGAAGVTGAGLATLAGGLQSHRPELVDGVGLIVGI"
                       "DRFMSEARALTNFAGNAVATVLVGTWTKEIDKARVTEVLAGNIPFDEKTLVDDHA"
                       "PVPVPDQRAEGGEEKARAGV")
        cds = helpers.DummyCDS(0, len(translation))
        cds.translation = translation
        results = smcogs.classify("test", [cds], get_config())
        assert results.best_hits[cds.get_name()].hit_id == "SMCOG1212:sodium:dicarboxylate symporter"
        record = helpers.DummyRecord(seq=translation)
        record.add_cds_feature(cds)
        record.add_protocluster(helpers.DummyProtocluster(0, len(translation)))

        # if we don't handle multiple semicolons right, this line will crash
        results.add_to_record(record)
        gene_functions = cds.gene_functions.get_by_tool("smcogs")
        assert len(gene_functions) == 1
        assert str(gene_functions[0]).startswith("transport (smcogs) SMCOG1212:sodium:dicarboxylate symporter"
                                                 " (Score: 416; E-value: 2.3e-126)")

    def test_results_reconstruction(self):
        results = smcogs.classify(self.record.id, self.record.get_cds_features(), self.options)
        assert results.tool == "smcogs"
        assert results.best_hits["nisB"].hit_id == 'SMCOG1155:Lantibiotic dehydratase domain protein'

        data = json.loads(json.dumps(results.to_json()))
        reconstructed = core.FunctionResults.from_json(data, self.record)
        assert reconstructed.tool == "smcogs"
        assert reconstructed.best_hits["nisB"].hit_id == 'SMCOG1155:Lantibiotic dehydratase domain protein'
        assert reconstructed.to_json() == data

    def test_annotations(self):
        results = smcogs.classify(self.record.id, self.record.get_cds_features(), self.options)
        results.add_to_record(self.record)

        for cds in self.record.get_cds_features():
            if cds.gene_functions.get_by_tool("rule-based-clusters"):
                continue
            assert cds.gene_function == results.function_mapping.get(cds.get_name(), GeneFunction.OTHER)
