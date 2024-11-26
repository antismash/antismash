# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from unittest import TestCase

from antismash.config import destroy_config, update_config
from antismash.common.fasta import read_fasta
from antismash.common.record_processing import parse_input_sequence, pre_process_sequences
from antismash.common.secmet.test.helpers import DummyRecord
from antismash.common.test.helpers import get_simple_options, get_path_to_nisin_fasta
from antismash.support import genefinding

# for any sequence padding added, ensure it won't contain gene calls by
# having stop codons in both forward and reverse strands
STOP_PADDING = "TAGCTAACTGACCGTCAGTTAGCTA"


class TestProdigal(TestCase):
    def setUp(self):
        self.options = update_config(get_simple_options(genefinding,
                                                        ['--genefinding-tool', 'prodigal',
                                                         '--cpus', '1']))

    def tearDown(self):
        destroy_config()

    def test_nisin(self):
        record = parse_input_sequence(get_path_to_nisin_fasta())[0]
        assert record.get_feature_count() == 0
        record = pre_process_sequences([record], self.options, genefinding)[0]
        assert record.get_feature_count() == 12
        # and make sure they're all CDS features
        assert len(record.get_cds_features()) == 12

    def test_rotated_nis_b(self):
        # test that cross-origin prodigal works as expected by taking the nisB gene
        # from nisin, putting the origin of the record in the middle of it, with some filler around it
        # i.e.  |AABB| -> AA|BB -> |BB....AA|

        # start with finding the gene with prodigal, without any cross-origin complication
        sequence = list(read_fasta(get_path_to_nisin_fasta()).values())[0]
        old_record = DummyRecord(seq=sequence)
        genefinding.run_prodigal(old_record)
        original_features = old_record.get_cds_features()
        # nisB is the longest gene, so take the sequence of that, along with its stop codon
        nis_b = sorted(original_features, key=lambda x: len(x.translation), reverse=True)[0]
        nis_b_seq = old_record.seq[nis_b.location.start:nis_b.location.end + 3]
        assert len(nis_b.location) == 2982
        assert len(nis_b.translation) == 993

        # rotate and add padding to test the cross-origin part
        # the padding helps avoid the creation of any genes which may cause prodigal
        # to give different results
        padding = STOP_PADDING * 200
        cut_point = len(nis_b_seq) // 3
        new_seq = nis_b_seq[cut_point:] + padding + nis_b_seq[:cut_point]

        new_linear_record = DummyRecord(seq=new_seq)
        assert not new_linear_record.is_circular()
        genefinding.run_prodigal(new_linear_record)
        found_linear = new_linear_record.get_cds_features()
        assert len(found_linear) == 2

        circular_record = DummyRecord(seq=new_seq, circular=True)
        assert circular_record.is_circular()
        genefinding.run_prodigal(circular_record)
        found_circular = circular_record.get_cds_features()
        assert len(found_circular) == 1
        cross_origin_cds = found_circular[0]
        assert cross_origin_cds.crosses_origin()
        assert cross_origin_cds.start > cross_origin_cds.end

        assert cross_origin_cds.translation == nis_b.translation
        assert found_linear[1].location.start == cross_origin_cds.start
        assert found_linear[0].location.end == cross_origin_cds.end
