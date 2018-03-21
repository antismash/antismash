# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

from antismash.common import hmmer
from antismash.common.test.helpers import FakeHSP, DummyRecord, DummyCDS


def _create_dummy_record(reverse=False):
    seq = Seq('GTGGAGCGGTACTAAATGTACTCCACTATCTGCTGATTGGAAACCACGGAGCGCTCTTAG',
              generic_dna)
    strand = 1
    if reverse:
        seq = seq.reverse_complement()
        strand = -1
    rec = DummyRecord(seq=str(seq))

    idx = 1
    for start, end in [(0, 15), (15, 36), (36, 60)]:
        if reverse:
            start, end = len(seq) - end + 3, len(seq) - start  # TODO: check this
        rec.add_cds_feature(DummyCDS(start, end, strand=strand,
                            locus_tag="orf%04d" % idx))
        idx += 1

    return rec
