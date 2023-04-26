# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring,consider-using-with

from unittest.mock import patch

from antismash.common.subprocessing import memesuite

from .helpers import DummyResult

# output format of versions up to 4.11.2
DATA_4_11_2 = """
#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence
a	gene1	1	12	+	24	5.62e-08	2.09e-07	acgtacgtacgt
b	gene9	1	12	+	24	5.62e-08	2.09e-07	acgtacgtacgt
"""

# versions from 4.11.3 to latest (5.5.2 at time of writing)
DATA_latest = """
motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
MA0060.1	NFYA	chr2	60221163	60221178	-	18.75	3.36e-09	0.00195	CTCGGCCAATCAGAGC
UP00093_1	Klf7_primary	chr2	172266018	172266033	+	19.4954	1.06e-08	0.00496	GCGACCCCGCCCCTTT
"""


def test_legacy():
    with patch.object(memesuite, "execute", return_value=DummyResult(DATA_4_11_2)):
        results = memesuite.run_fimo_simple("file", "queries")
    assert results == [
        memesuite.FIMOMotif("a", "", "gene1", 1, 12, "+", 24, 5.62e-08, 2.09e-07, "acgtacgtacgt"),
        memesuite.FIMOMotif("b", "", "gene9", 1, 12, "+", 24, 5.62e-08, 2.09e-07, "acgtacgtacgt"),
    ]
    assert results == memesuite.read_fimo_output(DATA_4_11_2)  # these should be equivalent


def test_latest():
    with patch.object(memesuite, "execute", return_value=DummyResult(DATA_latest)):
        results = memesuite.run_fimo_simple("file", "queries")
    assert results == [
        memesuite.FIMOMotif("MA0060.1", "NFYA", "chr2", 60221163, 60221178, "-",
                            18.75, 3.36e-09, 0.00195, "CTCGGCCAATCAGAGC"),
        memesuite.FIMOMotif("UP00093_1", "Klf7_primary", "chr2", 172266018, 172266033, "+",
                            19.4954, 1.06e-08, 0.00496, "GCGACCCCGCCCCTTT"),
    ]
    assert results == memesuite.read_fimo_output(DATA_latest)  # these should be equivalent
