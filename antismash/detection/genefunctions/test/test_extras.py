# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from unittest.mock import Mock, patch

from antismash.common.secmet.qualifiers import ECGroup, GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.detection.genefunctions.tools import extras


def test_load():
    entries = extras._load_metadata()
    for name, entry in entries.items():
        assert name == entry.identifier
        if entry.function == GeneFunction.ADDITIONAL:
            assert entry.groups


def build_dummy_entry(**kwargs):
    defaults = {
        "identifier": "dummy_id",
        "description": "desc",
        "cutoff": 25.,
        "function": GeneFunction.ADDITIONAL,
        "groups": [ECGroup.LYASES,],
        "source": "doi...",
        "subfunctions": ["some func"],
    }
    defaults.update(kwargs)
    return extras.Entry(**defaults)


def test_classify():
    hits = {
        "q1": extras.HMMHit(query_id="q1", reference_id="r1", evalue=1e-20, bitscore=50),
        "q2": extras.HMMHit(query_id="q2", reference_id="r2", evalue=1e-20, bitscore=250),
        "q3": extras.HMMHit(query_id="q2", reference_id="r2", evalue=1e-20, bitscore=150),
    }
    entries = {
        "r1": build_dummy_entry(identifier="r1", cutoff=20., function=GeneFunction.REGULATORY),
        "r2": build_dummy_entry(identifier="r2", cutoff=200., function=GeneFunction.ADDITIONAL),
    }
    with patch.object(extras, "_load_metadata", return_value=entries):
        with patch.object(extras, "scan_profiles_for_functions", return_value=hits) as patched_scan:
            results = extras.classify([DummyCDS(locus_tag="q1"), DummyCDS(locus_tag="q2")], Mock())
            patched_scan.assert_called_once()
    assert len(results.best_hits) == 2
    assert results.best_hits["q1"].reference_id == "r1"
    assert results.best_hits["q2"].reference_id == "r2"
    assert results.function_mapping == {
        "q1": entries["r1"].function,
        "q2": entries["r2"].function,
    }
    assert results.group_mapping == {
        "q1": entries["r1"].groups,
        "q2": entries["r2"].groups,
    }
