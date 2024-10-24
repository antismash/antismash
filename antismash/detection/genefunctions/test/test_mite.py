# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from unittest.mock import Mock, patch

from antismash.common.secmet.qualifiers import ECGroup, GeneFunction
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.common.utils import Hit
from antismash.detection.genefunctions.tools import mite


def build_dummy_hit(query, reference, identity=95, evalue=1e-15, bitscore=50.):
    return Hit(query_id=query, reference_id=reference, identity=identity,
               evalue=evalue, bitscore=bitscore)

def build_mite_hit_from_blast(blast, subfunctions=None, description=""):
    return mite.MiteHit(reference_id=blast.reference_id, query_id=blast.query_id,
                        identity=blast.identity, evalue=blast.evalue, bitscore=blast.bitscore,
                        subfunctions=subfunctions or [], description=description)

def test_classify():
    features = [DummyCDS(), DummyCDS(), DummyCDS()]
    hits = [
        build_dummy_hit(features[0].get_name(), "r2", identity=75.),
        build_dummy_hit(features[0].get_name(), "r1", identity=80.),
        build_dummy_hit(features[1].get_name(), "r2", identity=70.),
    ]
    entries = {
        "r1": mite.MiteEntry("ACC_1", "desc", GeneFunction.ADDITIONAL, groups=(ECGroup.LYASES,),
                             subfunctions=("func A",), version="1.many"),
        "r2": mite.MiteEntry("ACC_2", "desc", GeneFunction.ADDITIONAL, groups=(ECGroup.TRANSFERASES,),
                             subfunctions=("func B",), version="1.many"),
    }
    dataset = mite.Dataset(version="version", database_path="path", entries=entries,
                           url="dummy_url/{accession}")
    with patch.object(mite, "get_dataset", return_value=dataset):
        with patch.object(mite, "_get_blast_hits", return_value=hits) as patched_hits:
            result = mite.classify(features, Mock())
            patched_hits.assert_called_once()
    r1 = entries["r1"]
    r2 = entries["r2"]
    assert result.best_hits == {
        features[0].get_name(): build_mite_hit_from_blast(hits[1], subfunctions=list(r1.subfunctions),
                                                          description=", ".join(r1.subfunctions)),
        features[1].get_name(): build_mite_hit_from_blast(hits[2], subfunctions=list(r2.subfunctions),
                                                          description=", ".join(r2.subfunctions)),
    }
    assert all(func == GeneFunction.ADDITIONAL for func in result.function_mapping.values())
    assert result.group_mapping[features[0].get_name()] == entries["r1"].groups
