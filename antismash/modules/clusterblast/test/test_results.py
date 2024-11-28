# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import patch

from antismash.modules.clusterblast.data_structures import Protein, ReferenceCluster
from antismash.modules.clusterblast import results
from antismash.common.secmet.test.helpers import DummyRecord, DummyRegion


def create_reference(accession, proteins):
    return ReferenceCluster(
        accession,
        "c1",
        proteins,
        description="desc",
        cluster_type="type",
        tags=[protein.locus_tag for protein in proteins],
    )


def create_protein(name, tag):
    return Protein(name, tag, "1234-5678", "1", "annots")


class TestRegionResults(unittest.TestCase):
    def test_reference_protein_filtering(self):
        proteins = {f"t{c}": create_protein(c, f"t{c}") for c in "ABCDEFX"}
        scores = []
        for i in range(4, 1, -1):
            scores.append(results.Score())
            scores[-1].hits = i

        ranking = [
            (create_reference("same_accession", proteins=[proteins["tX"]]), scores[0]),
            (create_reference("acc1", proteins=[proteins["tA"]]), scores[0]),
            (create_reference("acc2", proteins=[proteins[f"t{c}"] for c in "BC"]), scores[2]),
            (create_reference("acc3", proteins=[proteins[f"t{c}"] for c in "DEF"]), scores[2]),
        ]

        region = DummyRegion()
        region.parent_record = DummyRecord()
        region.parent_record.id = "same_accession"
        with patch.object(results, "get_display_limit", return_value=2):
            result = results.RegionResult(region, ranking, proteins, "")
        # ensure self hit is removed entirely
        assert len(result.ranking) == len(ranking) - 1
        assert ranking[0] not in result.ranking
        assert len(proteins) == 7
        assert len(result.displayed_reference_proteins) == 3
