# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from dataclasses import FrozenInstanceError
import json
import unittest

from antismash.common.hmmer import HmmerHit


def create_hmmer_hit(location="[500:700]", label="ref_name", locus_tag="locus",
                     domain="hit", evalue=1e-10, score=40.5, translation=None,
                     identifier="TESTID_0005", description="some description",
                     protein_start=9, protein_end=None):
    if translation is None:
        if protein_end is None:
            protein_end = 13
        translation = "M" * (protein_end - protein_start)

    protein_end = protein_end or (protein_start + len(translation))

    return HmmerHit(
        location=location,
        label=label,
        locus_tag=locus_tag,
        domain=domain,
        evalue=evalue,
        score=score,
        translation=translation,
        identifier=identifier,
        description=description,
        protein_start=protein_start,
        protein_end=protein_end
    )


class TestHmmerHit(unittest.TestCase):
    def test_immutable(self):
        hit = create_hmmer_hit()
        with self.assertRaisesRegex(FrozenInstanceError, "cannot assign to field 'invalid'"):
            hit.invalid = "test"
        assert hit.protein_start
        with self.assertRaisesRegex(FrozenInstanceError, "cannot assign to field 'protein_start'"):
            hit.protein_start = 7

    def test_length(self):
        hit = create_hmmer_hit(protein_start=9, protein_end=13)
        assert len(hit) == 4
        hit = create_hmmer_hit(protein_start=1, protein_end=23)
        assert len(hit) == 22

    def test_invalid_positions(self):
        with self.assertRaisesRegex(ValueError, "inverted start and end"):
            create_hmmer_hit(protein_start=9, protein_end=9)
        with self.assertRaisesRegex(ValueError, "inverted start and end"):
            create_hmmer_hit(protein_start=5, protein_end=3)

    def test_mismatching_lengths(self):
        hit = create_hmmer_hit(protein_start=5, protein_end=15, translation="A"*10)
        assert len(hit) == len(hit.translation) == hit.protein_end - hit.protein_start
        with self.assertRaisesRegex(ValueError, "translation length does not match"):
            create_hmmer_hit(protein_start=5, protein_end=15, translation="A"*4)

    def test_json_conversion(self):
        hit = create_hmmer_hit()
        as_json = json.loads(json.dumps(hit.to_json()))
        assert hit == HmmerHit.from_json(as_json)
