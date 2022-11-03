# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import random
import unittest

from antismash.modules.nrps_pks.kr_analysis import kr_analysis


class TestActivity(unittest.TestCase):
    def test_is_active(self):
        # test all known actives
        good = set()
        for sig in ["E_HH", "K_YN", "K_YG"]:
            for char in "SAG":
                good.add(sig.replace("_", char))
        for sig in good:
            assert kr_analysis.is_active(sig)

        # a little fuzzing for inactives
        valid_chars = list("ARNDCQEGHILKMFPSTWYV")
        random.shuffle(valid_chars)
        chunks = "".join(valid_chars)
        for i in range(len(chunks) - 3):
            sig = chunks[i:i+4]
            assert sig in good or not kr_analysis.is_active(sig)


class TestStereochemistry(unittest.TestCase):
    def test_predicition(self):
        test_func = kr_analysis.predict_stereochemistry

        assert test_func("XXXWXYAN") == "A1"
        assert test_func("XXXWHYAN") == "A2"
        assert test_func("LDDWHYXN") == "B1"
        assert test_func("LDDWHYPN") == "B2"
        assert test_func("XXXXXXXX") == "C1"
        assert test_func("XXXXXYXX") == "C2"
        assert test_func("XXXXXYXN") is None
