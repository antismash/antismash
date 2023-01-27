# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

from collections import defaultdict
import unittest

from antismash.modules.nrps_pks import name_mappings


class TestNameMappings(unittest.TestCase):
    def setUp(self) -> None:
        self.known = name_mappings.KNOWN_SUBSTRATES

    def test_long_name_unique(self) -> None:
        long_names = set()
        counts: dict[str, int] = defaultdict(int)

        for substrate in self.known:
            long_names.add(substrate.long)
            counts[substrate.short] += 1

        assert len(self.known) == len(long_names), counts

    def test_short_name_unique(self) -> None:
        short_names = set()
        counts: dict[str, int] = defaultdict(int)

        for substrate in self.known:
            short_names.add(substrate.short)
            counts[substrate.short] += 1

        assert len(self.known) == len(short_names), [(name, count) for name, count in counts.items() if count > 1]
