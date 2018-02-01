# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import os
from shutil import copy

from antismash.common import path
from antismash.detection.cassis import runners

from .test_cassis import create_fake_record, convert_newline, CassisTestCore


def read_generated_expected_file(generated_file, expected_file):
    """Read generated and expected files and save to string variables"""
    generated_string = ""
    with open(generated_file) as handle:
        generated_string = handle.read()
    generated_string = convert_newline(generated_string.rstrip())

    expected_string = ""
    with open(path.get_full_path(__file__, "data", expected_file)) as handle:
        expected_string = handle.read()
    expected_string = convert_newline(expected_string.rstrip())

    return [generated_string, expected_string]


class TestCassisRunners(CassisTestCore):
    def test_run_fimo(self):
        seq_record = create_fake_record()
        meme_dir = os.path.join(self.options.output_dir, "meme")
        fimo_dir = os.path.join(self.options.output_dir, "fimo")
        fimo_subdirs = ["+00_-03", "+03_-00"]

        # fimo input (i)
        # --> sequences to search in
        copy(
            path.get_full_path(__file__, "data", "expected_promoter_sequences.fasta"),
            os.path.join(self.options.output_dir, seq_record.name + "_promoter_sequences.fasta")
        )
        # fimo input (ii)
        # --> meme output (meme.html; motifs to search for)
        # --> file derived from meme output (binding_sites.fasta;
        #      only additional information for users; still checked by fimo)
        for subdir in fimo_subdirs:
            if not os.path.exists(os.path.join(meme_dir, subdir)):
                os.makedirs(os.path.join(meme_dir, subdir))
            copy(path.get_full_path(__file__, "data", "fake_meme.html"),
                 os.path.join(meme_dir, subdir, "meme.html"))
            copy(path.get_full_path(__file__, "data", "fake_binding_sites.fasta"),
                 os.path.join(meme_dir, subdir, "binding_sites.fasta"))

        runners.run_fimo(meme_dir, fimo_dir, seq_record, self.options, verbose=True)

        for subdir in fimo_subdirs:
            fimo_result, expected_fimo_result = read_generated_expected_file(
                                                    os.path.join(self.options.output_dir, "fimo", subdir, "fimo.txt"),
                                                    "fake_long_fimo.txt")

            self.assertEqual(fimo_result, expected_fimo_result)
