# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
from shutil import copy

from antismash.common import path
from antismash.common.subprocessing.memesuite import read_fimo_output
from antismash.detection.cassis import runners

from .test_cassis import create_fake_record, CassisTestCore


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
            filename_base = os.path.join(self.options.output_dir, "fimo", subdir, "fimo")
            filename = ""
            for extension in ["tsv", "txt"]:  # newer fimo outputs to tsv instead of txt
                test_name = f"{filename_base}.{extension}"
                if os.path.exists(test_name):
                    filename = test_name
                    break
            assert filename, "no compatible files found"

            with open(filename, encoding="utf-8") as handle:
                generated_motifs = read_fimo_output(handle.read())

            # purge the fields old FIMO didn't have
            for motif in generated_motifs:
                motif.alternate_name = ""

            with open(path.get_full_path(__file__, "data", "fake_long_fimo.txt"), encoding="utf-8") as handle:
                expected_motifs = read_fimo_output(handle.read())

            assert len(generated_motifs) == len(expected_motifs)
            for generated, expected in zip(generated_motifs, expected_motifs):
                assert generated == expected
