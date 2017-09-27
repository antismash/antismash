# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from __future__ import print_function, division

import unittest
import os
from collections import defaultdict
from minimock import mock, restore

from antismash.common.secmet import CDSFeature, Record # CDSFeature mocked
from antismash.common.test.helpers import DummyCDS
import antismash.modules.clusterblast.core as core

class TestBlastParsing(unittest.TestCase):
    def setUp(self):
        #this means geneclustergenes inside blastparse will be empty
        mock('CDSFeature.get_accession', returns=None)
        #used by parse_subject, every sequence will be 100 long
        mock('Record.get_cds_mapping', returns=defaultdict(lambda: DummyCDS(1, 101)))
        mock('core.get_cds_lengths', returns={})
        self.sample_data = self.read_sample_data()
        self.sample_data_as_lists = self.file_data_to_lists(self.sample_data)

    def tearDown(self):
        restore()

    def parse_subject_wrapper(self, subject_line):
        geneclustergenes = {}
        seq_record = Record("dummy")
        seqlengths = {}
        return core.parse_subject(subject_line, seqlengths,
                                  geneclustergenes, seq_record)

    def read_sample_data(self, filename="data/diamond_output_sample.txt"):
        path = os.path.join(__file__.rsplit(os.sep, 1)[0], filename)
        return open(path, "r").read()

    def file_data_to_lists(self, data):
        return [line.split("\t") for line in data.rstrip().split("\n")]

    def test_unique_pairings_filter(self):
        data = self.file_data_to_lists(self.sample_data)
        sample = core.uniqueblasthitfilter(data)
        self.assertEqual(len(sample), len(data))
        self.assertEqual(sample, core.uniqueblasthitfilter(data*2))

        # test empty
        data = [[], ["a"], ["abc"]]
        results = core.uniqueblasthitfilter(data)
        self.assertEqual(results, [])

    def verify_subjects_and_clusters_represented(self, subjects,
                                                 cluster_name_to_queries):
        subject_clusters = set()
        for subject in subjects:
            self.assertTrue(subject.genecluster in cluster_name_to_queries)
            subject_clusters.add(subject.genecluster)
        self.assertEqual(sorted(subject_clusters), sorted(cluster_name_to_queries.keys()))

    def test_blastparse(self):
        queries, clusters = core.blastparse(self.sample_data, Record(), 0, 0)

        # check we process the right number of queries
        self.assertEqual(len(queries), len(set([i[0] for i in self.sample_data_as_lists])))

        # check we have entries for every gene_cluster we found
        subjects = [self.parse_subject_wrapper(i) for i in self.sample_data_as_lists]
        self.verify_subjects_and_clusters_represented(subjects, clusters)

        # test perc_coverage threshold (value arbitrary due to mocking)
        coverage_threshold = 650
        queries, clusters = core.blastparse(self.sample_data, Record(), coverage_threshold, 0)
        new_subjects = [s for s in subjects if s.perc_coverage > coverage_threshold]
        assert 0 < len(new_subjects) < len(subjects), "coverage test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test perc_identity threshold
        ident_threshold = 35
        queries, clusters = core.blastparse(self.sample_data, Record(), 0, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold]
        assert 0 < len(new_subjects) < len(subjects), "identity% test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test combo threshold
        queries, clusters = core.blastparse(self.sample_data, Record(), coverage_threshold, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold and s.perc_coverage > coverage_threshold]
        assert 0 < len(new_subjects) < len(subjects), "combo test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

    def test_blastparse_on_empty(self):
        for blast in ["", "\n", "\r\n", "\n\n"]:
            queries, clusters = core.blastparse(blast, Record(), 0, 0)
            self.assertEqual(len(queries), 0)
            self.assertEqual(len(clusters), 0)

    def test_parse_all_single_cluster(self):
        # single cluster to test the thresholds and content
        def parse_all_wrapper(coverage_threshold, ident_threshold):
            clusters_by_number, queries_by_number = core.parse_all_clusters(self.sample_data,
                                Record(), coverage_threshold, ident_threshold)
            # make sure we only found one cluster number
            self.assertEqual(len(clusters_by_number), 1)
            self.assertEqual(list(clusters_by_number.keys()), [24])
            self.assertEqual(len(queries_by_number), 1)
            self.assertEqual(list(queries_by_number.keys()), [24])

            # now test the values of those queries
            queries = queries_by_number[24]
            clusters = clusters_by_number[24]
            return queries, clusters

        queries, clusters = parse_all_wrapper(0, 0)

        # check we process the right number of queries
        self.assertEqual(len(queries), len(set([i[0] for i in self.sample_data_as_lists])))

        # check we have entries for every gene_cluster we found
        subjects = [self.parse_subject_wrapper(i) for i in self.sample_data_as_lists]
        self.verify_subjects_and_clusters_represented(subjects, clusters)

        # test perc_coverage threshold (value arbitrary due to mocking)
        coverage_threshold = 650
        queries, clusters = parse_all_wrapper(coverage_threshold, 0)
        new_subjects = [s for s in subjects if s.perc_coverage > coverage_threshold]
        assert 0 < len(new_subjects) < len(subjects), "coverage test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test perc_identity threshold
        ident_threshold = 35
        queries, clusters = parse_all_wrapper(0, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold]
        assert 0 < len(new_subjects) < len(subjects), "identity% test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test combo threshold
        queries, clusters = parse_all_wrapper(coverage_threshold, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold and s.perc_coverage > coverage_threshold]
        assert 0 < len(new_subjects) < len(subjects), "combo test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

    def test_parse_all_multi_cluster(self):
        # test we partition correctly by cluster number
        sample_data = self.read_sample_data("data/diamond_output_sample_multicluster.txt")
        clusters_by_number, queries_by_number = core.parse_all_clusters(sample_data, Record(), 0, 0)
        self.assertEqual(len(clusters_by_number), 3)
        self.assertEqual(sorted(clusters_by_number.keys()), [1, 2, 4])
        self.assertEqual(len(queries_by_number), 3)
        self.assertEqual(sorted(queries_by_number.keys()), [1, 2, 4])
        for i in [1, 2, 4]:
            self.assertEqual(len(clusters_by_number[i]), i)
            self.assertEqual(len(queries_by_number[i]), i)

    def test_parse_all_empty(self):
        for sample_data in ["", "\n", "\r\n", "\n\n"]:
            clusters, queries = core.parse_all_clusters(sample_data, Record(), 0, 0)
        self.assertEqual(len(clusters), 0)
        self.assertEqual(len(queries), 0)

class TestSubject(unittest.TestCase):
    def test_init(self):
        expected = ["a", "b", 1, 2, "+", "c", 5., 1, 5., 1e-8, "loc"]
        s = core.Subject(*expected)
        got = [s.name, s.genecluster, s.start, s.end, s.strand, s.annotation,
               s.perc_ident, s.blastscore, s.perc_coverage, s.evalue, s.locus_tag]
        for exp, val in zip(expected, got):
            self.assertEqual(exp, val)

class TestQuery(unittest.TestCase):
    def test_init(self):
        query_line = "input|c1|0-759|-|CAG25751.1|putative"
        q = core.Query(query_line, 0)
        self.assertEqual(q.entry, query_line)
        self.assertEqual(q.cluster_number, 1)
        self.assertEqual(q.id, "CAG25751.1")
        self.assertEqual(q.index, 0)
        for container in [q.cluster_name_to_subjects, q.subjects]:
            self.assertEqual(len(container), 0)

    def test_subject_tracking(self):
        q = core.Query("input|c1|0-759|-|CAG25751.1|putative", 0)
        s1 = core.Subject("a1", "a", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        q.add_subject(s1)
        containers = [q.cluster_name_to_subjects, q.subjects]
        #check it was properly added to the various containers
        for container in containers:
            self.assertEqual(len(container), 1)
        self.assertEqual(q.cluster_name_to_subjects["a"], [s1])
        assert list(q.subjects.keys()) == ["a1"]
        s2 = core.Subject("a2", "b", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        q.add_subject(s2)
        for container in containers:
            self.assertEqual(len(container), 2)
        self.assertEqual(q.cluster_name_to_subjects["a"], [s1])
        self.assertEqual(q.cluster_name_to_subjects["b"], [s2])

        # check we don't override when cluster names overlap
        s3 = core.Subject("a3", "a", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        q.add_subject(s3)
        self.assertEqual(len(q.subjects), 3)
        # check the new subject was properly added to the old list
        self.assertEqual(q.cluster_name_to_subjects["a"], [s1, s3])
        self.assertEqual(q.cluster_name_to_subjects["b"], [s2])
        # check ordering preserved on subject names
        self.assertEqual(list(q.subjects.keys()), ["a1", "a2", "a3"])

        # check the getter has the same results as direct access
        self.assertEqual(q.get_subjects_by_cluster("a"), [s1, s3])
        # check that an empty iterable is returned if cluster not known
        self.assertEqual(q.get_subjects_by_cluster("new_name"), [])

# pylint: disable=assigning-non-slot
class TestScore(unittest.TestCase):
    def test_score(self):
        score = core.Score()
        self.assertEqual(score.hits, 0)
        self.assertEqual(score.core_gene_hits, 0)
        self.assertEqual(score.blast_score, 0)
        self.assertEqual(score.synteny_score, 0)
        self.assertEqual(score.core_bonus, 0)
        self.assertEqual(score.sort_score(), 0)
        with self.assertRaises(AttributeError):
            score.other = 0

        score.hits = 1
        self.assertEqual(score.sort_score(), 1)

        score.synteny_score = 1
        self.assertEqual(score.sort_score(), 2)

        score.core_gene_hits = 1
        self.assertAlmostEqual(score.sort_score(), 5.01)
        self.assertEqual(score.core_bonus, 3)

        score.blast_score = 225
        self.assertAlmostEqual(score.sort_score(), 5.01000225)

class TestProtein(unittest.TestCase):
    def test_string_conversion(self):
        protein = core.Protein("n", "l", "a-b", "+", "anno")
        self.assertEqual(str(protein), "l\tn\ta\tb\t+\tanno\n")

        # test name is used when no locus tag
        protein = core.Protein("n", "no_locus_tag", "a-b", "+", "anno")
        self.assertEqual(str(protein), "n\tn\ta\tb\t+\tanno\n")

        # test location
        protein = core.Protein("n", "no_locus_tag", "abb", "+", "anno")
        with self.assertRaises(ValueError) as context:
            str(protein)
        self.assertTrue(str(context.exception).startswith("Invalid location in Protein"))
        protein = core.Protein("n", "no_locus_tag", 123., "+", "anno")
        with self.assertRaises(AttributeError) as context:
            str(protein)

    def test_members(self):
        protein = core.Protein("n", "no_locus_tag", "a-b", "+", "anno")
        protein.locus_tag = "l"
        #if this doesn't raise an exception, __slots__ was removed from Protein
        with self.assertRaises(AttributeError):
            protein.something = 1
# pylint: enable=assigning-non-slot

class TestSubjectParsing(unittest.TestCase):
    def setUp(self):
        self.geneclustergenes = {"CAG25752":""}
        self.seq_record = Record("dummy")
        self.seqlengths = {"CAG25751.1" : 253}
        #used by parse_subject, but only if locus tag not in seqlengths
        mock('core.get_cds_lengths', returns=self.seqlengths)
        mock('Record.get_cds_mapping', returns={"TEST" : DummyCDS(1, 301)})

    def tearDown(self):
        restore()

    def parse_subject_wrapper(self, subject_line):
        return core.parse_subject(subject_line, self.seqlengths,
                                          self.geneclustergenes, self.seq_record)

    def test_all_parsing(self):
        # test known good input
        subject_line = ["input|c1|0-759|-|CAG25751.1|putative",
                        "Y16952|c1|1-759|-|no_locus_tag"
                        "|putative_two-component_system_sensor_kinase|CAG25751",
                        "100.0", "253", "0", "0", "1", "253", "1", "253",
                        "7.2e-129", "465.0"]
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.name, "CAG25751")
        self.assertEqual(sub.blastscore, 465)
        self.assertEqual(sub.evalue, 7.2e-129)
        self.assertEqual(sub.locus_tag, "CAG25751")
        self.assertEqual(sub.genecluster, "Y16952_c1")
        self.assertEqual(sub.start, "1")
        self.assertEqual(sub.end, "759")
        self.assertEqual(sub.strand, "-")
        self.assertEqual(sub.perc_ident, 100)
        self.assertEqual(sub.perc_coverage, 100.)
        self.assertEqual(sub.annotation, "putative_two-component_system_sensor_kinase")

        # with different tag name to test branching
        subject_line[0] = "input|c1|0-759|-|TEST|putative"
        sub = self.parse_subject_wrapper(subject_line)
        # TEST has sequence length of 300 after our mocking
        self.assertAlmostEqual(sub.perc_coverage, (253./300)*100)

        # and check the name change for a name in our dummy geneclustergenes
        subject_line[1] = subject_line[1].replace("CAG25751", "CAG25752")
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.name, "h_CAG25752")

    def test_locus_tag_redirect(self):
        # no_locus_tag normally goes to position 6 instead, but
        subject_line = ["input|c1|0-759|-|CAG25751.1|putative",
                        "Y16952|c1|1-759|-|no_locus_tag"
                        "|putative_two-component_system_sensor_kinase|CAG25751",
                        "100.0", "253", "0", "0", "1", "253", "1", "253",
                        "7.2e-129", "465.0"]
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.name, "CAG25751")

    def test_rounding(self):
        # check the value rounding, this may not be desired, but this test
        # will at least detect change
        subject_line = ["input|c1|0-759|-|CAG25751.1|putative",
                "Y16952|c1|1-759|-|no_locus_tag"
                "|putative_two-component_system_sensor_kinase|CAG25751",
                "99.5", "253", "0", "0", "1", "253", "1", "253",
                "7.2e-129", "464.9"]
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.blastscore, 465)
        self.assertEqual(sub.perc_ident, 100)
        # rounding down instead of up
        subject_line[2] = "99.49"
        subject_line[11] = "465.4"
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.perc_ident, 99)
        self.assertEqual(sub.blastscore, 465)

    def test_short_locus(self):
        # if the
        subject_line = ["input|c1|0-759|-|CAG25751.1|putative",
                "Y16952|c1|1-759|-|some_locus_tag"
                "|putative_two-component_system_sensor_kinase",
                "100.0", "253", "0", "0", "1", "253", "1", "253",
                "7.2e-129", "464.9"]
        sub = self.parse_subject_wrapper(subject_line)
        self.assertEqual(sub.locus_tag, "")

    def test_invalid_input(self):
        original = ["input|c1|0-759|-|CAG25751.1|putative",
                    "Y16952|c1|1-759|-|no_locus_tag"
                    "|putative_two-component_system_sensor_kinase|CAG25751",
                    "100.0", "253", "0", "0", "1", "253", "1", "253",
                    "7.2e-129", "465.0"]
        # test for invalid input
        for i in [2, 3, 11]: #pid, length, blastscore
            subject_line = original[:]
            subject_line[i] = "text"
            with self.assertRaises(ValueError):
                self.parse_subject_wrapper(subject_line)

        # check it breaks if missing info
        for i in range(1, len(original)):
            subject_line = original[:i] + original[i + 1:]
            with self.assertRaises(IndexError):
                self.parse_subject_wrapper(subject_line)

        subject_line = original[:]
        # and if chunks are missing internal parts
        # missing c1 for name
        subject_line[1] = "Y16956|1-759|+|no_locus_tag|putative_two-component_system_sensor_kinase|CAG25751"
        with self.assertRaises(ValueError):
            self.parse_subject_wrapper(subject_line)
        #start missing end
        subject_line[1] = "Y16956|c1|1|-|no_locus_tag|putative_two-component_system_sensor_kinase|CAG25751"
        with self.assertRaises(ValueError):
            self.parse_subject_wrapper(subject_line)
