# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

from collections import OrderedDict
import os
import unittest
from unittest import mock as unittest_mock

from minimock import mock, restore
from helperlibs.wrappers.io import TemporaryDirectory

from antismash import config
from antismash.common.secmet import Record, Region
from antismash.common.test.helpers import DummyCDS, DummyCluster, DummySuperCluster
import antismash.modules.clusterblast.core as core


class TestBlastParsing(unittest.TestCase):
    def setUp(self):
        # used by parse_subject, every sequence will be 100 long
        mock('Record.get_cds_by_name', returns=DummyCDS(1, 101))
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
        sample = core.remove_duplicate_hits(data)
        self.assertEqual(len(sample), len(data))
        self.assertEqual(sample, core.remove_duplicate_hits(data*2))

        # test empty
        data = [[], ["a"], ["abc"]]
        results = core.remove_duplicate_hits(data)
        self.assertEqual(results, [])

    def verify_subjects_and_clusters_represented(self, subjects,
                                                 cluster_name_to_queries):
        subject_clusters = set()
        for subject in subjects:
            self.assertTrue(subject.genecluster in cluster_name_to_queries)
            subject_clusters.add(subject.genecluster)
        self.assertEqual(sorted(subject_clusters), sorted(cluster_name_to_queries))

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
        assert new_subjects and len(new_subjects) < len(subjects), "coverage test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test perc_identity threshold
        ident_threshold = 35
        queries, clusters = core.blastparse(self.sample_data, Record(), 0, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold]
        assert new_subjects and len(new_subjects) < len(subjects), "identity% test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test combo threshold
        queries, clusters = core.blastparse(self.sample_data, Record(), coverage_threshold, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold and s.perc_coverage > coverage_threshold]
        assert new_subjects and len(new_subjects) < len(subjects), "combo test has become meaningless"
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
            self.assertEqual(list(clusters_by_number), [24])
            self.assertEqual(len(queries_by_number), 1)
            self.assertEqual(list(queries_by_number), [24])

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
        assert new_subjects and len(new_subjects) < len(subjects), "coverage test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test perc_identity threshold
        ident_threshold = 35
        queries, clusters = parse_all_wrapper(0, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold]
        assert new_subjects and len(new_subjects) < len(subjects), "identity% test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

        # test combo threshold
        queries, clusters = parse_all_wrapper(coverage_threshold, ident_threshold)
        new_subjects = [s for s in subjects if s.perc_ident > ident_threshold and s.perc_coverage > coverage_threshold]
        assert new_subjects and len(new_subjects) < len(subjects), "combo test has become meaningless"
        self.verify_subjects_and_clusters_represented(new_subjects, clusters)

    def test_parse_all_multi_cluster(self):
        # test we partition correctly by cluster number
        sample_data = self.read_sample_data("data/diamond_output_sample_multicluster.txt")
        clusters_by_number, queries_by_number = core.parse_all_clusters(sample_data, Record(), 0, 0)
        self.assertEqual(len(clusters_by_number), 3)
        self.assertEqual(sorted(clusters_by_number), [1, 2, 4])
        self.assertEqual(len(queries_by_number), 3)
        self.assertEqual(sorted(queries_by_number), [1, 2, 4])
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
        s = core.Subject(*expected)  # pylint: disable=invalid-name
        got = [s.name, s.genecluster, s.start, s.end, s.strand, s.annotation,
               s.perc_ident, s.blastscore, s.perc_coverage, s.evalue, s.locus_tag]
        for exp, val in zip(expected, got):
            self.assertEqual(exp, val)


class TestQuery(unittest.TestCase):
    def test_init(self):
        query_line = "input|c1|0-759|-|CAG25751.1|putative"
        query = core.Query(query_line, 0)
        self.assertEqual(query.entry, query_line)
        self.assertEqual(query.cluster_number, 1)
        self.assertEqual(query.id, "CAG25751.1")
        self.assertEqual(query.index, 0)
        for container in [query.cluster_name_to_subjects, query.subjects]:
            self.assertEqual(len(container), 0)

    def test_subject_tracking(self):
        query = core.Query("input|c1|0-759|-|CAG25751.1|putative", 0)
        sub1 = core.Subject("a1", "a", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        query.add_subject(sub1)
        containers = [query.cluster_name_to_subjects, query.subjects]
        # check it was properly added to the various containers
        for container in containers:
            self.assertEqual(len(container), 1)
        self.assertEqual(query.cluster_name_to_subjects["a"], [sub1])
        assert list(query.subjects) == ["a1"]
        sub2 = core.Subject("a2", "b", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        query.add_subject(sub2)
        for container in containers:
            self.assertEqual(len(container), 2)
        self.assertEqual(query.cluster_name_to_subjects["a"], [sub1])
        self.assertEqual(query.cluster_name_to_subjects["b"], [sub2])

        # check we don't override when cluster names overlap
        sub3 = core.Subject("a3", "a", 1, 2, "+", "c", 0.5, 1, 0.5, 1e-8, "loc")
        query.add_subject(sub3)
        self.assertEqual(len(query.subjects), 3)
        # check the new subject was properly added to the old list
        self.assertEqual(query.cluster_name_to_subjects["a"], [sub1, sub3])
        self.assertEqual(query.cluster_name_to_subjects["b"], [sub2])
        # check ordering preserved on subject names
        self.assertEqual(list(query.subjects), ["a1", "a2", "a3"])

        # check the getter has the same results as direct access
        self.assertEqual(query.get_subjects_by_cluster("a"), [sub1, sub3])
        # check that an empty iterable is returned if cluster not known
        self.assertEqual(query.get_subjects_by_cluster("new_name"), [])


# pylint: disable=assigning-non-slot
class TestScore(unittest.TestCase):
    def test_score(self):
        score = core.Score()
        self.assertEqual(score.hits, 0)
        self.assertEqual(score.core_gene_hits, 0)
        self.assertEqual(score.blast_score, 0)
        self.assertEqual(score.synteny_score, 0)
        self.assertEqual(score.core_bonus, 0)
        self.assertEqual(score.sort_score(), (0, 0))
        with self.assertRaises(AttributeError):
            score.other = 0

        score.hits = 1
        assert score.score == 1
        assert score.sort_score() == (1, 0)

        score.synteny_score = 1
        assert score.score == 2
        assert score.sort_score() == (2, 0)

        score.core_gene_hits = 1
        assert score.score == 6
        assert score.sort_score() == (6, 0)

        score.blast_score = 225
        assert score.score == 6
        assert score.sort_score() == (6, 225)

        score.blast_score = 1225
        assert score.score == 6
        assert score.sort_score() == (6, 1225)


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
        # if this doesn't raise an exception, __slots__ was removed from Protein
        with self.assertRaises(AttributeError):
            protein.something = 1
# pylint: enable=assigning-non-slot


class TestSubjectParsing(unittest.TestCase):
    def setUp(self):
        self.geneclustergenes = {"CAG25752": ""}
        self.seq_record = Record("dummy")
        self.seqlengths = {"CAG25751.1": 253}
        # used by parse_subject, but only if locus tag not in seqlengths
        mock('core.get_cds_lengths', returns=self.seqlengths)
        mock('Record.get_cds_by_name', returns=DummyCDS(1, 301))

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
        self.assertEqual(sub.start, 1)
        self.assertEqual(sub.end, 759)
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
        for i in [2, 3, 11]:  # pid, length, blastscore
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
        # start missing end
        subject_line[1] = "Y16956|c1|1|-|no_locus_tag|putative_two-component_system_sensor_kinase|CAG25751"
        with self.assertRaises(ValueError):
            self.parse_subject_wrapper(subject_line)


class TestInputGeneration(unittest.TestCase):
    def setUp(self):
        self.index = 0
        self.old_blast_inputs = core.create_blast_inputs
        core.create_blast_inputs = self.dummy_blast_inputs
        self.dummy_cluster = DummyCluster(1, 100)
        self.supercluster = DummySuperCluster([self.dummy_cluster])
        self.region = Region([self.supercluster], [])
        self.regions = [self.region, self.region]

    def tearDown(self):
        core.create_blast_inputs = self.old_blast_inputs

    def dummy_blast_inputs(self, cluster):
        names = []
        seqs = []
        for _ in cluster.cds_children:
            index = self.index
            self.index += 1
            names.append("L%d" % index)
            seqs.append("S%d" % index)
        return names, seqs

    def add_cdses_to_region(self, cdses):
        for cds in cdses:
            self.region.add_cds(cds)

    def test_empty(self):
        with TemporaryDirectory(change=True):
            with self.assertRaisesRegex(ValueError, "Diamond search space contains no sequences"):
                core.write_fastas_with_all_genes(self.regions, "test")

    def test_bad_partitions(self):
        with TemporaryDirectory(change=True):
            for i in [-10, -1, 0]:
                with self.assertRaisesRegex(ValueError, "Partitions must be greater than 0"):
                    core.write_fastas_with_all_genes(self.regions, "test", partitions=i)
            for i in ["str", None, 1.5]:
                with self.assertRaisesRegex(TypeError, "Partitions must be an int greater than 0"):
                    core.write_fastas_with_all_genes(self.regions, "test", partitions=i)

    def test_single_file(self):
        self.add_cdses_to_region([DummyCDS(1, i) for i in range(3, 6)])
        with TemporaryDirectory(change=True):
            files = core.write_fastas_with_all_genes(self.regions, "test.fasta")
            assert files == ["test.fasta"]
            assert os.path.exists("test.fasta")
            expected = "".join(">L{0}\nS{0}\n".format(i) for i in range(len(self.regions)*3))
            assert open("test.fasta").read() == expected

    def test_single_partition(self):
        self.add_cdses_to_region([DummyCDS(1, i) for i in range(3, 6)])
        with TemporaryDirectory(change=True):
            files = core.write_fastas_with_all_genes(self.regions, "test.fasta", partitions=1)
            assert files == ["test.fasta"]
            assert os.path.exists("test.fasta")
            expected = "".join(">L{0}\nS{0}\n".format(i) for i in range(len(self.regions)*3))
            assert open("test.fasta").read() == expected

    def test_multiple_files(self):
        self.add_cdses_to_region([DummyCDS(1, i) for i in range(3, 6)])
        for partitions in [2, 3]:
            with TemporaryDirectory(change=True):
                self.index = 0
                chunk_size = (len(self.regions) * 3) // partitions
                files = core.write_fastas_with_all_genes(self.regions, "test.fasta", partitions=partitions)
                assert files == ["test%d.fasta" % i for i in range(partitions)]
                for index in range(partitions):
                    assert os.path.exists("test%d.fasta" % index)
                    print(index, chunk_size)
                    contents = open("test%d.fasta" % index).read()
                    assert contents.count(">") == chunk_size
                    expected = "".join(">L{0}\nS{0}\n".format(i + index * chunk_size) for i in range(chunk_size))
                    assert contents == expected


class TestOrthologousGroups(unittest.TestCase):
    def setUp(self):
        self.query_lines = ["input|c1|0-759|-|CAG25751.1|putative",
                            "input|c1|0-759|-|CAG25751.2|putative",
                            "input|c1|0-759|-|CAG25751.3|putative"]
        self.queries = OrderedDict()
        for line in self.query_lines:
            self.queries[line] = core.Query(line, 0)
            assert not self.queries[line].subjects

    def run_base_comparison(self, clusters):
        groups = core.find_internal_orthologous_groups(self.queries, clusters)
        assert groups and len(groups) <= len(clusters)
        return groups

    def set_query_subjects(self, index, subjects):
        self.queries[self.query_lines[index]].subjects = subjects

    def set_independent_subjects(self):
        self.set_query_subjects(0, ["a1", "a2", "a3"])
        self.set_query_subjects(1, ["b1", "b2", "b3"])
        self.set_query_subjects(2, ["c1", "c2", "c3"])

    def test_no_subjects_and_single_not_present(self):
        groups = self.run_base_comparison(["other|c1|0-759|-|other_1|putative"])
        assert groups == [["other_1"]]

    def test_no_subjects_and_multiple_not_present(self):
        groups = self.run_base_comparison(["other|c1|0-759|-|other_1|putative",
                                           "other|c1|0-759|-|other_2|putative"])
        assert groups == [["other_1"], ["other_2"]]

    def test_no_subjects_and_clusters_present(self):
        groups = self.run_base_comparison(self.query_lines[::2])
        assert groups == [["CAG25751.1"], ["CAG25751.3"]]

    def test_independent_subjects_and_not_present(self):
        self.set_independent_subjects()
        clusters = ["other|c1|0-759|-|other_1|putative", "other|c1|0-759|-|other_2|putative"]
        groups = self.run_base_comparison(clusters)
        assert groups == [["other_1"], ["other_2"]]

    def test_independent_subjects_and_clusters_present(self):  # pylint: disable=invalid-name
        self.set_independent_subjects()
        groups = self.run_base_comparison(self.query_lines)
        assert groups == [['CAG25751.1', 'a1', 'a2', 'a3'],
                          ['CAG25751.2', 'b1', 'b2', 'b3'],
                          ['CAG25751.3', 'c1', 'c2', 'c3']]

    def test_partial_overlapping_subjects_and_not_present(self):  # pylint: disable=invalid-name
        self.set_independent_subjects()
        self.set_query_subjects(1, ["b1", "a2", "b3"])
        clusters = ["other|c1|0-759|-|other_1|putative", "other|c1|0-759|-|other_2|putative"]
        groups = self.run_base_comparison(clusters)
        assert groups == [["other_1"], ["other_2"]]

    def test_partial_overlapping_subjects_and_present(self):  # pylint: disable=invalid-name
        self.set_independent_subjects()
        self.set_query_subjects(1, ["b1", "a2", "b3"])
        groups = self.run_base_comparison(self.query_lines)
        assert groups == [['CAG25751.1', 'CAG25751.2', 'a1', 'a2', 'a3', 'b1', 'b3'],
                          ['CAG25751.3', 'c1', 'c2', 'c3']]

    def test_multiple_overlapping_and_not_present(self):
        self.set_independent_subjects()
        self.set_query_subjects(1, ["b1", "a2", "b3"])
        clusters = ["other|c1|0-759|-|other_1|putative", "other|c1|0-759|-|other_2|putative"]
        groups = self.run_base_comparison(clusters)
        assert groups == [["other_1"], ["other_2"]]

    def test_multiple_overlapping_and_present(self):
        self.set_independent_subjects()
        self.set_query_subjects(1, ["b1", "a2", "c3"])
        groups = self.run_base_comparison(self.query_lines)
        assert groups == [['CAG25751.1', 'CAG25751.2', 'CAG25751.3', 'a1', 'a2', 'a3', 'b1', 'c1', 'c2', 'c3']]

    def test_with_uniqueness_modifier(self):
        # tests that the 'h_' prefix used to differentiate names is removed
        # if this test breaks because a better naming system is implemented, that's fine
        self.set_independent_subjects()
        self.set_query_subjects(2, ["_c1", "h_c2", "c_3"])
        groups = self.run_base_comparison(self.query_lines)
        assert groups == [['CAG25751.1', 'a1', 'a2', 'a3'],
                          ['CAG25751.2', 'b1', 'b2', 'b3'],
                          ['CAG25751.3', '_c1', 'c2', 'c_3']]


class TestReferenceProteinLoading(unittest.TestCase):
    def mock_with(self, content):
        # unittest_mock.mock_open doesn't handle the __iter__ case, so work around it
        mocked_open = unittest_mock.mock_open(read_data=content)
        mocked_open.return_value.__iter__ = lambda self: iter(self.readline, '')
        return mocked_open

    def setUp(self):
        config.build_config({})

    def tearDown(self):
        config.destroy_config()

    def test_standard(self):
        hit = ">CVNH01000008|c1|65549-69166|-|BN1184_AH_00620|Urea_carboxylase_{ECO:0000313}|CRH36422"
        with unittest_mock.patch("builtins.open", self.mock_with(hit)):
            proteins = core.load_reference_proteins(set(), "clusterblast")
        assert len(proteins) == 1
        protein = proteins["CRH36422"]
        assert protein.name == "CRH36422"
        assert protein.annotations == "Urea_carboxylase_{ECO:0000313}"

    def test_non_standard(self):
        hit = ">CVNH01000008|c1|65549-69166|-|BN1184_AH_00620|Urea_carboxylase_{ECO:0000313|EMBL:CCF11062.1}|CRH36422"
        with unittest_mock.patch("builtins.open", self.mock_with(hit)):
            proteins = core.load_reference_proteins(set(), "clusterblast")
        assert len(proteins) == 1
        protein = proteins["CRH36422"]
        assert protein.name == "CRH36422"
        assert protein.annotations == "Urea_carboxylase_{ECO:0000313|EMBL:CCF11062.1}"


class TestMissingProteinCleanup(unittest.TestCase):
    def test_missing_removed(self):
        proteins = {}
        for name in "ABDEFH":
            proteins[name] = core.Protein(name, "dummylocus", "1-5", "+", "annotation")
        clusters = {
            "1": core.ReferenceCluster("1", "c1", ["A", "B", "C"], "desc", "type", ["tag1"]),
            "2": core.ReferenceCluster("2", "c1", ["D", "E"], "desc", "type", ["tag3"]),
            "3": core.ReferenceCluster("3", "c1", ["F", "G"], "desc", "type", ["tag2"]),
            "4": core.ReferenceCluster("4", "c2", ["H"], "desc", "type", ["tag4"])
        }
        core.strip_clusters_missing_proteins(clusters, proteins)
        assert sorted(clusters) == ["2", "4"]
        assert sorted(proteins) == ["D", "E", "H"]
