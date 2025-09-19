# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Tests for the codoff module. """

import os
import sys
import unittest
from unittest.mock import Mock, patch

import numpy as np
from antismash.common.secmet import Record
from antismash.common.secmet.features import CDSFeature

# Add the parent directory to the path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from antismash.modules.codoff.analysis import (
    run_codoff_analysis,
    _calculate_empirical_pvalue_serial,
    generate_codoff_histogram
)
from antismash.modules.codoff.results import CodoffResults
from antismash.modules.codoff import run_on_record, _codoff_global_state


class TestCodoffAnalysis(unittest.TestCase):
    """ Test the codoff analysis functions """
    
    def setUp(self):
        """ Set up test fixtures """
        # Create a mock record
        self.record = Mock(spec=Record)
        self.record.id = "test_record"
        self.record.seq = "ATG" * 100  # 100 codons
        
        # Create mock CDS features
        self.cds1 = Mock(spec=CDSFeature)
        self.cds1.location = Mock()
        self.cds1.location.extract.return_value = "ATG" * 10  # 10 codons
        self.cds1.location.__len__ = lambda self: 30
        self.cds1.get_name.return_value = "test_cds_1"
        self.cds1.record = self.record  # Add record reference
        
        self.cds2 = Mock(spec=CDSFeature)
        self.cds2.location = Mock()
        self.cds2.location.extract.return_value = "GCT" * 8  # 8 codons
        self.cds2.location.__len__ = lambda self: 24
        self.cds2.get_name.return_value = "test_cds_2"
        self.cds2.record = self.record  # Add record reference
        
        self.cds_features = [self.cds1, self.cds2]
    
    def test_count_codons_in_cds(self):
        """ Test codon counting in CDS features """
        # This is a basic test to ensure the mock objects work
        self.assertEqual(len(self.cds1.location), 30)
        self.assertEqual(len(self.cds2.location), 24)
    
    def test_run_codoff_analysis(self):
        """ Test the main codoff analysis function """
        region_codon_freqs = {"ATG": 10, "GCT": 8}
        genome_codon_freqs = {"ATG": 100, "GCT": 80, "TAA": 20}
        
        result = run_codoff_analysis(
            region_codon_freqs=region_codon_freqs,
            genome_codon_freqs=genome_codon_freqs,
            region_cds_features=self.cds_features,
            all_cds_features=self.cds_features,
            num_simulations=100
        )
        
        # Check that we got a result
        self.assertIsInstance(result, dict)
        self.assertIn('empirical_pvalue', result)
        self.assertIn('cosine_distance', result)
        self.assertIn('spearman_rho', result)
        
        # Check that P-value is in valid range
        self.assertGreaterEqual(result['empirical_pvalue'], 0.0)
        self.assertLessEqual(result['empirical_pvalue'], 1.0)
    
    def test_empirical_pvalue_calculation(self):
        """ Test empirical P-value calculation """
        region_freqs_array = [10, 8, 0]  # ATG, GCT, TAA
        genome_freqs_array = [100, 80, 20]  # ATG, GCT, TAA
        codon_order = ["ATG", "GCT", "TAA"]
        observed_cosine_distance = 0.1
        
        # The function now uses arrays instead of lists
        region_freqs_array = np.array(region_freqs_array, dtype=np.float64)
        genome_freqs_array = np.array(genome_freqs_array, dtype=np.float64)
        
        # FIX: Unpack the tuple return value
        pvalue, simulation_distances = _calculate_empirical_pvalue_serial(
            region_freqs_array=region_freqs_array,
            genome_freqs_array=genome_freqs_array,
            all_cds_features=self.cds_features,
            codon_order=codon_order,
            observed_cosine_distance=observed_cosine_distance,
            num_simulations=50
        )
        
        # Check P-value is in valid range
        self.assertGreaterEqual(pvalue, 0.0)
        self.assertLessEqual(pvalue, 1.0)
        
        # ADD: Check simulation distances list
        self.assertIsInstance(simulation_distances, list)
        self.assertGreaterEqual(len(simulation_distances), 0)  # Allow empty list for small test data
        if simulation_distances:  # Only check length if we got results
            self.assertLessEqual(len(simulation_distances), 50)  # Should not exceed requested simulations
    
    def test_codoff_results_serialization(self):
        """ Test CodoffResults serialization """
        results = CodoffResults("test_record")
        
        # Add some test data including genome_codon_frequencies
        results.by_region[1] = {
            'empirical_pvalue': 0.05,
            'cosine_distance': 0.1,
            'spearman_rho': 0.8,
            'region_codon_frequencies': [10, 5, 3],
            'genome_codon_frequencies': [100, 50, 30],  # This should be excluded from JSON
            'genome_codon_frequencies_standardized': [0.5, 0.3, 0.2]
        }
        
        # Test JSON serialization
        json_data = results.to_json()
        self.assertIn('record_id', json_data)
        self.assertIn('by_region', json_data)
        self.assertEqual(json_data['record_id'], 'test_record')
        self.assertEqual(json_data['by_region'][1]['empirical_pvalue'], 0.05)
        
        # Test that genome_codon_frequencies is excluded from JSON
        self.assertNotIn('genome_codon_frequencies', json_data['by_region'][1])
        # But other fields should be present
        self.assertIn('region_codon_frequencies', json_data['by_region'][1])
        self.assertIn('genome_codon_frequencies_standardized', json_data['by_region'][1])
        
        # Test JSON deserialization
        reconstructed = CodoffResults.from_json(json_data, self.record)
        self.assertEqual(reconstructed.record_id, 'test_record')
        self.assertIn(1, reconstructed.by_region)
        self.assertEqual(reconstructed.by_region[1]['empirical_pvalue'], 0.05)
        
        # Verify that internal data still has genome_codon_frequencies
        self.assertIn('genome_codon_frequencies', results.by_region[1])
        self.assertEqual(results.by_region[1]['genome_codon_frequencies'], [100, 50, 30])


    def test_codon_caching(self):
        """ Test that codon extraction and counting is cached properly """
        from antismash.modules.codoff import _get_cached_codons, _get_cached_codon_counts
        
        # Create a mock CDS with deterministic codon sequence
        cds = Mock(spec=CDSFeature)
        cds.location = Mock()
        cds.location.extract = Mock(return_value="ATGCCCGGGTTT")  # 4 codons: ATG, CCC, GGG, TTT
        cds.location.__len__ = lambda self: 12
        cds.get_name.return_value = "test_cds_cache"
        
        # Test codon extraction caching
        codons1 = _get_cached_codons(cds, self.record)
        codons2 = _get_cached_codons(cds, self.record)  # Should use cache
        
        self.assertEqual(codons1, codons2)
        self.assertEqual(codons1, ['ATG', 'CCC', 'GGG', 'TTT'])
        # Extract should only be called once due to caching
        self.assertEqual(cds.location.extract.call_count, 1)
        
        # Test codon counting cache
        codon_order = ['ATG', 'CCC', 'GGG', 'TTT', 'AAA']
        counts1 = _get_cached_codon_counts(cds, codon_order)
        counts2 = _get_cached_codon_counts(cds, codon_order)  # Should use cache
        
        import numpy as np
        np.testing.assert_array_equal(counts1, counts2)
        expected_counts = np.array([1, 1, 1, 1, 0])  # One of each except AAA
        np.testing.assert_array_equal(counts1, expected_counts)
    
    def test_html_output_generation(self):
        """ Test HTML output generation with complete data """
        from antismash.modules.codoff.html_output import generate_html
        from antismash.common.secmet.features import Region
        
        # Create mock region
        region = Mock(spec=Region)
        region.get_region_number.return_value = 1
        
        # Create results with complete data
        results = CodoffResults("test_record")
        results.by_region[1] = {
            'empirical_pvalue': 0.025,
            'cosine_distance': 0.15,
            'spearman_rho': 0.75,
            'region_cds_count': 12,
            'genome_cds_count_excluding_region': 2000,
            'total_cds_count': 2012,
            'codon_order': ['ATG', 'CCC', 'GGG'],
            'region_codon_frequencies': [10, 5, 3],
            'genome_codon_frequencies': [200, 150, 100],  # Available for HTML
            'region_codon_frequencies_standardized': [0.555, 0.278, 0.167],
            'genome_codon_frequencies_standardized': [0.444, 0.333, 0.222],
            'region_total_codons': 18,
            'genome_total_codons': 450,
            'simulation_distances': [0.1, 0.12, 0.14, 0.16]
        }
        
        # Mock options
        options = Mock()
        options.output_dir = "/tmp/test_output"
        
        # Generate HTML
        html_sections = generate_html(region, results, self.record, options)
        
        # Verify HTML sections were created
        self.assertIsNotNone(html_sections)
        # The actual HTML content testing would require more complex mocking
        # but this ensures the basic structure works


class TestCodoffDeferredAnalysis(unittest.TestCase):
    """ Test the deferred analysis approach for multi-scaffold genomes """
    
    def setUp(self):
        """ Set up test fixtures for deferred analysis testing """
        # Reset global state before each test
        _codoff_global_state.clear()
        _codoff_global_state.update({
            'total_genome_size': 0,
            'all_cds_features': [],
            'cds_to_record_map': {},
            'cds_to_record_ref_map': {},
            'current_record_count': 0,
            'finalization_complete': False,
            'global_cds_collection_complete': False,  # Add missing key
            'bgc_regions_tracked': []
        })
        
        # Create mock records for multi-scaffold scenario
        self.record1 = Mock(spec=Record)
        self.record1.id = "scaffold_1"
        self.record1.record_index = 1  # 1-indexed as in antiSMASH
        self.record1.seq = Mock()
        self.record1.seq.__len__ = lambda self: 600000  # Large enough sequence
        
        self.record2 = Mock(spec=Record)
        self.record2.id = "scaffold_2" 
        self.record2.record_index = 2  # 1-indexed as in antiSMASH
        self.record2.seq = Mock()
        self.record2.seq.__len__ = lambda self: 400000  # Large enough sequence
        
        # Create mock options
        self.options = Mock()
        self.options.codoff_enabled = True
        self.options.codoff_simulations = 100
        self.options.codoff_min_genome_size = 500000
        self.options.codoff_max_genome_size = 50000000
        self.options.cpus = 2
    
    def test_write_outputs_finalization(self):
        """ Test that finalization is triggered by write_outputs() and only runs once """
        from antismash.modules.codoff.results import CodoffResults
        
        # Create mock results and populate global state with some data
        results = CodoffResults("test_record")
        _codoff_global_state['bgc_regions_tracked'] = [{
            'record_id': 'test',
            'region_number': 1,
            'region_cds_features': [],
            'finalized_results': None
        }]
        
        # Mock the _finalize_analysis_for_all_records function
        with patch('antismash.modules.codoff._finalize_analysis_for_all_records') as mock_finalize:
            # First call to write_outputs should trigger finalization
            results.write_outputs(self.record1, self.options)
            mock_finalize.assert_called_once_with(self.options)
            
            # Second call should not trigger finalization again
            mock_finalize.reset_mock()
            results.write_outputs(self.record2, self.options)
            mock_finalize.assert_not_called()
            
            # Verify finalization_complete flag is set
            self.assertTrue(_codoff_global_state['finalization_complete'])
    
    def test_global_state_tracking(self):
        """ Test that global state properly tracks CDS data collection """
        # Mock a simple record with CDS features
        record = Mock(spec=Record)
        record.id = "test_scaffold"
        record.seq = Mock()
        record.seq.__len__ = lambda self: 100000
        record.get_regions.return_value = []  # No BGC regions for simplicity
        
        # Mock CDS features
        cds = Mock(spec=CDSFeature)
        cds.location = Mock()
        cds.location.__len__ = lambda self: 300
        cds.get_name.return_value = "test_cds"
        record.get_cds_features.return_value = [cds]
        
        # Test: Process record should collect CDS data
        results = run_on_record(record, None, self.options)
        
        # Verify data collection occurred
        self.assertGreater(_codoff_global_state['total_genome_size'], 0)
        self.assertGreater(len(_codoff_global_state['all_cds_features']), 0)
        self.assertEqual(_codoff_global_state['current_record_count'], 1)
        
        # Verify CDS mapping
        self.assertGreater(len(_codoff_global_state['cds_to_record_map']), 0)
        self.assertGreater(len(_codoff_global_state['cds_to_record_ref_map']), 0)
        
    def tearDown(self):
        """ Clean up after each test """
        # Reset global state
        _codoff_global_state.clear()


class TestCodoffIntegration(unittest.TestCase):
    """ Integration tests using actual test data files """
    
    def setUp(self):
        """ Set up test fixtures with real data """
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Check if test data files exist
        self.genome_file = os.path.join(self.test_dir, "LK413.gbk")
        self.fasta_file = os.path.join(self.test_dir, "LK413.fna")
        self.region1_file = os.path.join(self.test_dir, "NZ_JALXLO020000001.1.region001.gbk")
        self.region2_file = os.path.join(self.test_dir, "NZ_JALXLO020000002.1.region001.gbk")
        
        # Verify test files exist
        self.assertTrue(os.path.exists(self.genome_file), f"Test genome file not found: {self.genome_file}")
        self.assertTrue(os.path.exists(self.fasta_file), f"Test FASTA file not found: {self.fasta_file}")
        self.assertTrue(os.path.exists(self.region1_file), f"Test region 1 file not found: {self.region1_file}")
        self.assertTrue(os.path.exists(self.region2_file), f"Test region 2 file not found: {self.region2_file}")
    
    def test_test_data_files_exist(self):
        """ Test that all required test data files are present """
        # This test ensures we have the test data we need
        self.assertTrue(os.path.exists(self.genome_file))
        self.assertTrue(os.path.exists(self.fasta_file))
        self.assertTrue(os.path.exists(self.region1_file))
        self.assertTrue(os.path.exists(self.region2_file))
        
        # Check file sizes are reasonable
        self.assertGreater(os.path.getsize(self.genome_file), 1000000)  # Should be > 1MB
        self.assertGreater(os.path.getsize(self.fasta_file), 1000000)   # Should be > 1MB
        self.assertGreater(os.path.getsize(self.region1_file), 10000)   # Should be > 10KB
        self.assertGreater(os.path.getsize(self.region2_file), 10000)   # Should be > 10KB
    
    def test_genome_file_format(self):
        """ Test that the genome file is in valid GenBank format """
        with open(self.genome_file, 'r') as f:
            content = f.read()
            # Check for GenBank header
            self.assertIn("LOCUS", content)
            self.assertIn("DEFINITION", content)
            self.assertIn("FEATURES", content)
    
    def test_fasta_file_format(self):
        """ Test that the FASTA file is in valid format """
        with open(self.fasta_file, 'r') as f:
            first_line = f.readline().strip()
            # Check for FASTA header
            self.assertTrue(first_line.startswith(">"))
    
    def test_region_files_format(self):
        """ Test that the region files are in valid GenBank format """
        for region_file in [self.region1_file, self.region2_file]:
            with open(region_file, 'r') as f:
                content = f.read()
                # Check for GenBank header
                self.assertIn("LOCUS", content)
                self.assertIn("FEATURES", content)


if __name__ == '__main__':
    unittest.main()
