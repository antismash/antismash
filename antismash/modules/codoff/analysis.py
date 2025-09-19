# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Core analysis functions for the codoff module. """

import logging
from typing import Any, Dict, List, Tuple
import numpy as np
from scipy import stats, spatial
from antismash.common.secmet.features import CDSFeature


def _prepare_codon_arrays(region_codon_freqs: Dict[str, int], genome_codon_freqs: Dict[str, int]) -> Tuple[List[str], np.ndarray, np.ndarray]:
    """ Prepare codon arrays for analysis """
    all_codons = set(region_codon_freqs.keys()) | set(genome_codon_freqs.keys())
    codon_order = sorted(all_codons)
    
    region_freqs_array = np.array([region_codon_freqs.get(codon, 0) for codon in codon_order], dtype=np.float64)
    genome_freqs_array = np.array([genome_codon_freqs.get(codon, 0) for codon in codon_order], dtype=np.float64)
    
    return codon_order, region_freqs_array, genome_freqs_array


def _check_analysis_viability(region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray) -> bool:
    """ Check if we have enough data for analysis """
    return not (np.sum(region_freqs_array) == 0 and np.sum(genome_freqs_array) == 0)


def _create_default_results(region_cds_features: List[CDSFeature], all_cds_features: List[CDSFeature], 
                           codon_order: List[str], region_freqs_array: np.ndarray, 
                           genome_freqs_array: np.ndarray, total_cds_features: int = 0) -> Dict[str, Any]:
    """ Create default results when analysis is not viable """
    return {
        'empirical_pvalue': 1.0,
        'cosine_distance': 1.0,
        'spearman_rho': 0.0,
        'region_cds_count': len(region_cds_features),
        'genome_cds_count_excluding_region': len(all_cds_features),  # CDS count excluding the current region
        'total_cds_count': total_cds_features,  # Total CDS count across all scaffolds
        'codon_order': codon_order,
        'region_codon_frequencies': region_freqs_array.tolist(),
        'genome_codon_frequencies': genome_freqs_array.tolist(),  # Keep for HTML, exclude from JSON
        'region_codon_frequencies_standardized': [0.0] * len(codon_order),
        'genome_codon_frequencies_standardized': [0.0] * len(codon_order),
        'region_total_codons': 0,
        'genome_total_codons': 0
    }


def _calculate_cosine_distance(region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray) -> float:
    """ Calculate cosine distance between region and genome codon frequencies """
    try:
        if np.std(region_freqs_array) == 0 and np.std(genome_freqs_array) == 0:
            return 0.0  # Both are constant, so they're identical
        
        # Add small epsilon to avoid division by zero
        eps = 1e-10
        cosine_distance = spatial.distance.cosine(region_freqs_array + eps, genome_freqs_array + eps)
        
        # Handle NaN/Inf values
        if np.isnan(cosine_distance) or np.isinf(cosine_distance):
            return 1.0
        return cosine_distance
    except (ValueError, RuntimeWarning):
        return 1.0


def _calculate_spearman_correlation(region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray) -> float:
    """ Calculate Spearman correlation between region and genome codon frequencies """
    try:
        if np.std(region_freqs_array) == 0 or np.std(genome_freqs_array) == 0:
            return 0.0  # Can't correlate with constant data
        
        spearman_rho, _ = stats.spearmanr(region_freqs_array, genome_freqs_array)
        # Handle NaN values
        if np.isnan(spearman_rho):
            return 0.0
        return spearman_rho
    except (ValueError, RuntimeWarning):
        return 0.0


def _calculate_standardized_frequencies(region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """ Calculate standardized (proportional) frequencies """
    region_total = np.sum(region_freqs_array)
    genome_total = np.sum(genome_freqs_array)
    
    if region_total > 0:
        region_standardized = region_freqs_array / region_total
    else:
        region_standardized = np.zeros_like(region_freqs_array)
    
    if genome_total > 0:
        genome_standardized = genome_freqs_array / genome_total
    else:
        genome_standardized = np.zeros_like(genome_freqs_array)
    
    return region_standardized, genome_standardized


def _prepare_final_results(emp_pval: float, cosine_distance: float, spearman_rho: float,
                          region_cds_features: List[CDSFeature], all_cds_features: List[CDSFeature],
                          codon_order: List[str], region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray,
                          region_standardized: np.ndarray, genome_standardized: np.ndarray,
                          simulation_distances: List[float], total_cds_features: int = 0) -> Dict[str, Any]:
    """ Prepare the final results dictionary.
    
        Returns:
            dictionary containing all analysis results
    """
    region_total = np.sum(region_freqs_array)
    genome_total = np.sum(genome_freqs_array)
    
    return {
        'empirical_pvalue': float(emp_pval),
        'cosine_distance': float(cosine_distance),
        'spearman_rho': float(spearman_rho),
        'region_cds_count': len(region_cds_features),
        'genome_cds_count_excluding_region': len(all_cds_features),  # CDS count excluding the current region
        'total_cds_count': total_cds_features,  # Total CDS count across all scaffolds
        'codon_order': codon_order,
        'region_codon_frequencies': [int(x) for x in region_freqs_array],
        'genome_codon_frequencies': [int(x) for x in genome_freqs_array],  # Keep for HTML, exclude from JSON
        'region_codon_frequencies_standardized': [float(x) for x in region_standardized],
        'genome_codon_frequencies_standardized': [float(x) for x in genome_standardized],
        'region_total_codons': int(region_total),
        'genome_total_codons': int(genome_total),
        'simulation_distances': [float(x) for x in simulation_distances] if simulation_distances else []
    }


def run_codoff_analysis(region_codon_freqs: Dict[str, int], 
                        genome_codon_freqs: Dict[str, int], region_cds_features: List[CDSFeature],
                        all_cds_features: List[CDSFeature],
                        num_simulations: int = 10000, total_cds_count: int = 0) -> Dict[str, Any]:
    """ Run the complete codoff analysis for a single region.
    
        Arguments:
            region_codon_freqs: codon frequencies for the region
            genome_codon_freqs: codon frequencies for the genome (excluding region)
            region_cds_features: list of CDS features in the region
            all_cds_features: list of all CDS features for genome-wide analysis
            num_simulations: number of Monte Carlo simulations
            total_cds_count: total CDS count across all scaffolds
            
        Returns:
            dictionary containing analysis results
    """
    codon_order, region_freqs_array, genome_freqs_array = _prepare_codon_arrays(region_codon_freqs, genome_codon_freqs)
    
    if not _check_analysis_viability(region_freqs_array, genome_freqs_array):
        logging.warning("Insufficient data for codoff analysis")
        return _create_default_results(region_cds_features, all_cds_features, codon_order,
                                     region_freqs_array, genome_freqs_array, len(all_cds_features))
    
    cosine_distance = _calculate_cosine_distance(region_freqs_array, genome_freqs_array)
    spearman_rho = _calculate_spearman_correlation(region_freqs_array, genome_freqs_array)
    
    emp_pval, simulation_distances = _calculate_empirical_pvalue_serial(
        region_freqs_array=region_freqs_array, genome_freqs_array=genome_freqs_array,
        all_cds_features=all_cds_features, codon_order=codon_order, 
        observed_cosine_distance=cosine_distance, num_simulations=num_simulations
    )
    
    region_standardized, genome_standardized = _calculate_standardized_frequencies(region_freqs_array, genome_freqs_array)
    
    # Note: total_cds_features represents the total CDS count across all scaffolds
    return _prepare_final_results(emp_pval, cosine_distance, spearman_rho, region_cds_features, 
                                 all_cds_features, codon_order, region_freqs_array, genome_freqs_array,
                                 region_standardized, genome_standardized, simulation_distances, 
                                 total_cds_count if total_cds_count > 0 else len(all_cds_features))


def _calculate_empirical_pvalue_serial(region_freqs_array: np.ndarray, genome_freqs_array: np.ndarray,
                                      all_cds_features: List[CDSFeature],
                                      codon_order: List[str], observed_cosine_distance: float,
                                      num_simulations: int) -> Tuple[float, List[float]]:
    """ Calculate empirical P-value using serial processing.
    
        Returns:
            tuple of (empirical_pvalue, simulation_distances)
    """
    # Set random seed for deterministic behavior
    np.random.seed(42)
    
    target_codon_count = int(np.sum(region_freqs_array))
    
    # Pre-populate codon count arrays for all CDS features
    from . import _get_cached_codon_counts
    
    # Build codon count arrays
    gene_counts_arrays = []
    gene_total_codons_list = []
    for cds in all_cds_features:
        codon_count_array = _get_cached_codon_counts(cds, codon_order)
        if np.sum(codon_count_array) > 0:  # Only include CDS with codons
            gene_counts_arrays.append(codon_count_array)
            gene_total_codons_list.append(int(np.sum(codon_count_array)))
    
    if not gene_counts_arrays:
        logging.warning("No valid CDS features found for simulations")
        return 1.0, []
    
    gene_total_codons_arr = np.array(gene_total_codons_list, dtype=np.int64)
    
    # Total codon counts vector for background computation
    # IMPORTANT: This includes ALL genes (focal + background)
    total_codon_counts_vec = np.zeros(len(codon_order), dtype=np.float64)
    for gene_counts in gene_counts_arrays:
        total_codon_counts_vec += gene_counts
    
    # Sequential simulation
    emp_pval = 0
    sim_cosine_distances = []
    
    
    for sim in range(num_simulations):
        # Shuffle gene indices
        gene_indices = np.arange(len(gene_counts_arrays))
        np.random.shuffle(gene_indices)
        
        # Initialize simulation focal frequencies
        sim_focal_vec = np.zeros(len(codon_order), dtype=np.float64)
        total_collected = 0
        
        # Sample codons from shuffled genes until target count reached
        for idx in gene_indices:
            cds_total = gene_total_codons_arr[idx]
            if total_collected + cds_total <= target_codon_count:
                sim_focal_vec += gene_counts_arrays[idx]
                total_collected += cds_total
            else:
                remaining = target_codon_count - total_collected
                if remaining > 0 and cds_total > 0:
                    prop = remaining / float(cds_total)
                    sim_focal_vec += gene_counts_arrays[idx] * prop
                break
        
        # Compute cosine distance on counts vectors (focal vs background)
        sim_background_vec = total_codon_counts_vec - sim_focal_vec
        # Add small epsilon to both to match reference behavior
        eps = 1e-10
        sim_cosine_distance = spatial.distance.cosine(sim_focal_vec + eps, sim_background_vec + eps)
        
        sim_cosine_distances.append(sim_cosine_distance)
        if sim_cosine_distance >= observed_cosine_distance:
            emp_pval += 1
    
    empirical_pvalue = (emp_pval + 1) / (num_simulations + 1)
    return empirical_pvalue, sim_cosine_distances


def generate_codoff_histogram(simulation_distances: List[float], observed_distance: float, 
                             output_dir: str, region_number: int, record_id: str = "") -> str:
    """ Generate a histogram of simulation distances with the observed distance marked.
    
        Arguments:
            simulation_distances: list of cosine distances from simulations
            observed_distance: the observed cosine distance
            output_dir: directory to save the plot
            region_number: region number for the filename
            
        Returns:
            path to the generated plot file
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        import os
        
        # Create plots directory if it doesn't exist
        plots_dir = os.path.join(output_dir, "codoff_plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        # Set up the plot with matplotlib styling
        plt.figure(figsize=(12, 8))
        plt.style.use('default')  # Use clean matplotlib style
        plt.gca().set_facecolor('white')
        plt.gca().grid(True, alpha=0.3, color='gray', linestyle='-', linewidth=0.5)
        
        # Create histogram
        if simulation_distances:
            plt.hist(simulation_distances, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
            plt.axvline(x=observed_distance, color='red', linestyle='--', linewidth=2, label=f'Observed: {observed_distance:.3f}')
        else:
            # Handle case with no simulation data
            plt.text(0.5, 0.5, 'No simulation data available', ha='center', va='center', transform=plt.gca().transAxes)
        
        plt.xlabel('Cosine Distance')
        plt.ylabel('Frequency')
        plt.title(f'codoff Simulation Results - Region {region_number}')
        plt.legend()
        
        # Save the plot - try SVG first, fall back to PNG
        # Include record_id to avoid filename collisions across multiple scaffolds
        safe_record_id = record_id.replace("/", "_").replace("\\", "_") if record_id else "unknown"
        plot_path = os.path.join(plots_dir, f"codoff_{safe_record_id}_region_{region_number}.svg")
        
        try:
            plt.savefig(plot_path, format='svg', dpi=300, bbox_inches='tight')
        except Exception as svg_error:
            logging.warning("SVG save failed, trying PNG: %s", svg_error)
            plot_path = os.path.join(plots_dir, f"codoff_{safe_record_id}_region_{region_number}.png")
            plt.savefig(plot_path, format='png', dpi=300, bbox_inches='tight')
        
        plt.close()
        return plot_path
        
    except ImportError as e:
        logging.warning("Could not generate histogram: %s", e)
        return ""
    except Exception as e:
        logging.error("Error generating histogram: %s", e)
        return ""