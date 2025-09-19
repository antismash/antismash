# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" codon usage assessment """

import logging
from typing import Any, Dict, List, Optional
from collections import Counter
from antismash.common.secmet import Record, CDSFeature
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs
from .results import CodoffResults
from .html_output import generate_html, will_handle

# Import scipy functions for statistics
try:
    from scipy.spatial.distance import cosine
    import numpy as np
except ImportError:
    # Fallback if scipy is not available
    cosine = lambda x, y: 0.5
    np = None

NAME = "codoff"
SHORT_DESCRIPTION = "codon usage assessment"

# Global state for multi-scaffold genome analysis
_codoff_global_state = {
    'total_genome_size': 0,
    'all_cds_features': [],
    'cds_to_record_map': {},  # Map CDS features to their record IDs
    'cds_to_record_ref_map': {},  # Map CDS features to their record references
    'current_record_count': 0,
    'finalization_complete': False,  # Flag to prevent multiple finalizations
    'bgc_regions_tracked': []  # Track all BGC regions for finalization
}

# Global codon cache to avoid repeated BioPython extractions
_codon_extraction_cache = {}

# Global cache for pre-computed codon count arrays
_codon_counts_cache = {}


def _get_cached_codons(cds: CDSFeature, record_ref: Record = None) -> List[str]:
    """ Extract codons from CDS feature with caching """
    cds_id = id(cds)
    
    # Return cached result if available
    if cds_id in _codon_extraction_cache:
        return _codon_extraction_cache[cds_id]
    
    try:
        # Get record reference
        if record_ref is None:
            if cds_id in _codoff_global_state.get('cds_to_record_ref_map', {}):
                record_ref = _codoff_global_state['cds_to_record_ref_map'][cds_id]
            elif hasattr(cds, 'record') and cds.record is not None:
                record_ref = cds.record
            else:
                _codon_extraction_cache[cds_id] = []
                return []
        
        # Extract sequence (this is the expensive BioPython operation)
        sequence = str(cds.location.extract(record_ref.seq)).upper()
        
        # Fast validation and codon extraction
        if len(sequence) % 3 != 0:
            _codon_extraction_cache[cds_id] = []
            return []
        
        # Vectorized codon extraction with validation
        codons = []
        valid_bases = set('ACGT')
        invalid_codon_count = 0
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]
            if len(codon) == 3 and all(base in valid_bases for base in codon):
                codons.append(codon)
            else:
                invalid_codon_count += 1
        
        
        # Cache the result
        _codon_extraction_cache[cds_id] = codons
        return codons
        
    except Exception:
        # Cache empty result for failed extractions
        _codon_extraction_cache[cds_id] = []
        return []


def _get_cached_codon_counts(cds: CDSFeature, codon_order: List[str], record_ref: Record = None) -> np.ndarray:
    """ Get codon count array for a CDS feature """
    cds_id = id(cds)
    
    # Check if already cached
    if cds_id in _codon_counts_cache:
        return _codon_counts_cache[cds_id]
    
    # Get the codon list (this will populate _codon_extraction_cache if needed)
    codons = _get_cached_codons(cds, record_ref)
    
    if not codons:
        # Cache empty array for failed extractions
        empty_array = np.zeros(len(codon_order), dtype=np.int32)
        _codon_counts_cache[cds_id] = empty_array
        return empty_array
    
    # Count codons and create array in codon_order
    from collections import Counter
    codon_counter = Counter(codons)
    codon_counts_array = np.array([codon_counter.get(codon, 0) for codon in codon_order], dtype=np.int32)
    
    # Cache the result
    _codon_counts_cache[cds_id] = codon_counts_array
    return codon_counts_array


def get_arguments() -> ModuleArgs:
    """ Build and return arguments for the module """
    args = ModuleArgs("codoff options", "codoff", enabled_by_default=False)
    args.add_analysis_toggle("enabled",
                             dest="enabled",
                             action="store_true",
                             default=False,
                             help="Run codon usage analysis for BGCs")
    args.add_option("simulations",
                    dest="simulations",
                    type=int,
                    default=10000,
                    help="Number of simulations for empirical P-value calculation (default: %(default)s)")
    args.add_option("min-genome-size",
                    dest="min_genome_size",
                    type=int,
                    default=500000,
                    help="Minimum genome size for codoff analysis (default: %(default)s)")
    args.add_option("max-genome-size",
                    dest="max_genome_size",
                    type=int,
                    default=50000000,
                    help="Maximum genome size for codoff analysis (default: %(default)s)")
    return args


def check_options(options: ConfigType) -> List[str]:
    """ Check module options for conflicts """
    errors = []
    if options.codoff_simulations < 100:
        errors.append("Number of simulations must be at least 100")
    if options.codoff_min_genome_size < 1000:
        errors.append("Minimum genome size must be at least 1,000 bp")
    if options.codoff_max_genome_size < options.codoff_min_genome_size:
        errors.append("Maximum genome size must be greater than minimum genome size")
    return errors


def prepare_data(_logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared """
    global _codoff_global_state, _codon_extraction_cache, _codon_counts_cache
    
    logging.info("Codoff: prepare_data called - initializing global state")
    
    # Reset global state for fresh run
    _codoff_global_state.clear()
    
    # Clear both caches to free memory from previous runs
    _codon_extraction_cache.clear()
    _codon_counts_cache.clear()
    _codoff_global_state.update({
        'total_genome_size': 0,
        'all_cds_features': [],
        'cds_to_record_map': {},
        'cds_to_record_ref_map': {},
        'current_record_count': 0,
        'finalization_complete': False,
        'finalized_results': {},  # Store finalized results for each record
        'bgc_regions_tracked': []  # Track all BGC regions for finalization
    })
    
    logging.info("Codoff: Global state initialized with %d keys", len(_codoff_global_state))
    
    # No external data required for codoff analysis
    return []


def check_prereqs(_options: ConfigType) -> List[str]:
    """ Check prerequisites are available """
    errors = []
    
    # Check required Python packages
    try:
        import matplotlib
    except ImportError:
        errors.append("matplotlib is required for plotting")
    
    try:
        import numpy
    except ImportError:
        errors.append("numpy is required for numerical operations")
    
    try:
        import scipy
    except ImportError:
        errors.append("scipy is required for statistical calculations")
    
    return errors


def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled by the given options """
    return options.codoff_enabled


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[CodoffResults]:
    """ Rebuild previous results from JSON """
    if not previous:
        return None
    return CodoffResults.from_json(previous, record)


def _ensure_json_serializable(value: Any) -> Any:
    """ Ensure a value is JSON serializable """
    if hasattr(value, 'item'):  # numpy scalar
        return value.item()
    elif hasattr(value, 'tolist'):  # numpy array
        return value.tolist()
    elif isinstance(value, (list, tuple)):
        return [_ensure_json_serializable(v) for v in value]
    elif isinstance(value, dict):
        return {k: _ensure_json_serializable(v) for k, v in value.items()}
    else:
        return value


def _calculate_region_codon_frequencies(region_cds_features: List[CDSFeature]) -> Dict[str, int]:
    """ Calculate codon frequencies for a specific region """
    # Use vectorized counting for processing
    all_codons = []
    cache_hits = 0
    cache_misses = 0
    
    for cds in region_cds_features:
        try:
            cds_id = id(cds)
            if cds_id in _codon_extraction_cache:
                cache_hits += 1
                cds_codons = _codon_extraction_cache[cds_id]
            else:
                cache_misses += 1
                # Get cached codons (should populate cache)
                cds_codons = _get_cached_codons(cds)
            
            all_codons.extend(cds_codons)
        except Exception:
            continue
    
    
    # Vectorized counting is much faster than manual loops
    return dict(Counter(all_codons))


def _calculate_genome_codon_frequencies(genome_cds_features: List[CDSFeature]) -> Dict[str, int]:
    """ Calculate codon frequencies for genome-wide CDS features """
    # Use vectorized counting for processing
    all_codons = []
    cache_hits = 0
    cache_misses = 0
    
    for cds in genome_cds_features:
        try:
            cds_id = id(cds)
            if cds_id in _codon_extraction_cache:
                cache_hits += 1
                cds_codons = _codon_extraction_cache[cds_id]
            else:
                cache_misses += 1
                # Get cached codons (should populate cache)
                cds_codons = _get_cached_codons(cds)
            
            all_codons.extend(cds_codons)
        except Exception:
            continue
    
    if cache_hits + cache_misses > 0:
        cache_hit_rate = cache_hits / (cache_hits + cache_misses) * 100
        logging.info("Codoff: Genome codon freq cache hit rate: %.1f%% (%d hits, %d misses, total cache size: %d)", 
                    cache_hit_rate, cache_hits, cache_misses, len(_codon_extraction_cache))
    
    # Vectorized counting is much faster than manual loops
    return dict(Counter(all_codons))


def _run_final_analysis_for_region(bgc_info: Dict[str, Any], region_codon_freqs: Dict[str, int], 
                                  genome_codon_freqs: Dict[str, int], final_total_cds: int, 
                                  options: ConfigType) -> None:
    """ Run the final codoff analysis for a single region """
    record_id = bgc_info['record_id']
    region_number = bgc_info['region_number']
    region_cds_features = bgc_info['region_cds_features']
    
    try:
        # Import here to avoid circular imports
        from .analysis import run_codoff_analysis
        
        # Get options from stored bgc_info
        stored_options = bgc_info.get('options', {})
        num_simulations = stored_options.get('num_simulations', options.codoff_simulations)
        num_cpus = stored_options.get('num_cpus', options.cpus)
        
        # Calculate final genome CDS features (excluding this region)
        final_genome_cds_features = [cds for cds in _codoff_global_state['all_cds_features']
                                   if not (cds in region_cds_features)]
        
        
        final_analysis = run_codoff_analysis(
            region_codon_freqs=region_codon_freqs,
            genome_codon_freqs=genome_codon_freqs,
            region_cds_features=region_cds_features,
            all_cds_features=final_genome_cds_features,
            num_simulations=num_simulations,
            total_cds_count=len(_codoff_global_state['all_cds_features'])
        )
        
        # Store the finalized results AND update the results object directly
        bgc_info['finalized_analysis'] = final_analysis
        
        # Update the results object directly
        if 'results_ref' in bgc_info and bgc_info['results_ref'] is not None:
            bgc_info['results_ref'].by_region[region_number] = final_analysis
        
    except Exception as e:
        logging.error("codoff: Error finalizing analysis for %s region %d: %s", record_id, region_number, e)
        # Store error info for debugging
        bgc_info['finalization_error'] = str(e)


def _finalize_analysis_for_all_records(options: ConfigType) -> None:
    """ Finalize the codoff analysis for all records """
    global _codoff_global_state
    
    if not _codoff_global_state.get('all_cds_features'):
        logging.warning("No CDS features collected, nothing to finalize")
        return
    
    final_total_cds = len(_codoff_global_state['all_cds_features'])
    # Count how many CDS features actually contribute codons across all records
    contributing_cds_total = 0
    for cds in _codoff_global_state['all_cds_features']:
        record_ref = _codoff_global_state['cds_to_record_ref_map'].get(id(cds))
        if record_ref:
            codons = _get_cached_codons(cds, record_ref)
            if len(codons) > 0:
                contributing_cds_total += 1
    
    # Process all tracked BGC regions to complete the analysis
    for bgc_info in _codoff_global_state['bgc_regions_tracked']:
        record_id = bgc_info['record_id']
        region_number = bgc_info['region_number']
        
        # Get the region CDS features
        region_cds_features = bgc_info['region_cds_features']
        
        # Calculate final genome-wide codon frequencies EXCLUDING this region
        # This uses the complete global CDS database
        final_genome_cds_features = [cds for cds in _codoff_global_state['all_cds_features']
                                   if not (cds in region_cds_features)]
        
        # Calculate codon frequencies using helper functions
        region_codon_freqs = _calculate_region_codon_frequencies(region_cds_features)
        
        genome_codon_freqs = _calculate_genome_codon_frequencies(final_genome_cds_features)
        
        
        # Run the final analysis for this region
        _run_final_analysis_for_region(bgc_info, region_codon_freqs, genome_codon_freqs, 
                                     final_total_cds, options)
    
    # Update global state
    _codoff_global_state['final_total_cds_count'] = final_total_cds
    _codoff_global_state['global_cds_collection_complete'] = True


def _collect_global_cds_info(record: Record, options: ConfigType) -> None:
    """ Collect CDS information from a record and add to global state """
    global _codoff_global_state
    
    # Update global state
    _codoff_global_state['total_genome_size'] += len(record.seq)
    _codoff_global_state['current_record_count'] += 1

    # Collect all CDS features from this record
    # Match standalone codoff: process all CDS features, filter during codon extraction
    record_cds_features = record.get_cds_features()
    valid_cds_features = record_cds_features  # No pre-filtering, let codon extraction handle it
    
    
    # Store CDS features with their record ID for later identification
    for cds in valid_cds_features:
        _codoff_global_state['cds_to_record_map'][id(cds)] = record.id
        _codoff_global_state['cds_to_record_ref_map'][id(cds)] = record
    _codoff_global_state['all_cds_features'].extend(valid_cds_features)
    

def run_on_record(record: Record, results: Optional[CodoffResults],
                  options: ConfigType) -> CodoffResults:
    """ Run codoff analysis on a record """
    global _codoff_global_state
    
    # Ensure global state is properly initialized
    if 'bgc_regions_tracked' not in _codoff_global_state:
        _codoff_global_state['bgc_regions_tracked'] = []
    if 'final_total_cds_count' not in _codoff_global_state:
        _codoff_global_state['final_total_cds_count'] = None

    if results is None:
        results = CodoffResults(record.id)

    # Collect genome-wide CDS information from this record
    _collect_global_cds_info(record, options)
    
    # Check genome size requirements using total genome size
    total_genome_size = _codoff_global_state['total_genome_size']
    if total_genome_size < options.codoff_min_genome_size:
        logging.warning("Total genome size too small for codoff analysis (%d bp < %d bp). Skipping.", 
                       total_genome_size, options.codoff_min_genome_size)
        return results
    if total_genome_size > options.codoff_max_genome_size:
        logging.warning("Total genome size too large for codoff analysis (%d bp > %d bp). Skipping.", 
                       total_genome_size, options.codoff_max_genome_size)
        return results
    
    # Store BGC region information for later analysis
    all_cds_features = _codoff_global_state['all_cds_features']
    cds_to_record_map = _codoff_global_state['cds_to_record_map']
    
    if not all_cds_features:
        logging.warning("No CDS features found that meet minimum length requirement")
        return results

    for region in record.get_regions():
        region_number = region.get_region_number()
        
        # Skip if already tracked
        region_key = f"{record.id}_{region_number}"
        existing_regions = [bgc['record_id'] + '_' + str(bgc['region_number']) 
                           for bgc in _codoff_global_state.get('bgc_regions_tracked', [])]
        if region_key in existing_regions:
            continue
            
        # Get CDS features in this region
        # Use region.cds_children to get CDS features (found this works in debug logs)
        
        # Use region.get_cds_features() if available, otherwise fall back to overlap detection
        if hasattr(region, 'get_cds_features'):
            region_cds_features = region.get_cds_features()
        elif hasattr(region, 'cds_children'):
            region_cds_features = region.cds_children
        else:
            # Fall back to overlap detection with stricter criteria
            # Only include CDS features that have significant overlap with the region
            region_cds_features = []
            for cds in all_cds_features:
                if cds_to_record_map.get(id(cds)) != record.id:
                    continue
                if not cds.overlaps_with(region):
                    continue
                    
                # Calculate overlap percentage
                cds_start = cds.location.start
                cds_end = cds.location.end
                region_start = region.location.start
                region_end = region.location.end
                
                # Find overlap boundaries
                overlap_start = max(cds_start, region_start)
                overlap_end = min(cds_end, region_end)
                overlap_length = max(0, overlap_end - overlap_start)
                cds_length = cds_end - cds_start
                
                # Include CDS if it has >50% overlap with the region
                if cds_length > 0 and (overlap_length / cds_length) > 0.5:
                    region_cds_features.append(cds)
            
        if not region_cds_features:
            logging.warning("No CDS features found in region %d", region_number)
            continue
           
        # Track this BGC region for final analysis - ONLY essential data
        bgc_info = {
            'record_id': record.id,
            'region_number': region_number,
            'region_start': region.location.start,
            'region_end': region.location.end,
            'region_cds_features': region_cds_features,
            'results_ref': results,  # Store reference to results object for updating
            'options': {
                'num_simulations': options.codoff_simulations,
                'num_cpus': options.cpus
            }
        }
        
        # Ensure bgc_regions_tracked exists
        if 'bgc_regions_tracked' not in _codoff_global_state:
            _codoff_global_state['bgc_regions_tracked'] = []
        
        _codoff_global_state['bgc_regions_tracked'].append(bgc_info)

    # Statistical analysis will be performed by write_outputs() after all records are processed
    return results