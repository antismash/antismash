# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" HTML output generation for the codoff module. """

import logging
import os
from typing import List
from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.secmet import Record
from antismash.common.secmet.features import Region
from antismash.config import ConfigType
from .analysis import generate_codoff_histogram
from . import CodoffResults


def will_handle(_products: List[str], _categories: List[str]) -> bool:
    """ Returns True if this module can handle the given products/categories """
    return True


def generate_html(region: Region, results: CodoffResults, record: Record, options: ConfigType) -> HTMLSections:
    """ Generate HTML content for the codoff module """
    html_sections = HTMLSections("codoff")
    
    # Safety check: ensure we have results and the region exists
    if not results or not hasattr(results, 'by_region') or region.get_region_number() not in results.by_region:
        return html_sections
    
    region_result = results.by_region[region.get_region_number()]
    
    # Safety check: ensure region_result is valid
    if not region_result or not isinstance(region_result, dict):
        logging.warning("Invalid region result for region %d, skipping HTML generation", region.get_region_number())
        return html_sections
    
    # Generate histogram if we have simulation data
    output_dir = options.output_dir
    region_number = region.get_region_number()
    
    simulation_distances = region_result.get('simulation_distances')
    observed_distance = region_result.get('cosine_distance')
    
    if simulation_distances and observed_distance is not None:
        histogram_path = generate_codoff_histogram(
            simulation_distances=simulation_distances,
            observed_distance=observed_distance,
            output_dir=output_dir,
            region_number=region_number,
            record_id=record.id
        )
        
        if histogram_path:
            # Make path relative to output directory for HTML
            try:
                relative_path = os.path.relpath(histogram_path, output_dir)
                region_result['histogram_path'] = relative_path
            except ValueError as e:
                logging.warning("Could not calculate relative path for histogram: %s", e)
                region_result['histogram_path'] = histogram_path
        else:
            region_result['histogram_path'] = ""
            logging.warning("Failed to generate histogram for region %d", region_number)
    else:
        region_result['histogram_path'] = ""
    
    # Format results for template
    formatted_results = {
        'empirical_pvalue': region_result.get('empirical_pvalue'),
        'cosine_distance': region_result.get('cosine_distance', 1.0),
        'spearman_rho': region_result.get('spearman_rho', 0.0),
        'region_cds_count': region_result.get('region_cds_count', 0),
        'genome_cds_count_excluding_region': region_result.get('genome_cds_count_excluding_region', 0),
        'total_cds_count': region_result.get('total_cds_count', 0),
        'codon_order': region_result.get('codon_order', []),
        'region_codon_frequencies': region_result.get('region_codon_frequencies', []),
        'genome_codon_frequencies': region_result.get('genome_codon_frequencies', []),
        'region_codon_frequencies_standardized': region_result.get('region_codon_frequencies_standardized', []),
        'genome_codon_frequencies_standardized': region_result.get('genome_codon_frequencies_standardized', []),
        'region_total_codons': region_result.get('region_total_codons', 0),
        'genome_total_codons': region_result.get('genome_total_codons', 0),
        'histogram_path': region_result.get('histogram_path', "")
    }
    
    # Create the HTML content
    template = FileTemplate(os.path.join(os.path.dirname(__file__), "templates", "codoff.html"))
    html_content = template.render(results=formatted_results)
    
    # Add the HTML section
    html_sections.add_detail_section(
        "codoff Analysis",
        html_content,
        "codoff"
    )
    
    return html_sections