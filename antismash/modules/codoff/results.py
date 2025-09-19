# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Results class for the codoff module. """

from typing import Any, Dict
from antismash.common.secmet import Record
from antismash.common.module_results import ModuleResults
from antismash.config import ConfigType


class CodoffResults(ModuleResults):
    """ Results class for codoff analysis. """
    
    schema_version = 1
    
    def __init__(self, record_id: str):
        super().__init__(record_id)
        self.by_region: Dict[int, Dict[str, Any]] = {}
    
    def to_json(self) -> Dict[str, Any]:
        """ Convert results to JSON-serializable format. """
        # Create a copy of by_region with genome_codon_frequencies excluded
        by_region_for_json = {}
        for region_num, region_data in self.by_region.items():
            if isinstance(region_data, dict):
                # Copy all fields except genome_codon_frequencies
                region_data_copy = {k: v for k, v in region_data.items() if k != 'genome_codon_frequencies'}
                by_region_for_json[region_num] = region_data_copy
            else:
                by_region_for_json[region_num] = region_data
        
        return {
            'schema_version': self.schema_version,
            'record_id': self.record_id,
            'by_region': by_region_for_json
        }
    
    @staticmethod
    def from_json(json_data: Dict[str, Any], record: Record) -> 'CodoffResults':
        """ Create CodoffResults from JSON data. """
        results = CodoffResults(json_data['record_id'])
        results.by_region = json_data.get('by_region', {})
        return results
    
    def add_to_record(self, record: Record) -> None:
        """ Add results to the record. """
        # Results are stored in the module results, not directly on the record
        pass
    
    def write_outputs(self, record: Record, options: ConfigType) -> None:
        """ Write output files and trigger final analysis """
        # Import here to avoid circular imports
        from . import _finalize_analysis_for_all_records, _codoff_global_state
        
        # Only run finalization once (on first write_outputs call)
        if not _codoff_global_state.get('finalization_complete', False):
            import logging
            _finalize_analysis_for_all_records(options)
            _codoff_global_state['finalization_complete'] = True
