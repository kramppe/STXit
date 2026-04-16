"""
Prophage Prediction Module
Load and parse PHASTEST prophage regions
"""

from typing import Dict, List, Optional
import csv

class PhagePredictor:
    """Load and analyze prophage prediction results"""
    
    def load_phastest_results(self, phastest_file: str) -> Dict:
        """Load PHASTEST results from file"""
        
        prophage_data = {
            'source': 'PHASTEST',
            'regions': []
        }
        
        try:
            with open(phastest_file, 'r') as f:
                # Try to parse as CSV/TSV
                content = f.read()
                
                if '\t' in content:
                    delimiter = '\t'
                else:
                    delimiter = ','
                
                f.seek(0)
                reader = csv.DictReader(f, delimiter=delimiter)
                
                for row in reader:
                    region = {
                        'region_id': row.get('region_id', ''),
                        'start': int(row.get('start', 0)),
                        'end': int(row.get('end', 0)),
                        'length': int(row.get('length', 0)),
                        'completeness': row.get('completeness', 'unknown'),
                        'score': float(row.get('score', 0)),
                        'insertion_site': row.get('insertion_site', 'unknown'),
                        'phage_name': row.get('phage_name', 'unknown')
                    }
                    prophage_data['regions'].append(region)
        
        except Exception as e:
            print(f"Error loading PHASTEST results: {e}")
            return {}
        
        return prophage_data
    
    def check_stx_in_prophage(self, stx_loci: List, prophage_regions: List[Dict]) -> Dict:
        """Check if STX loci overlap with prophage regions"""
        
        overlaps = {}
        
        for locus in stx_loci:
            locus_overlaps = []
            
            for region in prophage_regions:
                # Check for overlap
                if (locus.start <= region['end'] and locus.end >= region['start']):
                    overlap_info = {
                        'region_id': region['region_id'],
                        'completeness': region['completeness'],
                        'insertion_site': region['insertion_site'],
                        'phage_name': region['phage_name']
                    }
                    locus_overlaps.append(overlap_info)
                    
                    # Mark locus as in prophage
                    locus.in_phage = True
                    locus.prophage_region = f"Region_{region['region_id']}"
                    locus.insertion_site = region['insertion_site']
            
            if not locus_overlaps:
                locus.in_phage = False
                locus.prophage_region = ''
                locus.insertion_site = ''
            
            overlaps[f"STX_{locus.query_id}_{locus.start}"] = locus_overlaps
        
        return overlaps
