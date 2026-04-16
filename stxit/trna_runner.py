"""
tRNA Scanner Module
Run and parse tRNAscan-SE for insertion site analysis
"""

import subprocess
import tempfile
import os
from typing import Dict, List

class tRNARunner:
    """Run tRNAscan-SE and parse results"""
    
    def __init__(self):
        self.trnascan_cmd = "tRNAscan-SE"
    
    def scan_genome(self, genome_file: str) -> Dict:
        """Run tRNAscan-SE on genome"""
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = os.path.join(temp_dir, "trna_results.txt")
            
            # Run tRNAscan-SE
            cmd = [
                self.trnascan_cmd,
                "-B",  # Bacterial mode
                "-o", output_file,
                genome_file
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print("tRNAscan-SE completed successfully")
                
                # Parse results
                return self._parse_trnascan_output(output_file)
                
            except subprocess.CalledProcessError as e:
                print(f"tRNAscan-SE failed: {e}")
                return {}
            except FileNotFoundError:
                print("tRNAscan-SE not found. Please install tRNAscan-SE.")
                return {}
    
    def load_results(self, results_file: str) -> Dict:
        """Load pre-computed tRNAscan-SE results"""
        
        if not os.path.exists(results_file):
            raise FileNotFoundError(f"tRNA results file not found: {results_file}")
        
        return self._parse_trnascan_output(results_file)
    
    def _parse_trnascan_output(self, output_file: str) -> Dict:
        """Parse tRNAscan-SE output file"""
        
        trna_data = {}
        
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines
        data_lines = []
        for line in lines:
            if not line.startswith('#') and not line.startswith('Sequence') and line.strip():
                data_lines.append(line.strip())
        
        # Parse tRNA entries
        for line in data_lines[2:]:  # Skip first two lines (headers)
            if not line.strip():
                continue
            
            fields = line.split()
            if len(fields) >= 9:
                contig = fields[0]
                trna_num = fields[1]
                start = int(fields[2])
                end = int(fields[3])
                aa = fields[4]
                anticodon = fields[5]
                score = float(fields[8]) if fields[8] != '---' else 0.0
                
                if contig not in trna_data:
                    trna_data[contig] = []
                
                trna_data[contig].append({
                    'trna_id': f"tRNA-{aa}{trna_num}",
                    'start': start,
                    'end': end,
                    'amino_acid': aa,
                    'anticodon': anticodon,
                    'score': score,
                    'type': f"tRNA-{aa}"
                })
        
        return trna_data
