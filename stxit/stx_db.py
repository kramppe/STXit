"""\nSTX Database Module\nFASTA database loader with metadata and citations\n"""

import os
import csv
from pathlib import Path
from typing import Dict, List, Optional
from Bio import SeqIO

class STXDatabase:
    """STX reference database manager"""
    
    def __init__(self, db_path: Optional[str] = None):
        """Initialize STX database
        
        Args:
            db_path: Path to custom STX database (FASTA). If None, uses bundled database.
        """
        if db_path is None:
            # Use bundled database
            module_dir = Path(__file__).parent
            self.fasta_path = module_dir.parent / "databases" / "stx_references.fasta"
            self.metadata_path = module_dir.parent / "databases" / "stx_references.tsv"
        else:
            self.fasta_path = Path(db_path)
            # Look for accompanying metadata file
            self.metadata_path = self.fasta_path.with_suffix('.tsv')
            if not self.metadata_path.exists():
                self.metadata_path = self.fasta_path.with_name('stx_references.tsv')
        
        # Check database exists
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"STX database not found: {self.fasta_path}")
        
        # Load metadata
        self.metadata = self._load_metadata()
        
        # Load sequences
        self.sequences = self._load_sequences()
        
        # Build reference mapping
        self.references = self._build_reference_map()
    
    def _load_metadata(self) -> Dict[str, dict]:
        """Load STX reference metadata from TSV file"""
        metadata = {}
        
        if not self.metadata_path.exists():
            # Return empty metadata if file doesn't exist
            return metadata
        
        with open(self.metadata_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                ref_id = row.get('reference_id', row.get('accession', ''))
                metadata[ref_id] = row
        
        return metadata
    
    def _load_sequences(self) -> Dict[str, str]:
        """Load STX reference sequences from FASTA file"""
        sequences = {}
        
        for record in SeqIO.parse(self.fasta_path, 'fasta'):
            sequences[record.id] = str(record.seq)
        
        return sequences
    
    def _build_reference_map(self) -> Dict[str, dict]:
        """Build mapping of reference IDs to comprehensive info"""
        references = {}
        
        for seq_id in self.sequences:
            ref_info = {
                'sequence': self.sequences[seq_id],
                'length': len(self.sequences[seq_id])
            }
            
            # Add metadata if available
            if seq_id in self.metadata:
                ref_info.update(self.metadata[seq_id])
            
            # Parse header for additional info
            header_info = self._parse_header(seq_id)
            ref_info.update(header_info)
            
            references[seq_id] = ref_info
        
        return references
    
    def _parse_header(self, seq_id: str) -> Dict[str, str]:
        """Parse FASTA header for subtype, gene info, etc."""
        info = {}
        
        # Handle pipe-separated format: accession|subtype|gene|description
        if '|' in seq_id:
            parts = seq_id.split('|')
            if len(parts) >= 2:
                info['subtype'] = parts[1]
            if len(parts) >= 3:
                info['gene'] = parts[2]
            if len(parts) >= 4:
                info['description'] = parts[3]
        else:
            # Plain subtype name (e.g. stx1a, stx2c) — use directly
            info['subtype'] = seq_id.split('_')[0] if '_' in seq_id else seq_id
        
        return info
    
    def get_reference_info(self, ref_id: str) -> Dict[str, str]:
        """Get comprehensive reference information"""
        return self.references.get(ref_id, {})
    
    def get_sequence(self, ref_id: str) -> Optional[str]:
        """Get sequence for reference ID"""
        return self.sequences.get(ref_id)
    
    def list_subtypes(self) -> set:
        """Get set of all STX subtypes in database"""
        subtypes = set()
        
        for ref_info in self.references.values():
            if 'subtype' in ref_info:
                subtypes.add(ref_info['subtype'])
        
        return subtypes
    
    def get_subtype_references(self, subtype: str) -> Dict[str, dict]:
        """Get all references for a specific subtype"""
        subtype_refs = {}
        
        for ref_id, ref_info in self.references.items():
            if ref_info.get('subtype') == subtype:
                subtype_refs[ref_id] = ref_info
        
        return subtype_refs
    
    def get_protein_accessions(self, subtype: str) -> Dict[str, str]:
        """Get protein accessions for A and B subunits"""
        accessions = {'A': '', 'B': ''}
        
        subtype_refs = self.get_subtype_references(subtype)
        
        for ref_id, ref_info in subtype_refs.items():
            gene = ref_info.get('gene', '').upper()
            if 'A' in gene:
                accessions['A'] = ref_info.get('protein_accession', '')
            elif 'B' in gene:
                accessions['B'] = ref_info.get('protein_accession', '')
        
        return accessions
    
    def validate_database(self) -> List[str]:
        """Validate database integrity and return any issues"""
        issues = []
        
        # Check for required subtypes
        required_subtypes = [
            'stx1a', 'stx1c', 'stx1d',
            'stx2a', 'stx2b', 'stx2c', 'stx2d', 'stx2e', 'stx2f', 'stx2g', 'stx2h'
        ]
        
        found_subtypes = self.list_subtypes()
        missing = set(required_subtypes) - found_subtypes
        
        if missing:
            issues.append(f"Missing subtypes: {', '.join(sorted(missing))}")
        
        # Check sequence lengths (nucleotide sequences: 200–5000 bp)
        for ref_id, ref_info in self.references.items():
            seq_len = ref_info.get('length', 0)
            if seq_len < 200:
                issues.append(f"Short sequence ({seq_len} bp): {ref_id}")
            elif seq_len > 5000:
                issues.append(f"Long sequence ({seq_len} bp): {ref_id}")
        
        return issues
