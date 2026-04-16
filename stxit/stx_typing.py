"""
STX Typing Module
BLAST-based STX subtype detection and calling
"""

import os
import tempfile
import subprocess
from dataclasses import dataclass
from typing import List, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

@dataclass
class LocusCall:
    """STX locus detection result"""
    query_id: str
    start: int
    end: int
    strand: str
    subtype: str
    identity: float
    coverage: float
    ref_accession: str
    ref_length: int
    query_sequence: str
    
    # Reference database info
    gene_A: str = ""
    gene_B: str = ""
    ref_protein_accession_A: str = ""
    ref_protein_accession_B: str = ""
    nearest_reference_strain: str = ""
    
    # Query genome annotation
    stxA_locus_tag: str = ""
    stxA_protein_id: str = ""
    stxA_product: str = ""
    stxB_locus_tag: str = ""
    stxB_protein_id: str = ""
    stxB_product: str = ""

class STXTyper:
    """STX detection and typing using BLAST"""
    
    def __init__(self, stx_database, min_identity=80.0, min_coverage=70.0):
        self.stx_db = stx_database
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        
    def detect_stx(self, genome_file: str, threads: int = 4) -> List[LocusCall]:
        """Detect STX genes in genome using BLAST"""
        
        # Create temporary BLAST database from genome
        with tempfile.TemporaryDirectory() as temp_dir:
            genome_db = os.path.join(temp_dir, "genome")
            
            # Make BLAST database
            cmd = [
                "makeblastdb",
                "-in", genome_file,
                "-dbtype", "nucl",
                "-out", genome_db
            ]
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Run BLAST search
            blast_out = os.path.join(temp_dir, "blast_results.txt")
            cmd = [
                "blastn",
                "-query", self.stx_db.fasta_path,
                "-db", genome_db,
                "-out", blast_out,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
                "-num_threads", str(threads),
                "-max_target_seqs", "10",
                "-evalue", "1e-10"
            ]
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Parse BLAST results
            hits = self._parse_blast_results(blast_out, genome_file)
            
            # Filter and call loci
            loci = self._call_loci(hits)
            
            return loci
    
    def _parse_blast_results(self, blast_file: str, genome_file: str) -> List[dict]:
        """Parse BLAST tabular output"""
        hits = []
        
        with open(blast_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 14:
                    continue
                
                hit = {
                    'qseqid': fields[0],      # Reference STX gene
                    'sseqid': fields[1],      # Genome contig
                    'pident': float(fields[2]),
                    'length': int(fields[3]),
                    'mismatch': int(fields[4]),
                    'gapopen': int(fields[5]),
                    'qstart': int(fields[6]),
                    'qend': int(fields[7]),
                    'sstart': int(fields[8]),
                    'send': int(fields[9]),
                    'evalue': float(fields[10]),
                    'bitscore': float(fields[11]),
                    'qlen': int(fields[12]),   # Reference length
                    'slen': int(fields[13])    # Contig length
                }
                
                # Calculate coverage
                hit['coverage'] = (hit['length'] / hit['qlen']) * 100
                
                # Apply filters
                if (hit['pident'] >= self.min_identity and 
                    hit['coverage'] >= self.min_coverage and
                    hit['length'] >= 0.7 * hit['qlen']):  # Length filter for spurious hits
                    
                    hits.append(hit)
        
        return hits
    
    def _call_loci(self, hits: List[dict]) -> List[LocusCall]:
        """Convert BLAST hits to STX locus calls"""
        loci = []
        
        for hit in hits:
            # Extract subtype from reference ID
            ref_id = hit['qseqid']
            subtype = self._extract_subtype(ref_id)
            
            # Get reference info
            ref_info = self.stx_db.get_reference_info(ref_id)
            
            # Determine coordinates and strand
            start = min(hit['sstart'], hit['send'])
            end = max(hit['sstart'], hit['send'])
            strand = '+' if hit['sstart'] < hit['send'] else '-'
            
            # Extract query sequence (placeholder - would need genome parsing)
            query_seq = "SEQUENCE_PLACEHOLDER"
            
            locus = LocusCall(
                query_id=hit['sseqid'],
                start=start,
                end=end,
                strand=strand,
                subtype=subtype,
                identity=hit['pident'],
                coverage=hit['coverage'],
                ref_accession=ref_id,
                ref_length=hit['qlen'],
                query_sequence=query_seq,
                gene_A=ref_info.get('gene_A', f'{subtype}A'),
                gene_B=ref_info.get('gene_B', f'{subtype}B'),
                ref_protein_accession_A=ref_info.get('protein_A', ''),
                ref_protein_accession_B=ref_info.get('protein_B', ''),
                nearest_reference_strain=ref_info.get('strain', '')
            )
            
            loci.append(locus)
        
        return loci
    
    def _extract_subtype(self, ref_id: str) -> str:
        """Extract STX subtype from reference ID"""
        # Handle different ID formats
        if '|' in ref_id:
            # Format: accession|subtype|...
            parts = ref_id.split('|')
            if len(parts) >= 2:
                return parts[1]
        
        # Fallback: look for stx pattern
        ref_lower = ref_id.lower()
        for subtype in ['stx1a', 'stx1c', 'stx1d', 'stx2a', 'stx2b', 'stx2c', 
                       'stx2d', 'stx2e', 'stx2f', 'stx2g', 'stx2h']:
            if subtype in ref_lower:
                return subtype
        
        return 'unknown'
