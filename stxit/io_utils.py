"""
I/O Utilities Module
Query parsing and genome annotation utilities
"""

import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_query_genomes(query_string: str) -> List[Tuple[str, str]]:
    """Parse query genome string into (file, strain_name) pairs
    
    Formats supported:
    - single_file.fasta
    - file1.fasta:strain1,file2.fasta:strain2
    - file.fasta:strain_name
    """
    
    queries = []
    
    if ',' in query_string:
        # Multiple files
        for item in query_string.split(','):
            item = item.strip()
            if ':' in item:
                file_path, strain_name = item.split(':', 1)
                queries.append((file_path.strip(), strain_name.strip()))
            else:
                # Use filename as strain name
                file_path = item
                strain_name = Path(file_path).stem
                queries.append((file_path, strain_name))
    else:
        # Single file
        if ':' in query_string:
            file_path, strain_name = query_string.split(':', 1)
            queries.append((file_path.strip(), strain_name.strip()))
        else:
            file_path = query_string.strip()
            strain_name = Path(file_path).stem
            queries.append((file_path, strain_name))
    
    return queries

def load_genome_sequences(genome_file: str) -> Dict[str, SeqRecord]:
    """Load genome sequences from FASTA file"""
    
    sequences = {}
    
    for record in SeqIO.parse(genome_file, 'fasta'):
        sequences[record.id] = record
    
    return sequences

def find_genome_annotations(genome_file: str, feature_type: str = "CDS") -> Dict[str, List[dict]]:
    """Extract annotations from GenBank file if available
    
    Returns dictionary mapping contig IDs to lists of feature annotations
    """
    
    annotations = {}
    
    # Check if GenBank file exists
    gb_file = Path(genome_file).with_suffix('.gb')
    if not gb_file.exists():
        gb_file = Path(genome_file).with_suffix('.gbk')
    
    if not gb_file.exists():
        return annotations
    
    try:
        for record in SeqIO.parse(gb_file, 'genbank'):
            contig_features = []
            
            for feature in record.features:
                if feature.type == feature_type:
                    feature_info = {
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': feature.location.strand,
                        'locus_tag': feature.qualifiers.get('locus_tag', [''])[0],
                        'protein_id': feature.qualifiers.get('protein_id', [''])[0],
                        'product': feature.qualifiers.get('product', ['hypothetical protein'])[0],
                        'gene': feature.qualifiers.get('gene', [''])[0]
                    }
                    contig_features.append(feature_info)
            
            annotations[record.id] = contig_features
    
    except Exception as e:
        print(f"Warning: Could not parse GenBank file {gb_file}: {e}")
    
    return annotations

def find_overlapping_genes(position: int, contig_id: str, 
                          annotations: Dict[str, List[dict]],
                          tolerance: int = 100) -> List[dict]:
    """Find genes that overlap with a given position"""
    
    overlapping = []
    
    if contig_id not in annotations:
        return overlapping
    
    for gene in annotations[contig_id]:
        gene_start = gene['start']
        gene_end = gene['end']
        
        # Check for overlap with tolerance
        if (gene_start - tolerance <= position <= gene_end + tolerance):
            overlapping.append(gene)
    
    return overlapping

def extract_sequence_region(record: SeqRecord, start: int, end: int, 
                           strand: str = '+') -> str:
    """Extract sequence region from SeqRecord"""
    
    # Ensure coordinates are within bounds
    start = max(0, start)
    end = min(len(record.seq), end)
    
    if start >= end:
        return ""
    
    sequence = str(record.seq[start:end])
    
    if strand == '-':
        from Bio.Seq import Seq
        sequence = str(Seq(sequence).reverse_complement())
    
    return sequence

def format_coordinates(start: int, end: int, strand: str = '+') -> str:
    """Format genomic coordinates for display"""
    
    coord_str = f"{start:,}–{end:,}"
    
    if strand == '+':
        coord_str += " (+)"
    elif strand == '-':
        coord_str += " (−)"
    
    return coord_str

def validate_fasta_file(file_path: str) -> Tuple[bool, str]:
    """Validate that a file is a proper FASTA file"""
    
    if not Path(file_path).exists():
        return False, f"File does not exist: {file_path}"
    
    try:
        records = list(SeqIO.parse(file_path, 'fasta'))
        
        if not records:
            return False, "No sequences found in file"
        
        # Check for valid sequences
        for record in records:
            if len(record.seq) == 0:
                return False, f"Empty sequence found: {record.id}"
            
            # Check for valid nucleotide characters
            valid_chars = set('ATCGNatcgn-')
            seq_chars = set(str(record.seq))
            
            invalid = seq_chars - valid_chars
            if invalid:
                return False, f"Invalid nucleotide characters in {record.id}: {invalid}"
        
        return True, f"Valid FASTA file with {len(records)} sequences"
    
    except Exception as e:
        return False, f"Error reading FASTA file: {str(e)}"

def clean_sequence_id(seq_id: str) -> str:
    """Clean sequence ID for safe use in filenames"""
    
    # Remove problematic characters
    clean_id = re.sub(r'[^a-zA-Z0-9_.-]', '_', seq_id)
    
    # Limit length
    if len(clean_id) > 50:
        clean_id = clean_id[:50]
    
    return clean_id
