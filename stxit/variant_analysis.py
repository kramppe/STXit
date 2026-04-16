"""
Variant Analysis Module
Codon-aware variant analysis and phylogenetic trees
"""

from typing import List, Dict
from .stx_typing import LocusCall
from .stx_db import STXDatabase

class VariantAnalyzer:
    """Analyze STX sequence variants with codon awareness"""
    
    def __init__(self):
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
    
    def analyze_variants(self, locus: LocusCall, stx_db: STXDatabase):
        """Analyze sequence variants for an STX locus"""
        
        # Get reference sequence
        ref_seq = stx_db.get_sequence(locus.ref_accession)
        
        if not ref_seq:
            return
        
        # Placeholder: Perform alignment and variant calling
        # In real implementation, would use pairwise alignment
        variants = self._identify_variants(locus.query_sequence, ref_seq)
        
        # Analyze codon impact
        analyzed_variants = []
        synonymous_count = 0
        nonsynonymous_count = 0
        
        for variant in variants:
            analysis = self._analyze_codon_impact(variant, ref_seq)
            analyzed_variants.append(analysis)
            
            if analysis['variant_type'] == 'synonymous':
                synonymous_count += 1
            elif analysis['variant_type'] == 'nonsynonymous':
                nonsynonymous_count += 1
        
        # Store results in locus object
        locus.variants = analyzed_variants
        locus.variant_count = len(variants)
        locus.synonymous_count = synonymous_count
        locus.nonsynonymous_count = nonsynonymous_count
    
    def _identify_variants(self, query_seq: str, ref_seq: str) -> List[Dict]:
        """Identify sequence variants between query and reference"""
        variants = []
        
        # Placeholder implementation
        # Real version would use proper alignment (e.g., Biopython pairwise2)
        min_len = min(len(query_seq), len(ref_seq))
        
        for i in range(min_len):
            if query_seq[i] != ref_seq[i]:
                variants.append({
                    'position': i + 1,
                    'ref_nt': ref_seq[i],
                    'query_nt': query_seq[i]
                })
        
        return variants
    
    def _analyze_codon_impact(self, variant: Dict, ref_seq: str) -> Dict:
        """Analyze the impact of a variant on codon and amino acid"""
        
        pos = variant['position'] - 1  # Convert to 0-based
        
        # Find codon position
        codon_start = (pos // 3) * 3
        codon_pos = pos % 3
        
        # Extract reference codon
        ref_codon = ref_seq[codon_start:codon_start + 3].upper()
        
        # Create query codon
        query_codon = ref_codon
        if codon_pos < len(query_codon):
            query_codon = (query_codon[:codon_pos] + 
                          variant['query_nt'].upper() + 
                          query_codon[codon_pos + 1:])
        
        # Translate
        ref_aa = self.genetic_code.get(ref_codon, 'X')
        query_aa = self.genetic_code.get(query_codon, 'X')
        
        # Determine variant type
        if ref_aa == query_aa:
            variant_type = 'synonymous'
        elif query_aa == '*' or ref_aa == '*':
            variant_type = 'nonsense'
        else:
            variant_type = 'nonsynonymous'
        
        return {
            **variant,
            'ref_codon': ref_codon,
            'query_codon': query_codon,
            'ref_aa': ref_aa,
            'query_aa': query_aa,
            'variant_type': variant_type,
            'aa_position': (codon_start // 3) + 1,
            'codon_position': codon_pos + 1
        }
