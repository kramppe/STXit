"""
STX Reporting Module
Generate comprehensive output files and summaries
"""

import csv
import json
from pathlib import Path
from typing import List, Optional, Dict
from .stx_typing import LocusCall

class STXReporter:
    """Generate STX analysis reports"""
    
    def __init__(self, prefix: str = "stx"):
        self.prefix = prefix
    
    def write_results(self, loci: List[LocusCall], output_dir: Path,
                     prophage_data: Optional[Dict] = None,
                     trna_data: Optional[Dict] = None):
        """Write all output files"""
        
        # Main results files
        self._write_calls_long(loci, output_dir)
        self._write_summary(loci, output_dir)
        self._write_variant_differences(loci, output_dir)
        
        # Optional analyses
        if any(hasattr(locus, 'indels') and locus.indels for locus in loci):
            self._write_indel_report(loci, output_dir)
        
        # JSON export for programmatic access
        self._write_json_export(loci, output_dir, prophage_data, trna_data)
        
        # Analysis summary log
        self._write_analysis_log(loci, output_dir, prophage_data, trna_data)
    
    def _write_calls_long(self, loci: List[LocusCall], output_dir: Path):
        """Write detailed STX calls with full annotation"""
        
        output_file = output_dir / f"{self.prefix}_calls.long.tsv"
        
        fieldnames = [
            'locus_id', 'contig', 'start', 'end', 'strand', 'length',
            'subtype', 'identity_pct', 'coverage_pct', 'ref_accession',
            'nearest_reference_strain', 
            'gene_A', 'gene_B',
            'ref_protein_accession_A', 'ref_protein_accession_B',
            'stxA_locus_tag', 'stxA_protein_id', 'stxA_product',
            'stxB_locus_tag', 'stxB_protein_id', 'stxB_product',
            'variant_count', 'synonymous_variants', 'nonsynonymous_variants',
            'stx_in_phage', 'prophage_region', 'insertion_site',
            'region_tools'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for i, locus in enumerate(loci, 1):
                row = {
                    'locus_id': f"STX_{i:03d}",
                    'contig': locus.query_id,
                    'start': locus.start,
                    'end': locus.end,
                    'strand': locus.strand,
                    'length': locus.end - locus.start + 1,
                    'subtype': locus.subtype,
                    'identity_pct': f"{locus.identity:.2f}",
                    'coverage_pct': f"{locus.coverage:.2f}",
                    'ref_accession': locus.ref_accession,
                    'nearest_reference_strain': locus.nearest_reference_strain,
                    'gene_A': locus.gene_A,
                    'gene_B': locus.gene_B,
                    'ref_protein_accession_A': locus.ref_protein_accession_A,
                    'ref_protein_accession_B': locus.ref_protein_accession_B,
                    'stxA_locus_tag': locus.stxA_locus_tag,
                    'stxA_protein_id': locus.stxA_protein_id,
                    'stxA_product': locus.stxA_product,
                    'stxB_locus_tag': locus.stxB_locus_tag,
                    'stxB_protein_id': locus.stxB_protein_id,
                    'stxB_product': locus.stxB_product,
                    'variant_count': getattr(locus, 'variant_count', 0),
                    'synonymous_variants': getattr(locus, 'synonymous_count', 0),
                    'nonsynonymous_variants': getattr(locus, 'nonsynonymous_count', 0),
                    'stx_in_phage': getattr(locus, 'in_phage', False),
                    'prophage_region': getattr(locus, 'prophage_region', ''),
                    'insertion_site': getattr(locus, 'insertion_site', ''),
                    'region_tools': getattr(locus, 'region_tools', '')
                }
                
                writer.writerow(row)
    
    def _write_summary(self, loci: List[LocusCall], output_dir: Path):
        """Write concise summary table"""
        
        output_file = output_dir / f"{self.prefix}_summary.tsv"
        
        fieldnames = [
            'locus_id', 'subtype', 'contig', 'coordinates', 'strand',
            'identity_pct', 'status', 'in_phage'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            for i, locus in enumerate(loci, 1):
                # Determine status
                status = "Exact match" if locus.identity >= 99.9 else "Subtype-like variant"
                
                row = {
                    'locus_id': f"STX_{i:03d}",
                    'subtype': locus.subtype,
                    'contig': locus.query_id,
                    'coordinates': f"{locus.start:,}–{locus.end:,}",
                    'strand': '(−)' if locus.strand == '-' else '(+)',
                    'identity_pct': f"{locus.identity:.1f}%",
                    'status': status,
                    'in_phage': getattr(locus, 'in_phage', False)
                }
                
                writer.writerow(row)
    
    def _write_variant_differences(self, loci: List[LocusCall], output_dir: Path):
        """Write variant differences table"""
        
        output_file = output_dir / f"{self.prefix}_variant_differences.tsv"
        
        fieldnames = [
            'locus_id', 'subtype', 'position', 'ref_nt', 'query_nt',
            'ref_codon', 'query_codon', 'ref_aa', 'query_aa',
            'variant_type', 'subunit', 'aa_position'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Placeholder for variant analysis results
            # Would be populated by VariantAnalyzer
            for i, locus in enumerate(loci, 1):
                if hasattr(locus, 'variants'):
                    for variant in locus.variants:
                        row = {
                            'locus_id': f"STX_{i:03d}",
                            'subtype': locus.subtype,
                            'position': variant.get('position', ''),
                            'ref_nt': variant.get('ref_nt', ''),
                            'query_nt': variant.get('query_nt', ''),
                            'ref_codon': variant.get('ref_codon', ''),
                            'query_codon': variant.get('query_codon', ''),
                            'ref_aa': variant.get('ref_aa', ''),
                            'query_aa': variant.get('query_aa', ''),
                            'variant_type': variant.get('type', ''),
                            'subunit': variant.get('subunit', ''),
                            'aa_position': variant.get('aa_position', '')
                        }
                        writer.writerow(row)
    
    def _write_indel_report(self, loci: List[LocusCall], output_dir: Path):
        """Write indel analysis report"""
        
        output_file = output_dir / f"{self.prefix}_indel_report.tsv"
        
        fieldnames = [
            'locus_id', 'subtype', 'indel_type', 'ref_pos', 'query_genome_pos',
            'indel_size', 'indel_sequence', 'codon_impact', 'subunit_affected',
            'zone', 'disruption_consequence',
            'potential_mobile_element', 'potential_mobile_element_reason',
            'tn_blast_run', 'tn_family', 'tn_pct_identity', 'tn_pct_coverage',
            'ref_flank_left', 'ref_flank_right'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Placeholder for indel analysis results
            # Would be populated by Indel analyzer
            for i, locus in enumerate(loci, 1):
                if hasattr(locus, 'indels'):
                    for indel in locus.indels:
                        row = {
                            'locus_id': f"STX_{i:03d}",
                            'subtype': locus.subtype,
                            **indel  # Spread indel analysis results
                        }
                        writer.writerow(row)
    
    def _write_json_export(self, loci: List[LocusCall], output_dir: Path,
                          prophage_data: Optional[Dict] = None,
                          trna_data: Optional[Dict] = None):
        """Write machine-readable JSON export"""
        
        output_file = output_dir / f"{self.prefix}_results.json"
        
        export_data = {
            'analysis_type': 'STXit',
            'version': '1.0.3',
            'stx_loci_count': len(loci),
            'stx_loci': [],
            'prophage_analysis': prophage_data,
            'trna_analysis': trna_data
        }
        
        for i, locus in enumerate(loci, 1):
            locus_data = {
                'locus_id': f"STX_{i:03d}",
                'contig': locus.query_id,
                'start': locus.start,
                'end': locus.end,
                'strand': locus.strand,
                'subtype': locus.subtype,
                'identity': locus.identity,
                'coverage': locus.coverage,
                'reference': {
                    'accession': locus.ref_accession,
                    'strain': locus.nearest_reference_strain,
                    'gene_A': locus.gene_A,
                    'gene_B': locus.gene_B,
                    'protein_A': locus.ref_protein_accession_A,
                    'protein_B': locus.ref_protein_accession_B
                },
                'annotation': {
                    'stxA_locus_tag': locus.stxA_locus_tag,
                    'stxA_product': locus.stxA_product,
                    'stxB_locus_tag': locus.stxB_locus_tag,
                    'stxB_product': locus.stxB_product
                }
            }
            
            # Add optional analysis results
            if hasattr(locus, 'variants'):
                locus_data['variants'] = locus.variants
            if hasattr(locus, 'indels'):
                locus_data['indels'] = locus.indels
            if hasattr(locus, 'in_phage'):
                locus_data['prophage_context'] = {
                    'in_phage': locus.in_phage,
                    'region': getattr(locus, 'prophage_region', ''),
                    'insertion_site': getattr(locus, 'insertion_site', '')
                }
            
            export_data['stx_loci'].append(locus_data)
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
    
    def _write_analysis_log(self, loci: List[LocusCall], output_dir: Path,
                           prophage_data: Optional[Dict] = None,
                           trna_data: Optional[Dict] = None):
        """Write analysis summary log"""
        
        output_file = output_dir / f"{self.prefix}_analysis.log"
        
        with open(output_file, 'w') as f:
            f.write("STXit Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"STX Loci Detected: {len(loci)}\n")
            
            if loci:
                # Subtype summary
                subtypes = {}
                for locus in loci:
                    subtypes[locus.subtype] = subtypes.get(locus.subtype, 0) + 1
                
                f.write("\nSubtype Distribution:\n")
                for subtype, count in sorted(subtypes.items()):
                    f.write(f"  {subtype}: {count}\n")
                
                # Identity distribution
                exact_matches = sum(1 for locus in loci if locus.identity >= 99.9)
                variants = len(loci) - exact_matches
                
                f.write(f"\nSequence Identity:\n")
                f.write(f"  Exact matches (≥99.9%): {exact_matches}\n")
                f.write(f"  Subtype-like variants: {variants}\n")
                
                # Prophage context
                if prophage_data:
                    phage_associated = sum(1 for locus in loci 
                                         if getattr(locus, 'in_phage', False))
                    f.write(f"\nProphage Context:\n")
                    f.write(f"  STX loci in prophage regions: {phage_associated}\n")
                    f.write(f"  STX loci in chromosome: {len(loci) - phage_associated}\n")
                
                # Variant analysis
                total_variants = sum(getattr(locus, 'variant_count', 0) for locus in loci)
                if total_variants > 0:
                    synonymous = sum(getattr(locus, 'synonymous_count', 0) for locus in loci)
                    nonsynonymous = sum(getattr(locus, 'nonsynonymous_count', 0) for locus in loci)
                    
                    f.write(f"\nVariant Analysis:\n")
                    f.write(f"  Total variants: {total_variants}\n")
                    f.write(f"  Synonymous: {synonymous}\n")
                    f.write(f"  Non-synonymous: {nonsynonymous}\n")
            
            f.write("\nAnalysis completed successfully.\n")
