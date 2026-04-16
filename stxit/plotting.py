"""
STX Plotting Module
Generate genome maps showing STX loci and prophage context
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from typing import List, Optional, Dict
from Bio import SeqIO
from .stx_typing import LocusCall

class GenomePlotter:
    """Generate genome visualization plots"""
    
    def __init__(self):
        plt.style.use('default')
        
    def plot_genome(self, genome_file: str, loci: List[LocusCall],
                   prophage_data: Optional[Dict] = None,
                   trna_data: Optional[Dict] = None,
                   output_dir: Path = Path('.')):
        """Generate comprehensive genome plot"""
        
        # Load genome
        contigs = {}
        for record in SeqIO.parse(genome_file, 'fasta'):
            contigs[record.id] = len(record)
        
        if not contigs:
            print("Warning: No contigs found in genome")
            return
        
        # Generate plots
        self._plot_linear_map(contigs, loci, prophage_data, trna_data, output_dir)
        
        # Generate circular plot for single chromosome
        if len(contigs) == 1:
            contig_id = list(contigs.keys())[0]
            if contigs[contig_id] > 1000000:  # Only for large chromosomes
                self._plot_circular_map(contig_id, contigs[contig_id], loci, 
                                      prophage_data, trna_data, output_dir)
    
    def _plot_linear_map(self, contigs: Dict[str, int], loci: List[LocusCall],
                        prophage_data: Optional[Dict], trna_data: Optional[Dict],
                        output_dir: Path):
        """Generate linear genome map"""
        
        n_contigs = len(contigs)
        fig, axes = plt.subplots(n_contigs, 1, figsize=(14, 2 * n_contigs))
        
        if n_contigs == 1:
            axes = [axes]
        
        for i, (contig_id, length) in enumerate(contigs.items()):
            ax = axes[i]
            
            # Draw chromosome/contig backbone
            ax.plot([0, length], [0, 0], 'k-', linewidth=3, alpha=0.7)
            
            # Add prophage regions
            if prophage_data and contig_id in prophage_data:
                self._add_prophage_regions(ax, prophage_data[contig_id], length)
            
            # Add tRNA sites
            if trna_data and contig_id in trna_data:
                self._add_trna_sites(ax, trna_data[contig_id])
            
            # Add STX loci
            contig_loci = [locus for locus in loci if locus.query_id == contig_id]
            self._add_stx_loci(ax, contig_loci)
            
            # Formatting
            ax.set_xlim(-length * 0.05, length * 1.05)
            ax.set_ylim(-0.5, 0.5)
            ax.set_xlabel(f'{contig_id} ({length:,} bp)')
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        
        plt.tight_layout()
        
        # Save PNG and SVG
        plt.savefig(output_dir / 'stx_genome_linear.png', dpi=300, bbox_inches='tight')
        plt.savefig(output_dir / 'stx_genome_linear.svg', bbox_inches='tight')
        plt.close()
    
    def _plot_circular_map(self, contig_id: str, length: int, loci: List[LocusCall],
                          prophage_data: Optional[Dict], trna_data: Optional[Dict],
                          output_dir: Path):
        """Generate circular genome map for single chromosome"""
        
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
        
        # Convert coordinates to angles
        theta = [2 * 3.14159 * pos / length for pos in range(0, length, 10000)]
        r = [1] * len(theta)
        
        # Draw chromosome circle
        ax.plot(theta, r, 'k-', linewidth=2, alpha=0.7)
        
        # Add prophage regions
        contig_loci = [locus for locus in loci if locus.query_id == contig_id]
        
        for locus in contig_loci:
            start_angle = 2 * 3.14159 * locus.start / length
            end_angle = 2 * 3.14159 * locus.end / length
            
            # STX locus arc
            angles = [start_angle + (end_angle - start_angle) * i / 50 for i in range(51)]
            radii = [1.1] * len(angles)
            
            color = self._get_subtype_color(locus.subtype)
            ax.plot(angles, radii, color=color, linewidth=6, alpha=0.8)
            
            # Label
            mid_angle = (start_angle + end_angle) / 2
            ax.text(mid_angle, 1.2, locus.subtype, 
                   ha='center', va='center', fontsize=10, weight='bold')
        
        # Formatting
        ax.set_ylim(0, 1.5)
        ax.set_rticks([])
        ax.set_thetagrids(range(0, 360, 30), 
                         [f'{int(i * length / 360):,}' for i in range(0, 360, 30)])
        ax.set_title(f'{contig_id} Circular Map\nSTX Loci Distribution', 
                    pad=20, fontsize=14, weight='bold')
        
        # Legend
        subtypes = list(set(locus.subtype for locus in contig_loci))
        for i, subtype in enumerate(sorted(subtypes)):
            color = self._get_subtype_color(subtype)
            ax.plot([], [], color=color, linewidth=6, label=subtype)
        
        if subtypes:
            ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))
        
        plt.savefig(output_dir / 'stx_genome_circular.png', dpi=300, bbox_inches='tight')
        plt.savefig(output_dir / 'stx_genome_circular.svg', bbox_inches='tight')
        plt.close()
    
    def _add_prophage_regions(self, ax, prophage_regions: List[dict], genome_length: int):
        """Add prophage regions to plot"""
        for region in prophage_regions:
            start = region.get('start', 0)
            end = region.get('end', 0)
            
            # Prophage region rectangle
            rect = patches.Rectangle((start, -0.15), end - start, 0.3,
                                   facecolor='lightblue', alpha=0.6,
                                   edgecolor='blue', linewidth=1)
            ax.add_patch(rect)
            
            # Label
            mid = (start + end) / 2
            ax.text(mid, 0.25, 'Prophage', ha='center', va='bottom',
                   fontsize=8, style='italic')
    
    def _add_trna_sites(self, ax, trna_sites: List[dict]):
        """Add tRNA insertion sites to plot"""
        for site in trna_sites:
            pos = site.get('position', 0)
            trna_type = site.get('type', 'tRNA')
            
            # tRNA marker
            ax.plot([pos, pos], [-0.3, -0.1], 'g-', linewidth=2, alpha=0.8)
            ax.plot(pos, -0.3, 'go', markersize=6, alpha=0.8)
            
            # Label
            ax.text(pos, -0.35, trna_type, ha='center', va='top',
                   fontsize=7, color='green')
    
    def _add_stx_loci(self, ax, loci: List[LocusCall]):
        """Add STX loci to plot"""
        for locus in loci:
            color = self._get_subtype_color(locus.subtype)
            
            # STX locus rectangle
            height = 0.2 if locus.strand == '+' else -0.2
            y_pos = 0.1 if locus.strand == '+' else -0.3
            
            rect = patches.Rectangle((locus.start, y_pos), 
                                   locus.end - locus.start, height,
                                   facecolor=color, alpha=0.8,
                                   edgecolor='black', linewidth=1)
            ax.add_patch(rect)
            
            # Arrow for directionality
            arrow_x = locus.start + (locus.end - locus.start) * 0.8
            arrow_y = y_pos + height / 2
            
            if locus.strand == '+':
                ax.annotate('', xy=(locus.end, arrow_y), xytext=(arrow_x, arrow_y),
                           arrowprops=dict(arrowstyle='->', color='white', lw=2))
            else:
                ax.annotate('', xy=(locus.start, arrow_y), xytext=(arrow_x, arrow_y),
                           arrowprops=dict(arrowstyle='->', color='white', lw=2))
            
            # Label
            mid = (locus.start + locus.end) / 2
            label_y = 0.35 if locus.strand == '+' else -0.45
            ax.text(mid, label_y, locus.subtype, ha='center', va='center',
                   fontsize=10, weight='bold', color=color)
    
    def _get_subtype_color(self, subtype: str) -> str:
        """Get color for STX subtype"""
        color_map = {
            'stx1a': '#FF6B6B',  # Red
            'stx1c': '#FF8E53',  # Orange-red
            'stx1d': '#FF9F43',  # Orange
            'stx2a': '#4ECDC4',  # Teal
            'stx2b': '#45B7D1',  # Blue
            'stx2c': '#96CEB4',  # Green
            'stx2d': '#FFEAA7',  # Yellow
            'stx2e': '#DDA0DD',  # Plum
            'stx2f': '#F8C291',  # Peach
            'stx2g': '#A8E6CF',  # Mint
            'stx2h': '#FFB7B7'   # Pink
        }
        
        return color_map.get(subtype, '#CCCCCC')  # Default gray
