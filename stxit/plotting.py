from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import matplotlib.pyplot as plt
from .io_utils import QueryGenome

def _draw_trna_marker(ax, x: int, y: float = 0.62):
    ax.plot([x], [y], marker="^", markersize=6)

def build_closed_genome_plots(queries: List[QueryGenome], phase2, outdir: Path, enabled: bool = False) -> Dict[str, dict]:
    plot_files={}
    if not enabled:
        return plot_files
    phage_df=getattr(phase2, "phage_regions", None)
    overlap_df=getattr(phase2, "stx_phage_overlap", None)
    context_df=getattr(phase2, "phage_insertion_context", None)
    if phage_df is None or phage_df.empty:
        return plot_files
    for q in queries:
        if q.num_contigs != 1:
            continue
        contig=q.contig_order[0]; genome_len=q.contig_lengths[contig]
        regions=phage_df[(phage_df["sample"]==q.sample)&(phage_df["contig"]==contig)]
        if regions.empty: continue
        overlaps=overlap_df[(overlap_df["sample"]==q.sample)&(overlap_df["contig"]==contig)] if overlap_df is not None and not overlap_df.empty else None
        contexts=context_df[(context_df["sample"]==q.sample)&(context_df["contig"]==contig)] if context_df is not None and not context_df.empty else None
        fig, ax=plt.subplots(figsize=(13,3.8))
        ax.hlines(1.0, 1, genome_len, linewidth=3)
        for _, reg in regions.iterrows():
            start=int(reg["start"]); end=int(reg["end"]); width=end-start+1
            ax.broken_barh([(start,width)], (0.82,0.28))
            label=str(reg.get("region_id",""))
            if label: ax.text(start+width/2, 1.17, label, ha="center", va="bottom", fontsize=8)
        if overlaps is not None and not overlaps.empty:
            for _, row in overlaps.iterrows():
                start=int(row["stx_start"]); end=int(row["stx_end"]); width=end-start+1
                ax.broken_barh([(start,width)], (0.42,0.14))
                subtype=str(row.get("assigned_subtype",""))
                if subtype: ax.text(start+width/2, 0.33, subtype, ha="center", va="top", fontsize=8)
        if contexts is not None and not contexts.empty:
            used_positions=[]
            for _, row in contexts.iterrows():
                start=int(row["phage_start"]); end=int(row["phage_end"]); mid=(start+end)/2
                named_site=str(row.get("known_insertion_site_name","") or "")
                backbone_disrupted=str(row.get("backbone_disrupted_feature","") or "")
                left_trna=str(row.get("nearest_trna_left","") or "")
                right_trna=str(row.get("nearest_trna_right","") or "")
                top_label_parts=[]
                if named_site: top_label_parts.append(named_site)
                if backbone_disrupted and backbone_disrupted != named_site: top_label_parts.append(backbone_disrupted)
                top_label=" | ".join([x for x in top_label_parts if x])
                if top_label:
                    y=1.33
                    if used_positions and abs(mid-used_positions[-1]) < genome_len*0.08:
                        y=1.43
                    used_positions.append(mid)
                    ax.text(mid, y, top_label, ha="center", va="bottom", fontsize=8)
                if left_trna:
                    _draw_trna_marker(ax, start, y=0.64); ax.text(start, 0.70, left_trna, ha="right", va="bottom", fontsize=7, rotation=45)
                if right_trna:
                    _draw_trna_marker(ax, end, y=0.64); ax.text(end, 0.70, right_trna, ha="left", va="bottom", fontsize=7, rotation=45)
                left_flank=str(row.get("left_flank_gene","") or row.get("left_flank_locus_tag","") or "")
                right_flank=str(row.get("right_flank_gene","") or row.get("right_flank_locus_tag","") or "")
                backbone_left=str(row.get("backbone_left_gene","") or row.get("backbone_left_locus_tag","") or "")
                backbone_right=str(row.get("backbone_right_gene","") or row.get("backbone_right_locus_tag","") or "")
                if left_flank: ax.text(start, 0.13, f"L:{left_flank}", ha="right", va="top", fontsize=7)
                if right_flank: ax.text(end, 0.13, f"R:{right_flank}", ha="left", va="top", fontsize=7)
                if backbone_left: ax.text(start, 0.03, f"BL:{backbone_left}", ha="right", va="top", fontsize=7)
                if backbone_right: ax.text(end, 0.03, f"BR:{backbone_right}", ha="left", va="top", fontsize=7)
        ax.set_xlim(1, genome_len); ax.set_ylim(-0.05,1.55); ax.set_yticks([])
        ax.set_xlabel(f"{q.sample} genome position (bp)")
        ax.set_title(f"{q.sample}: PHASTEST regions, stx overlap, insertion context")
        for spine in ["left","right","top"]:
            ax.spines[spine].set_visible(False)
        fig.tight_layout()
        png_out=outdir/f"{q.sample}_phage_map.png"; svg_out=outdir/f"{q.sample}_phage_map.svg"
        fig.savefig(png_out, dpi=200); fig.savefig(svg_out); plt.close(fig)
        plot_files[q.sample]={"png":png_out.name,"svg":svg_out.name}
    return plot_files
