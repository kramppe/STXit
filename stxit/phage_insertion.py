from __future__ import annotations
from dataclasses import dataclass
import pandas as pd
from .io_utils import find_flanking_features, BackboneAnnotation
from .stx_db import KnownSite
from .stx_typing import TypingResults

@dataclass
class Phase2Context:
    phage_regions: pd.DataFrame
    stx_phage_overlap: pd.DataFrame
    phage_insertion_context: pd.DataFrame
    phage_insertion_coords: pd.DataFrame

def _find_nearest_trna(trnas, sample, contig, start, end):
    sub=trnas[(trnas["sample"]==sample)&(trnas["contig"]==contig)].copy()
    if sub.empty:
        return ("","","","")
    left=sub[sub["end"]<start].copy(); right=sub[sub["start"]>end].copy()
    left_row=None if left.empty else left.assign(dist=start-left["end"]).sort_values("dist").iloc[0]
    right_row=None if right.empty else right.assign(dist=right["start"]-end).sort_values("dist").iloc[0]
    left_name=left_dist=right_name=right_dist=""
    if left_row is not None:
        left_name=f"{left_row['trna_type']}({left_row['anticodon']})"; left_dist=int(left_row["dist"])
    if right_row is not None:
        right_name=f"{right_row['trna_type']}({right_row['anticodon']})"; right_dist=int(right_row["dist"])
    return left_name,left_dist,right_name,right_dist

def _match_known_site(left_gene, right_gene, left_locus, right_locus, known_sites):
    if not known_sites:
        return ("","","","")
    candidates={x for x in [left_gene,right_gene,left_locus,right_locus] if x}
    for site in known_sites:
        if site.match_gene and site.match_gene in candidates:
            return site.site_name, site.category, site.backbone_accession, site.backbone_label
        if site.match_locus_tag and site.match_locus_tag in candidates:
            return site.site_name, site.category, site.backbone_accession, site.backbone_label
    return ("","","","")

def _backbone_feature_lookup(backbone: BackboneAnnotation | None, gene: str, locus_tag: str):
    if backbone is None:
        return ("","","","")
    feat=None
    if gene and gene in backbone.gene_index:
        feat=backbone.gene_index[gene]
    elif locus_tag and locus_tag in backbone.locus_index:
        feat=backbone.locus_index[locus_tag]
    if feat is None:
        return ("","","","")
    return feat.gene, feat.locus_tag, feat.feature_type, feat.product

def _infer_backbone_interval(backbone_left, backbone_right):
    left_gene=backbone_left[0] or backbone_left[1]
    right_gene=backbone_right[0] or backbone_right[1]
    if left_gene and right_gene:
        return f"intergenic_between_{left_gene}_and_{right_gene}"
    if left_gene:
        return f"near_{left_gene}"
    if right_gene:
        return f"near_{right_gene}"
    return ""

def build_phase2_context(typing_results: TypingResults, phage_regions, trna_sites, known_sites, compare_backbone="default_k12", backbone: BackboneAnnotation | None = None) -> Phase2Context:
    if phage_regions is None:
        phage_regions=pd.DataFrame(columns=["sample","contig","region_id","start","end","length","completeness","score","attl","attr","note"])
    if trna_sites is None:
        trna_sites=pd.DataFrame(columns=["sample","contig","start","end","trna_type","anticodon","strand"])
    calls_df=pd.DataFrame([c.__dict__ for c in typing_results.calls])
    if calls_df.empty:
        return Phase2Context(phage_regions, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())
    overlaps=[]; contexts=[]; coords=[]
    query_map={q.sample:q for q in typing_results.queries}
    for _, call in calls_df.iterrows():
        regs=phage_regions[(phage_regions["sample"]==call["sample"])&(phage_regions["contig"]==call["contig"])]
        for _, reg in regs.iterrows():
            ov=max(0, min(int(call["end"]), int(reg["end"]))-max(int(call["start"]), int(reg["start"]))+1)
            if ov>0:
                overlaps.append({"sample":call["sample"],"contig":call["contig"],"stx_start":int(call["start"]), "stx_end":int(call["end"]), "assigned_subtype":call["assigned_subtype"], "phage_region_id":reg["region_id"], "phage_start":int(reg["start"]), "phage_end":int(reg["end"]), "phage_completeness":reg["completeness"], "overlap_bp":ov})
    overlap_df=pd.DataFrame(overlaps)
    for _, reg in phage_regions.iterrows():
        sample=str(reg["sample"]); contig=str(reg["contig"]); q=query_map.get(sample)
        left_feat=right_feat=None
        if q:
            left_feat,right_feat=find_flanking_features(q, contig, int(reg["start"]), int(reg["end"]))
        left_name,left_dist,right_name,right_dist=_find_nearest_trna(trna_sites, sample, contig, int(reg["start"]), int(reg["end"]))
        left_gene=left_feat.gene if left_feat else ""; left_locus=left_feat.locus_tag if left_feat else ""; left_product=left_feat.product if left_feat else ""; left_ftype=left_feat.feature_type if left_feat else ""
        right_gene=right_feat.gene if right_feat else ""; right_locus=right_feat.locus_tag if right_feat else ""; right_product=right_feat.product if right_feat else ""; right_ftype=right_feat.feature_type if right_feat else ""
        known_name,known_cat,backbone_acc,backbone_label=_match_known_site(left_gene,right_gene,left_locus,right_locus,known_sites)
        backbone_left=_backbone_feature_lookup(backbone, left_gene, left_locus)
        backbone_right=_backbone_feature_lookup(backbone, right_gene, right_locus)
        backbone_interval=_infer_backbone_interval(backbone_left, backbone_right)
        backbone_disrupted=known_name or backbone_interval or (backbone_left[0] or backbone_right[0] or "")
        region_overlaps=overlap_df[(overlap_df["sample"]==sample)&(overlap_df["phage_region_id"]==reg["region_id"])]
        stx_overlap=not region_overlaps.empty
        stx_subtypes=";".join(sorted(set(region_overlaps["assigned_subtype"].astype(str).tolist()))) if stx_overlap else ""
        insertion_site_class="named_backbone_site" if known_name else ("trna_associated" if ("tRNA" in [left_ftype,right_ftype] or "trna" in [left_ftype.lower(), right_ftype.lower()]) else "intergenic_or_gene_boundary")
        if backbone_interval and insertion_site_class=="intergenic_or_gene_boundary":
            insertion_site_class="backbone_interval_match"
        contexts.append({
            "sample":sample,"contig":contig,"region_id":reg["region_id"],"phage_start":int(reg["start"]),"phage_end":int(reg["end"]), "phage_length":int(reg["length"]), "phage_completeness":reg["completeness"], "stx_overlap":stx_overlap, "stx_subtypes":stx_subtypes,
            "left_flank_feature_type":left_ftype,"left_flank_gene":left_gene,"left_flank_locus_tag":left_locus,"left_flank_product":left_product,
            "right_flank_feature_type":right_ftype,"right_flank_gene":right_gene,"right_flank_locus_tag":right_locus,"right_flank_product":right_product,
            "nearest_trna_left":left_name,"distance_to_left_trna":left_dist,"nearest_trna_right":right_name,"distance_to_right_trna":right_dist,
            "known_insertion_site_match":"yes" if known_name else "no","known_insertion_site_name":known_name,"known_insertion_site_category":known_cat,
            "backbone_accession":backbone_acc,"backbone_label":backbone.label if backbone else (backbone_label or compare_backbone),
            "backbone_left_gene":backbone_left[0],"backbone_left_locus_tag":backbone_left[1],"backbone_left_feature_type":backbone_left[2],"backbone_left_product":backbone_left[3],
            "backbone_right_gene":backbone_right[0],"backbone_right_locus_tag":backbone_right[1],"backbone_right_feature_type":backbone_right[2],"backbone_right_product":backbone_right[3],
            "backbone_disrupted_feature":backbone_disrupted,"backbone_interval":backbone_interval,
            "insertion_site_class":insertion_site_class,"attL":reg.get("attl",""),"attR":reg.get("attr","")
        })
        coords.append({
            "sample":sample,"contig":contig,"phage_start":int(reg["start"]),"phage_end":int(reg["end"]),
            "left_flank_feature":left_gene or left_locus,"right_flank_feature":right_gene or right_locus,
            "left_flank_start":left_feat.start if left_feat else "","left_flank_end":left_feat.end if left_feat else "",
            "right_flank_start":right_feat.start if right_feat else "","right_flank_end":right_feat.end if right_feat else "",
            "backbone_left_feature":backbone_left[0] or backbone_left[1],"backbone_right_feature":backbone_right[0] or backbone_right[1],
            "backbone_disrupted_feature":backbone_disrupted,"backbone_interval":backbone_interval,
            "known_insertion_site_name":known_name,"insertion_site_class":insertion_site_class,
            "nearest_trna_left":left_name,"nearest_trna_right":right_name,"attL":reg.get("attl",""),"attR":reg.get("attr","")
        })
    return Phase2Context(phage_regions, overlap_df, pd.DataFrame(contexts), pd.DataFrame(coords))
