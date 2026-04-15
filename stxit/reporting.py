from __future__ import annotations
from pathlib import Path
from typing import Any, Dict
import json, platform
from datetime import datetime
import pandas as pd

def _calls_df(typing_results):
    rows=[]
    for c in typing_results.calls:
        rows.append({"sample":c.sample,"assembly_type":c.assembly_type,"contig":c.contig,"contig_num":c.contig_num,"start":c.start,"end":c.end,"strand":c.strand,"stx_family":c.stx_family,"assigned_subtype":c.assigned_subtype,"nearest_reference_id":c.nearest_reference_id,"nearest_reference_accession":c.nearest_reference_accession,"nearest_reference_strain":c.nearest_reference_strain,"nearest_reference_citation":c.nearest_reference_citation,"pct_identity":round(c.pct_identity,3),"pct_coverage":round(c.pct_coverage,3),"hit_length":c.hit_length,"complete_partial":c.complete_partial,"exact_match":c.exact_match,"variant_status":c.variant_status,"annotation_source":c.annotation_source,"query_feature_type":c.query_feature_type,"query_gene":c.query_gene,"query_locus_tag":c.query_locus_tag,"query_product":c.query_product})
    return pd.DataFrame(rows)

def write_outputs(typing_results, db, phase2, phase3, plot_files, phastest_jobs, trna_runs, outdir: Path, cli_args: Dict[str, Any], version: str):
    outdir.mkdir(parents=True, exist_ok=True)
    calls_df=_calls_df(typing_results)
    calls_df.to_csv(outdir/"stx_calls.long.tsv", sep="\t", index=False)
    phase2.phage_regions.to_csv(outdir/"phage_regions.tsv", sep="\t", index=False)
    phase2.stx_phage_overlap.to_csv(outdir/"stx_phage_overlap.tsv", sep="\t", index=False)
    phase2.phage_insertion_context.to_csv(outdir/"phage_insertion_context.tsv", sep="\t", index=False)
    phase2.phage_insertion_coords.to_csv(outdir/"phage_insertion_coords.tsv", sep="\t", index=False)
    phase3.variant_differences.to_csv(outdir/"stx_variant_differences.tsv", sep="\t", index=False)
    phase3.alignment_overview.to_csv(outdir/"stx_alignment_overview.tsv", sep="\t", index=False)
    if phastest_jobs: pd.DataFrame(phastest_jobs).to_csv(outdir/"phastest_jobs.tsv", sep="\t", index=False)
    if trna_runs: pd.DataFrame(trna_runs).to_csv(outdir/"trnascan_runs.tsv", sep="\t", index=False)

    tree_index=[]
    for k, v in phase3.tree_fasta_paths.items():
        sample, subtype = k.split("::",1)
        tree_index.append({"sample":sample,"subtype":subtype,"tree_fasta":v,"tree_newick":phase3.tree_newick_paths.get(k,"")})
    pd.DataFrame(tree_index).to_csv(outdir/"stx_variant_tree.index.tsv", sep="\t", index=False)

    plot_index=[]
    for sample, files in plot_files.items():
        plot_index.append({"sample":sample,"png":files.get("png",""),"svg":files.get("svg","")})
    pd.DataFrame(plot_index).to_csv(outdir/"plot_index.tsv", sep="\t", index=False)

    summary_rows=[]; stats_rows=[]
    calls_by_sample={s:g.copy() for s,g in calls_df.groupby("sample")} if not calls_df.empty else {}
    overlap_by_sample={s:g.copy() for s,g in phase2.stx_phage_overlap.groupby("sample")} if not phase2.stx_phage_overlap.empty else {}
    context_by_sample={s:g.copy() for s,g in phase2.phage_insertion_context.groupby("sample")} if not phase2.phage_insertion_context.empty else {}
    variant_by_sample={s:g.copy() for s,g in phase3.alignment_overview.groupby("sample")} if not phase3.alignment_overview.empty else {}

    for query in typing_results.queries:
        grp=calls_by_sample.get(query.sample, pd.DataFrame())
        og=overlap_by_sample.get(query.sample, pd.DataFrame())
        cg=context_by_sample.get(query.sample, pd.DataFrame())
        vg=variant_by_sample.get(query.sample, pd.DataFrame())
        subtypes=sorted({str(x) for x in grp["assigned_subtype"].tolist() if str(x)}) if not grp.empty else []
        exact_count=int(grp["exact_match"].sum()) if not grp.empty else 0
        variant_count=int((grp["variant_status"]=="subtype_like_variant").sum()) if not grp.empty else 0
        novel_count=int((grp["variant_status"]=="unresolved_or_novel").sum()) if not grp.empty else 0
        phage_ids=sorted({str(x) for x in og["phage_region_id"].tolist()}) if not og.empty else []
        named_sites=sorted({str(x) for x in cg["known_insertion_site_name"].tolist() if str(x)}) if not cg.empty else []
        trna_feats=sorted({str(x) for x in list(cg["nearest_trna_left"].tolist())+list(cg["nearest_trna_right"].tolist()) if str(x)}) if not cg.empty else []
        plot_generated=query.sample in plot_files
        summary_rows.append({"sample":query.sample,"stx_present":not grp.empty,"stx_family_calls":";".join(sorted({str(x) for x in grp["stx_family"].tolist()})) if not grp.empty else "","stx_subtype_calls":";".join(subtypes),"num_stx_copies":int(len(grp)),"exact_match_count":exact_count,"variant_match_count":variant_count,"novel_candidate_count":novel_count,"num_stx_phages":int(len(set(phage_ids))),"stx_in_phage":not og.empty,"phage_region_ids":";".join(phage_ids),"named_insertion_sites":";".join(named_sites),"nearest_trna_features":";".join(trna_feats),"variant_tree_generated":not vg.empty,"closed_genome_plot_generated":plot_generated})
        stats_rows.append({"sample":query.sample,"assembly_type":query.assembly_type,"num_contigs":query.num_contigs,"genome_size":query.genome_size,"stx_present":not grp.empty,"stx_subtype_calls":";".join(subtypes),"num_stx_copies":int(len(grp)),"num_stx_phages":int(len(set(phage_ids))),"stx_in_phage":not og.empty,"insertion_site_classes":";".join(sorted({str(x) for x in cg["insertion_site_class"].tolist()})) if not cg.empty else "","named_insertion_sites":";".join(named_sites),"variant_tree_generated":not vg.empty,"closed_genome_plot_generated":plot_generated,"notes":"no_stx_hit" if grp.empty else ""})

    pd.DataFrame(summary_rows).to_csv(outdir/"stx_summary.tsv", sep="\t", index=False)
    pd.DataFrame(stats_rows).to_csv(outdir/"STX_stats.tsv", sep="\t", index=False)
    md=pd.DataFrame([{"field":"tool","value":"STXit"},{"field":"version","value":version},{"field":"timestamp","value":datetime.utcnow().isoformat()+"Z"},{"field":"platform","value":platform.platform()},{"field":"db_fasta","value":str(db.fasta_path)},{"field":"db_reference_count","value":len(db.references)},{"field":"phastest_job_count","value":len(phastest_jobs)},{"field":"trnascan_run_count","value":len(trna_runs)},{"field":"cli_args","value":json.dumps(cli_args, sort_keys=True)}])
    md.to_csv(outdir/"run_metadata.tsv", sep="\t", index=False)
