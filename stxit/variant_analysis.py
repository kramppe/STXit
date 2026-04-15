from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import shutil, subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .stx_typing import TypingResults
from .stx_db import StxDatabase

AA_CHARGE={"A":"neutral","V":"neutral","L":"neutral","I":"neutral","M":"neutral","F":"neutral","W":"neutral","P":"neutral","G":"neutral","S":"neutral","T":"neutral","Y":"neutral","C":"neutral","N":"neutral","Q":"neutral","D":"negative","E":"negative","K":"positive","R":"positive","H":"positive","*":"stop","X":"unknown"}

@dataclass
class Phase3VariantOutputs:
    variant_differences: pd.DataFrame
    alignment_overview: pd.DataFrame
    tree_fasta_paths: Dict[str, str]
    tree_newick_paths: Dict[str, str]

def _load_db_sequences(db: StxDatabase) -> Dict[str, str]:
    return {rec.id.split("|")[0]: str(rec.seq).upper() for rec in SeqIO.parse(str(db.fasta_path), "fasta")}

def _load_query_contigs(path: Path, file_type: str) -> Dict[str, str]:
    fmt="genbank" if file_type=="genbank" else "fasta"
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(path), fmt)}

def _slice_query_sequence(seqs, contig, start, end, strand):
    seq=seqs.get(contig,"")
    if not seq: return ""
    subseq=seq[start-1:end]
    if strand=="-": subseq=str(Seq(subseq).reverse_complement())
    return subseq.upper()

def _translate_safe(nt: str) -> str:
    nt=nt.upper().replace("-","N")
    if len(nt)!=3: return "X"
    try: return str(Seq(nt).translate(table=11))
    except Exception: return "X"

def _codon_aware_diffs(ref_seq: str, query_seq: str):
    diffs=[]; max_len=min(len(ref_seq), len(query_seq))
    for i in range(max_len):
        if ref_seq[i] != query_seq[i]:
            codon_index=i//3; codon_start=codon_index*3
            ref_codon=ref_seq[codon_start:codon_start+3] if codon_start+3 <= len(ref_seq) else ""
            query_codon=query_seq[codon_start:codon_start+3] if codon_start+3 <= len(query_seq) else ""
            ref_aa=_translate_safe(ref_codon) if ref_codon else ""
            query_aa=_translate_safe(query_codon) if query_codon else ""
            syn_class="synonymous" if ref_aa==query_aa else ("nonsense" if query_aa=="*" else "nonsynonymous")
            diffs.append({
                "ref_pos":i+1,"codon_position":(i%3)+1,"ref_base":ref_seq[i],"query_base":query_seq[i],"base_change":f"{ref_seq[i]}>{query_seq[i]}",
                "ref_codon":ref_codon,"query_codon":query_codon,"ref_aa":ref_aa,"query_aa":query_aa,
                "ref_aa_charge":AA_CHARGE.get(ref_aa,"unknown"),"alt_aa_charge":AA_CHARGE.get(query_aa,"unknown"),
                "charge_change":"same" if AA_CHARGE.get(ref_aa,"")==AA_CHARGE.get(query_aa,"") else "changed",
                "syn_nonsyn_nonsense":syn_class,"stop_introduced":query_aa=="*" and ref_aa!="*"
            })
    if len(ref_seq)!=len(query_seq):
        diffs.append({"ref_pos":max_len+1,"codon_position":"","ref_base":f"len={len(ref_seq)}","query_base":f"len={len(query_seq)}","base_change":"length_difference","ref_codon":"","query_codon":"","ref_aa":"","query_aa":"","ref_aa_charge":"","alt_aa_charge":"","charge_change":"","syn_nonsyn_nonsense":"indel_or_length_difference","stop_introduced":False})
    return diffs

def _write_tree_and_alignment(sample, subtype, query_seq, ref_id, ref_seq, outdir: Path, threads: int, build_trees: bool):
    fasta_path=outdir/f"{sample}__{subtype}__stx_variant_tree.fasta"
    tree_path=outdir/f"{sample}__{subtype}__stx_variant_tree.nwk"
    SeqIO.write([SeqRecord(Seq(query_seq), id=f"{sample}|query", description=""), SeqRecord(Seq(ref_seq), id=f"{ref_id}|reference", description="")], str(fasta_path), "fasta")
    if build_trees:
        iqtree=shutil.which("iqtree2") or shutil.which("iqtree")
        fasttree=shutil.which("FastTree") or shutil.which("fasttree")
        if iqtree:
            prefix=fasta_path.with_suffix("")
            subprocess.run([iqtree, "-s", str(fasta_path), "-nt", str(threads), "-quiet", "-pre", str(prefix)], check=False)
            produced=prefix.with_suffix(".treefile")
            if produced.exists():
                produced.replace(tree_path)
        elif fasttree:
            with open(tree_path,"w",encoding="utf-8") as handle:
                subprocess.run([fasttree, "-nt", str(fasta_path)], stdout=handle, check=False)
        else:
            tree_path.write_text(f"({sample}|query,{ref_id}|reference);\\n", encoding="utf-8")
    else:
        tree_path.write_text(f"({sample}|query,{ref_id}|reference);\\n", encoding="utf-8")
    return fasta_path.name, tree_path.name

def build_variant_outputs(typing_results: TypingResults, db: StxDatabase, outdir: Path, threads: int = 1, build_trees: bool = False) -> Phase3VariantOutputs:
    db_seqs=_load_db_sequences(db)
    query_map={q.sample:q for q in typing_results.queries}
    diff_rows=[]; overview_rows=[]; tree_fastas={}; tree_newicks={}
    for call in typing_results.calls:
        if call.exact_match: continue
        query=query_map.get(call.sample)
        if query is None: continue
        query_seqs=_load_query_contigs(query.path, query.file_type)
        query_seq=_slice_query_sequence(query_seqs, call.contig, call.start, call.end, call.strand)
        ref_seq=db_seqs.get(call.nearest_reference_id,"")
        if not query_seq or not ref_seq: continue
        diffs=_codon_aware_diffs(ref_seq, query_seq)
        for d in diffs:
            diff_rows.append({"sample":call.sample,"contig":call.contig,"assigned_subtype":call.assigned_subtype,"nearest_reference_id":call.nearest_reference_id,"nearest_reference_accession":call.nearest_reference_accession,"variant_status":call.variant_status, **d})
        overview_rows.append({"sample":call.sample,"contig":call.contig,"assigned_subtype":call.assigned_subtype,"nearest_reference_id":call.nearest_reference_id,"nearest_reference_accession":call.nearest_reference_accession,"pct_identity":round(call.pct_identity,3),"pct_coverage":round(call.pct_coverage,3),"query_length":len(query_seq),"reference_length":len(ref_seq),"mismatch_count":len([d for d in diffs if d["base_change"]!="length_difference"]),"nonsynonymous_count":len([d for d in diffs if d.get("syn_nonsyn_nonsense")=="nonsynonymous"]),"nonsense_count":len([d for d in diffs if d.get("syn_nonsyn_nonsense")=="nonsense"]),"length_difference":len(query_seq)-len(ref_seq),"exact_match":call.exact_match,"variant_status":call.variant_status})
        fasta_name, tree_name=_write_tree_and_alignment(call.sample, call.assigned_subtype, query_seq, call.nearest_reference_id, ref_seq, outdir, threads, build_trees)
        key=f"{call.sample}::{call.assigned_subtype}"; tree_fastas[key]=fasta_name; tree_newicks[key]=tree_name
    return Phase3VariantOutputs(pd.DataFrame(diff_rows), pd.DataFrame(overview_rows), tree_fastas, tree_newicks)
