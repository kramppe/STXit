from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List
import shutil, subprocess, tempfile
import pandas as pd
from .io_utils import QueryGenome, contig_number, find_overlapping_feature
from .stx_db import StxDatabase

@dataclass
class LocusCall:
    sample: str
    assembly_type: str
    contig: str
    contig_num: int
    start: int
    end: int
    strand: str
    stx_family: str
    assigned_subtype: str
    nearest_reference_id: str
    nearest_reference_accession: str
    nearest_reference_strain: str
    nearest_reference_citation: str
    pct_identity: float
    pct_coverage: float
    hit_length: int
    complete_partial: str
    exact_match: bool
    variant_status: str
    annotation_source: str
    query_feature_type: str
    query_gene: str
    query_locus_tag: str
    query_product: str

@dataclass
class TypingResults:
    calls: List[LocusCall]
    queries: List[QueryGenome]
    blast_tables: Dict[str, Path]

def _run(cmd, cwd=None):
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)

def _check_executable(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")

def _cluster_hits(df: pd.DataFrame) -> List[pd.DataFrame]:
    if df.empty:
        return []
    df = df.sort_values(["sseqid","strand_norm","sstart_norm","send_norm","bitscore"], ascending=[True,True,True,True,False]).copy()
    clusters=[]; current=[]; current_key=None; current_start=None; current_end=None
    for _, row in df.iterrows():
        key=(row["sseqid"], row["strand_norm"]); start=int(row["sstart_norm"]); end=int(row["send_norm"])
        if current_key is None:
            current=[row]; current_key=key; current_start=start; current_end=end; continue
        overlaps = key == current_key and start <= current_end + 200 and end >= current_start - 200
        if overlaps:
            current.append(row); current_start=min(current_start, start); current_end=max(current_end, end)
        else:
            clusters.append(pd.DataFrame(current)); current=[row]; current_key=key; current_start=start; current_end=end
    if current:
        clusters.append(pd.DataFrame(current))
    return clusters

def run_typing(queries: List[QueryGenome], db: StxDatabase, outdir: Path, threads: int, min_identity: float, min_coverage: float, blast_task: str, keep_temp: bool = False) -> TypingResults:
    _check_executable("makeblastdb"); _check_executable("blastn")
    calls=[]; blast_tables={}
    temp_root_ctx = tempfile.TemporaryDirectory(prefix="stxit_merged_", dir=str(outdir))
    temp_root = Path(temp_root_ctx.name)
    for query in queries:
        qdir=temp_root / query.sample; qdir.mkdir(parents=True, exist_ok=True)
        db_prefix=qdir / f"{query.sample}.db"
        _run(["makeblastdb","-in",str(query.path),"-dbtype","nucl","-out",str(db_prefix)])
        out=qdir / f"{query.sample}.blast.tsv"
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sstrand"
        _run(["blastn","-task",blast_task,"-query",str(db.fasta_path),"-db",str(db_prefix),"-out",str(out),"-outfmt",outfmt,"-num_threads",str(threads)])
        blast_tables[query.sample]=out
        if out.stat().st_size == 0:
            continue
        cols=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs","sstrand"]
        df=pd.read_csv(out, sep="\t", names=cols)
        if df.empty:
            continue
        df["pident"]=df["pident"].astype(float); df["qcovs"]=df["qcovs"].astype(float); df["bitscore"]=df["bitscore"].astype(float)
        df=df[(df["pident"]>=min_identity)&(df["qcovs"]>=min_coverage)].copy()
        if df.empty:
            continue
        df["strand_norm"]=df["sstrand"].replace({"plus":"+","minus":"-"})
        df["sstart_norm"]=df[["sstart","send"]].min(axis=1)
        df["send_norm"]=df[["sstart","send"]].max(axis=1)
        for cluster in _cluster_hits(df):
            best=cluster.sort_values(["bitscore","pident","qcovs","length"], ascending=[False,False,False,False]).iloc[0]
            ref=db.references.get(best["qseqid"])
            if ref is None:
                continue
            start=int(best["sstart_norm"]); end=int(best["send_norm"]); strand=str(best["strand_norm"])
            feat=find_overlapping_feature(query, str(best["sseqid"]), start, end)
            exact=float(best["pident"])>=100.0 and float(best["qcovs"])>=100.0
            complete_partial="complete" if float(best["qcovs"])>=95.0 else "partial"
            variant_status="exact_match" if exact else ("subtype_like_variant" if float(best["qcovs"])>=90.0 and float(best["pident"])>=95.0 else "unresolved_or_novel")
            calls.append(LocusCall(
                sample=query.sample, assembly_type=query.assembly_type, contig=str(best["sseqid"]), contig_num=contig_number(query, str(best["sseqid"])),
                start=start, end=end, strand=strand, stx_family=ref.stx_family, assigned_subtype=ref.stx_subtype,
                nearest_reference_id=ref.reference_id, nearest_reference_accession=ref.accession, nearest_reference_strain=ref.strain,
                nearest_reference_citation=ref.citation_key, pct_identity=float(best["pident"]), pct_coverage=float(best["qcovs"]),
                hit_length=int(best["length"]), complete_partial=complete_partial, exact_match=exact, variant_status=variant_status,
                annotation_source=query.file_type if feat else "", query_feature_type=feat.feature_type if feat else "",
                query_gene=feat.gene if feat else "", query_locus_tag=feat.locus_tag if feat else "", query_product=feat.product if feat else ""
            ))
    return TypingResults(calls=calls, queries=queries, blast_tables=blast_tables)
