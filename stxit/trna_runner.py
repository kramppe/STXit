from __future__ import annotations
from pathlib import Path
from typing import List, Tuple
import shutil, subprocess
import pandas as pd
from .io_utils import QueryGenome
from .trna_scan import load_trna_sites

def _ensure_binary(binary: str) -> str:
    found=shutil.which(binary)
    if not found:
        raise RuntimeError(f"tRNAscan-SE executable not found: {binary}")
    return found

def _run_single(binary: str, query: QueryGenome, outdir: Path):
    txt_out=outdir/f"{query.sample}__trnascan.out.tsv"
    stats_out=outdir/f"{query.sample}__trnascan.stats.txt"
    subprocess.run([binary, "-o", str(txt_out), "-m", str(stats_out), str(query.path)], check=True)
    return txt_out, stats_out

def run_trnascan_for_queries(queries: List[QueryGenome], outdir: Path, binary: str = "tRNAscan-SE") -> Tuple[pd.DataFrame, list[dict]]:
    binary=_ensure_binary(binary)
    trna_dir=outdir/"trnascan_auto"; trna_dir.mkdir(parents=True, exist_ok=True)
    dfs=[]; run_rows=[]
    for query in queries:
        txt_out, stats_out = _run_single(binary, query, trna_dir)
        df=load_trna_sites(txt_out)
        df["sample"]=query.sample
        dfs.append(df)
        run_rows.append({"sample":query.sample, "trnascan_output":txt_out.name, "trnascan_stats":stats_out.name})
    trna_sites=pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=["sample","contig","start","end","trna_type","anticodon","strand"])
    return trna_sites, run_rows
