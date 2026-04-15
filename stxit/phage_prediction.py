from __future__ import annotations
from pathlib import Path
import pandas as pd

def _first_match(cols, candidates):
    norm={c.lower().strip():c for c in cols}
    for cand in candidates:
        if cand.lower() in norm:
            return norm[cand.lower()]
    for c in cols:
        lc=c.lower().strip()
        for cand in candidates:
            if cand.lower() in lc:
                return c
    return None

def load_phastest_regions(path: Path) -> pd.DataFrame:
    df=pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols=list(df.columns)
    sample_col=_first_match(cols,["sample","genome","strain"])
    contig_col=_first_match(cols,["contig","scaffold","sequence","seqid"])
    region_col=_first_match(cols,["region","region_id","phage_region"])
    start_col=_first_match(cols,["start"])
    end_col=_first_match(cols,["end","stop"])
    comp_col=_first_match(cols,["completeness","status","class"])
    score_col=_first_match(cols,["score"])
    attl_col=_first_match(cols,["attl"])
    attr_col=_first_match(cols,["attr"])
    note_col=_first_match(cols,["note","notes"])
    out=pd.DataFrame()
    out["sample"]=df[sample_col] if sample_col else ""
    out["contig"]=df[contig_col] if contig_col else ""
    out["region_id"]=df[region_col] if region_col else [f"region_{i+1}" for i in range(len(df))]
    out["start"]=pd.to_numeric(df[start_col], errors="coerce").fillna(0).astype(int) if start_col else 0
    out["end"]=pd.to_numeric(df[end_col], errors="coerce").fillna(0).astype(int) if end_col else 0
    out["length"]=(out["end"]-out["start"]+1).clip(lower=0)
    out["completeness"]=df[comp_col] if comp_col else ""
    out["score"]=df[score_col] if score_col else ""
    out["attl"]=df[attl_col] if attl_col else ""
    out["attr"]=df[attr_col] if attr_col else ""
    out["note"]=df[note_col] if note_col else ""
    return out
