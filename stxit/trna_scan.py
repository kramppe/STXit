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

def load_trna_sites(path: Path) -> pd.DataFrame:
    df=pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols=list(df.columns)
    sample_col=_first_match(cols,["sample","genome","strain"])
    contig_col=_first_match(cols,["contig","scaffold","sequence","seqid"])
    start_col=_first_match(cols,["start","begin"])
    end_col=_first_match(cols,["end","stop"])
    type_col=_first_match(cols,["trna_type","isotype","type"])
    anti_col=_first_match(cols,["anticodon"])
    strand_col=_first_match(cols,["strand"])
    out=pd.DataFrame()
    out["sample"]=df[sample_col] if sample_col else ""
    out["contig"]=df[contig_col] if contig_col else ""
    out["start"]=pd.to_numeric(df[start_col], errors="coerce").fillna(0).astype(int) if start_col else 0
    out["end"]=pd.to_numeric(df[end_col], errors="coerce").fillna(0).astype(int) if end_col else 0
    out["trna_type"]=df[type_col] if type_col else ""
    out["anticodon"]=df[anti_col] if anti_col else ""
    out["strand"]=df[strand_col] if strand_col else ""
    return out
