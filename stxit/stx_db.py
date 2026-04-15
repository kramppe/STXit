from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional
import pandas as pd
from Bio import SeqIO

@dataclass
class StxReference:
    reference_id: str
    accession: str
    stx_family: str
    stx_subtype: str
    variant_name: str
    strain: str
    citation_key: str
    length: int

@dataclass
class KnownSite:
    site_name: str
    feature_type: str
    match_gene: str
    match_locus_tag: str
    backbone_accession: str
    backbone_label: str
    category: str
    note: str

@dataclass
class StxDatabase:
    fasta_path: Path
    references: Dict[str, StxReference]
    metadata: pd.DataFrame
    citations: Optional[pd.DataFrame] = None

def _parse_header(record_id: str, seq_len: int) -> StxReference:
    parts = record_id.split("|")
    parts += [""] * (7 - len(parts))
    return StxReference(parts[0] or record_id, parts[1] or "NA", parts[2] or "", parts[3] or "", parts[4] or "", parts[5] or "", parts[6] or "", seq_len)

def load_stx_database(fasta_path: Path, metadata_path: Optional[Path] = None, citations_path: Optional[Path] = None) -> StxDatabase:
    refs = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        parsed = _parse_header(rec.id, len(rec.seq))
        refs[parsed.reference_id] = parsed
    md_rows = [{"reference_id": r.reference_id, "accession": r.accession, "stx_family": r.stx_family, "stx_subtype": r.stx_subtype, "variant_name": r.variant_name, "strain": r.strain, "length": r.length, "citation_key": r.citation_key} for r in refs.values()]
    metadata = pd.DataFrame(md_rows)
    if metadata_path and metadata_path.exists():
        external = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
        metadata = external.merge(metadata, on="reference_id", how="outer", suffixes=("", "_fasta"))
        for col in ["accession", "stx_family", "stx_subtype", "variant_name", "strain", "length", "citation_key"]:
            if f"{col}_fasta" in metadata.columns:
                metadata[col] = metadata[col].replace("", pd.NA).fillna(metadata[f"{col}_fasta"])
        metadata = metadata[[c for c in metadata.columns if not c.endswith("_fasta")]]
        refs = {}
        for _, row in metadata.iterrows():
            rid = str(row.get("reference_id", ""))
            refs[rid] = StxReference(rid, str(row.get("accession", "")), str(row.get("stx_family", "")), str(row.get("stx_subtype", "")), str(row.get("variant_name", "")), str(row.get("strain", "")), str(row.get("citation_key", "")), int(float(row.get("length", 0) or 0)))
    citations = pd.read_csv(citations_path, sep="\t", dtype=str).fillna("") if citations_path and citations_path.exists() else None
    return StxDatabase(fasta_path=fasta_path, references=refs, metadata=metadata, citations=citations)

def load_known_sites(path: Path) -> list[KnownSite]:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return [KnownSite(str(r.get("site_name","")), str(r.get("feature_type","")), str(r.get("match_gene","")), str(r.get("match_locus_tag","")), str(r.get("backbone_accession","")), str(r.get("backbone_label","")), str(r.get("category","")), str(r.get("note",""))) for _, r in df.iterrows()]
