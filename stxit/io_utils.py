from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO

@dataclass
class FeatureRecord:
    start: int
    end: int
    strand: str
    feature_type: str
    gene: str
    locus_tag: str
    product: str

@dataclass
class QueryGenome:
    path: Path
    sample: str
    file_type: str
    assembly_type: str
    contig_lengths: Dict[str, int]
    contig_order: List[str]
    features: Dict[str, List[FeatureRecord]] = field(default_factory=dict)

    @property
    def num_contigs(self) -> int:
        return len(self.contig_order)

    @property
    def genome_size(self) -> int:
        return sum(self.contig_lengths.values())

@dataclass
class BackboneAnnotation:
    label: str
    path: Path
    features: List[FeatureRecord]
    gene_index: Dict[str, FeatureRecord]
    locus_index: Dict[str, FeatureRecord]

def _parse_query_spec(spec: str) -> Tuple[Path, Optional[str]]:
    if ":" in spec:
        path_str, label = spec.rsplit(":", 1)
        if path_str:
            return Path(path_str), label
    return Path(spec), None

def _detect_file_type(path: Path) -> str:
    return "genbank" if path.suffix.lower() in {".gb", ".gbk", ".genbank", ".gbff"} else "fasta"

def _parse_genbank(path: Path):
    contig_lengths, contig_order, features = {}, [], {}
    for record in SeqIO.parse(str(path), "genbank"):
        contig_id = record.id
        contig_order.append(contig_id)
        contig_lengths[contig_id] = len(record.seq)
        feats = []
        for feat in record.features:
            if feat.type not in {"CDS", "gene", "tRNA", "rRNA"}:
                continue
            q = feat.qualifiers
            feats.append(FeatureRecord(
                start=int(feat.location.start)+1,
                end=int(feat.location.end),
                strand="+" if feat.location.strand != -1 else "-",
                feature_type=feat.type,
                gene=";".join(q.get("gene", [""])),
                locus_tag=";".join(q.get("locus_tag", [""])),
                product=";".join(q.get("product", [""])),
            ))
        feats.sort(key=lambda x: (x.start, x.end))
        features[contig_id] = feats
    return contig_lengths, contig_order, features

def _parse_fasta(path: Path):
    contig_lengths, contig_order = {}, []
    for record in SeqIO.parse(str(path), "fasta"):
        contig_order.append(record.id)
        contig_lengths[record.id] = len(record.seq)
    return contig_lengths, contig_order

def parse_query_inputs(query_specs: List[str]) -> List[QueryGenome]:
    queries = []
    for spec in query_specs:
        path, label = _parse_query_spec(spec)
        if not path.exists():
            raise FileNotFoundError(f"Query file not found: {path}")
        file_type = _detect_file_type(path)
        sample = label or path.stem
        if file_type == "genbank":
            contig_lengths, contig_order, features = _parse_genbank(path)
            assembly_type = "closed" if len(contig_order) == 1 else "draft"
            queries.append(QueryGenome(path, sample, file_type, assembly_type, contig_lengths, contig_order, features))
        else:
            contig_lengths, contig_order = _parse_fasta(path)
            assembly_type = "closed" if len(contig_order) == 1 else "draft"
            queries.append(QueryGenome(path, sample, file_type, assembly_type, contig_lengths, contig_order, {}))
    return queries

def load_backbone_annotation(path: Path, label: str = "default_k12") -> BackboneAnnotation:
    _, _, features_by_contig = _parse_genbank(path)
    features = []
    for feats in features_by_contig.values():
        features.extend(feats)
    gene_index = {f.gene: f for f in features if f.gene}
    locus_index = {f.locus_tag: f for f in features if f.locus_tag}
    return BackboneAnnotation(label=label, path=path, features=features, gene_index=gene_index, locus_index=locus_index)

def contig_number(query: QueryGenome, contig_id: str) -> int:
    try:
        return query.contig_order.index(contig_id)+1
    except ValueError:
        return -1

def find_overlapping_feature(query: QueryGenome, contig_id: str, start: int, end: int):
    feats = query.features.get(contig_id, [])
    best, best_overlap = None, 0
    for feat in feats:
        ov = max(0, min(end, feat.end)-max(start, feat.start)+1)
        if ov > best_overlap:
            best, best_overlap = feat, ov
    return best

def find_flanking_features(query: QueryGenome, contig_id: str, start: int, end: int):
    feats = query.features.get(contig_id, [])
    left = right = None
    for feat in feats:
        if feat.end < start:
            left = feat
        elif feat.start > end:
            right = feat
            break
    return left, right
