from __future__ import annotations
import argparse
from pathlib import Path

from . import __version__
from .io_utils import parse_query_inputs, load_backbone_annotation
from .stx_db import load_stx_database, load_known_sites
from .stx_typing import run_typing
from .phage_prediction import load_phastest_regions
from .phastest_client import run_phastest_for_queries
from .trna_scan import load_trna_sites
from .trna_runner import run_trnascan_for_queries
from .phage_insertion import build_phase2_context
from .variant_analysis import build_variant_outputs
from .plotting import build_closed_genome_plots
from .reporting import write_outputs

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="stxit", description="STXit merged implementation.")
    p.add_argument("--query", nargs="+", required=True, help="Query genome files. Supports FASTA/GenBank. Use path:SampleName to override sample label.")
    p.add_argument("--db-fasta", required=True, help="STX subtype reference FASTA.")
    p.add_argument("--db-metadata", default=None, help="Optional TSV metadata.")
    p.add_argument("--db-citations", default=None, help="Optional TSV citations.")
    p.add_argument("--phastest-tsv", default=None, help="Manual parsed PHASTEST region table.")
    p.add_argument("--run-phastest", action="store_true", help="Submit queries to PHASTEST API.")
    p.add_argument("--phastest-email", default="", help="Email address for PHASTEST.")
    p.add_argument("--phastest-input-mode", default="file", choices=["file", "accession"], help="PHASTEST submission mode.")
    p.add_argument("--phastest-contigs", action="store_true", help="Submit multi-contig FASTA to PHASTEST.")
    p.add_argument("--phastest-poll-seconds", type=int, default=30, help="PHASTEST polling interval.")
    p.add_argument("--trnascan-tsv", default=None, help="Manual parsed tRNAscan-SE table.")
    p.add_argument("--run-trnascan", action="store_true", help="Run tRNAscan-SE locally.")
    p.add_argument("--trnascan-binary", default="tRNAscan-SE", help="Path or name of tRNAscan-SE executable.")
    p.add_argument("--known-sites", default=None, help="Known insertion-site catalog TSV.")
    p.add_argument("--compare-backbone", default="default_k12", help="Backbone comparison label.")
    p.add_argument("--backbone-comparison-genbank", default=None, help="Optional backbone GenBank for annotation-aware comparison.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--threads", type=int, default=1, help="Threads for BLAST and trees.")
    p.add_argument("--min-identity", type=float, default=90.0, help="Minimum BLAST identity threshold.")
    p.add_argument("--min-coverage", type=float, default=60.0, help="Minimum BLAST query coverage threshold.")
    p.add_argument("--blast-task", default="blastn", choices=["blastn", "blastn-short"], help="blastn task.")
    p.add_argument("--keep-temp", action="store_true", help="Keep temporary files.")
    p.add_argument("--plot-closed-genome-maps", action="store_true", help="Generate closed-genome phage/STX plots.")
    p.add_argument("--build-variant-trees", action="store_true", help="Generate subtype-reference trees for non-exact matches.")
    p.add_argument("--version", action="version", version=f"STXit {__version__}")
    return p

def main() -> None:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    queries = parse_query_inputs(args.query)
    db = load_stx_database(
        fasta_path=Path(args.db_fasta),
        metadata_path=Path(args.db_metadata) if args.db_metadata else None,
        citations_path=Path(args.db_citations) if args.db_citations else None,
    )
    known_sites = load_known_sites(Path(args.known_sites)) if args.known_sites else None
    backbone = load_backbone_annotation(Path(args.backbone_comparison_genbank), label=args.compare_backbone) if args.backbone_comparison_genbank else None

    typing = run_typing(
        queries=queries, db=db, outdir=outdir, threads=args.threads,
        min_identity=args.min_identity, min_coverage=args.min_coverage,
        blast_task=args.blast_task, keep_temp=args.keep_temp
    )

    phage_regions = None
    phastest_jobs = []
    if args.phastest_tsv:
        phage_regions = load_phastest_regions(Path(args.phastest_tsv))
    elif args.run_phastest:
        phage_regions, phastest_jobs = run_phastest_for_queries(
            queries=queries, outdir=outdir, email=args.phastest_email,
            input_mode=args.phastest_input_mode, contigs=args.phastest_contigs,
            poll_seconds=args.phastest_poll_seconds,
        )

    trna_sites = None
    trna_runs = []
    if args.trnascan_tsv:
        trna_sites = load_trna_sites(Path(args.trnascan_tsv))
    elif args.run_trnascan:
        trna_sites, trna_runs = run_trnascan_for_queries(
            queries=queries, outdir=outdir, binary=args.trnascan_binary
        )

    phase2 = build_phase2_context(
        typing_results=typing, phage_regions=phage_regions, trna_sites=trna_sites,
        known_sites=known_sites, compare_backbone=args.compare_backbone, backbone=backbone
    )

    phase3 = build_variant_outputs(
        typing_results=typing, db=db, outdir=outdir, threads=args.threads,
        build_trees=args.build_variant_trees
    )

    plot_files = build_closed_genome_plots(
        queries=queries, phase2=phase2, outdir=outdir, enabled=args.plot_closed_genome_maps
    )

    write_outputs(
        typing_results=typing, db=db, phase2=phase2, phase3=phase3,
        plot_files=plot_files, phastest_jobs=phastest_jobs, trna_runs=trna_runs,
        outdir=outdir, cli_args=vars(args), version=__version__
    )

if __name__ == "__main__":
    main()
