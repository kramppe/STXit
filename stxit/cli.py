#!/usr/bin/env python3
"""
STXit Command Line Interface
Pipeline orchestration and argument parsing
"""

import argparse
import sys
import os
from pathlib import Path
from .stx_typing import STXTyper
from .stx_db import STXDatabase
from .phastest_client import PhastestClient
from .trna_runner import tRNARunner
from .variant_analysis import VariantAnalyzer
from .plotting import GenomePlotter
from .reporting import STXReporter
from .__init__ import __version__

def create_parser():
    """Create command line argument parser"""
    parser = argparse.ArgumentParser(
        prog='stxit',
        description='STXit - Shiga Toxin Detection and Analysis Pipeline v' + __version__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic STX detection
  stxit --genome genome.fasta --output results/

  # With prophage analysis
  stxit --genome genome.fasta --output results/ --run-phastest

  # Manual PHASTEST results
  stxit --genome genome.fasta --output results/ --phastest-results phastest_output.txt

  # Full analysis with tRNA scanning
  stxit --genome genome.fasta --output results/ --run-phastest --run-trna

  # Custom STX database
  stxit --genome genome.fasta --output results/ --stx-db custom_stx.fasta
        '''
    )

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '--genome', '-g',
        required=True,
        help='Input genome file (FASTA format)'
    )
    required.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory'
    )

    # Database arguments
    db_group = parser.add_argument_group('database options')
    db_group.add_argument(
        '--stx-db',
        default=None,
        help='STX reference database (default: bundled database)'
    )

    # Analysis options
    analysis_group = parser.add_argument_group('analysis options')
    analysis_group.add_argument(
        '--run-phastest',
        action='store_true',
        help='Run PHASTEST for prophage detection'
    )
    analysis_group.add_argument(
        '--phastest-results',
        help='Pre-computed PHASTEST results file'
    )
    analysis_group.add_argument(
        '--run-trna',
        action='store_true',
        help='Run tRNAscan-SE for insertion site analysis'
    )
    analysis_group.add_argument(
        '--trna-results',
        help='Pre-computed tRNAscan-SE results file'
    )

    # BLAST options
    blast_group = parser.add_argument_group('BLAST parameters')
    blast_group.add_argument(
        '--min-identity',
        type=float,
        default=80.0,
        help='Minimum BLAST identity threshold (default: 80.0)'
    )
    blast_group.add_argument(
        '--min-coverage',
        type=float,
        default=70.0,
        help='Minimum BLAST coverage threshold (default: 70.0)'
    )

    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument(
        '--prefix',
        default='stx',
        help='Output file prefix (default: stx)'
    )
    output_group.add_argument(
        '--no-plots',
        action='store_true',
        help='Skip genome plot generation'
    )

    # Misc options
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=4,
        help='Number of threads (default: 4)'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Quiet mode - minimal output'
    )
    parser.add_argument(
        '--version', '-v',
        action='version',
        version=f'STXit {__version__}'
    )

    return parser

def validate_args(args):
    """Validate command line arguments"""
    errors = []

    # Check input genome file
    if not os.path.exists(args.genome):
        errors.append(f"Genome file not found: {args.genome}")

    # Check PHASTEST inputs
    if args.run_phastest and args.phastest_results:
        errors.append("Cannot specify both --run-phastest and --phastest-results")

    # Check tRNA inputs
    if args.run_trna and args.trna_results:
        errors.append("Cannot specify both --run-trna and --trna-results")

    # Check external result files
    if args.phastest_results and not os.path.exists(args.phastest_results):
        errors.append(f"PHASTEST results file not found: {args.phastest_results}")

    if args.trna_results and not os.path.exists(args.trna_results):
        errors.append(f"tRNA results file not found: {args.trna_results}")

    # Validate thresholds
    if not 0 <= args.min_identity <= 100:
        errors.append("Min identity must be between 0-100")

    if not 0 <= args.min_coverage <= 100:
        errors.append("Min coverage must be between 0-100")

    if args.threads < 1:
        errors.append("Threads must be >= 1")

    return errors

def main():
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args()

    # Validate arguments
    errors = validate_args(args)
    if errors:
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.quiet:
        print(f"STXit v{__version__}")
        print(f"Genome: {args.genome}")
        print(f"Output: {args.output}")
        print("-" * 50)

    try:
        # Initialize components
        if not args.quiet:
            print("Loading STX database...")
        
        stx_db = STXDatabase(args.stx_db)
        typer = STXTyper(stx_db, args.min_identity, args.min_coverage)
        reporter = STXReporter(args.prefix)

        # Run STX detection
        if not args.quiet:
            print("Running STX detection...")
        
        stx_results = typer.detect_stx(args.genome, args.threads)

        if not stx_results:
            if not args.quiet:
                print("No STX genes detected.")
            # Still create empty output files
            reporter.write_results([], output_dir)
            return

        if not args.quiet:
            print(f"Detected {len(stx_results)} STX loci")

        # Prophage analysis
        prophage_data = None
        if args.run_phastest:
            if not args.quiet:
                print("Running PHASTEST...")
            client = PhastestClient()
            prophage_data = client.submit_genome(args.genome)
        elif args.phastest_results:
            if not args.quiet:
                print("Loading PHASTEST results...")
            from .phage_prediction import PhagePredictor
            predictor = PhagePredictor()
            prophage_data = predictor.load_phastest_results(args.phastest_results)

        # tRNA analysis
        trna_data = None
        if args.run_trna:
            if not args.quiet:
                print("Running tRNAscan-SE...")
            trna_runner = tRNARunner()
            trna_data = trna_runner.scan_genome(args.genome)
        elif args.trna_results:
            if not args.quiet:
                print("Loading tRNA results...")
            trna_runner = tRNARunner()
            trna_data = trna_runner.load_results(args.trna_results)

        # Variant analysis
        if not args.quiet:
            print("Analyzing variants...")
        analyzer = VariantAnalyzer()
        for result in stx_results:
            analyzer.analyze_variants(result, stx_db)

        # Generate plots
        if not args.no_plots:
            if not args.quiet:
                print("Generating genome plots...")
            plotter = GenomePlotter()
            plotter.plot_genome(args.genome, stx_results, 
                              prophage_data, trna_data, output_dir)

        # Write results
        if not args.quiet:
            print("Writing results...")
        reporter.write_results(stx_results, output_dir, 
                             prophage_data, trna_data)

        if not args.quiet:
            print(f"Analysis complete. Results in: {output_dir}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        if not args.quiet:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
