# Changelog

All notable changes to STXit will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.4] - 2026-04-16

### Fixed
- **Critical**: `databases/stx_references.fasta` contained protein (amino acid) sequences
  but the pipeline runs `blastn` (nucleotide BLAST), causing zero STX loci to be detected
  on any genome with the default database. Replaced with 145 real nucleotide sequences
  covering all 11 canonical subtypes (stx1a, stx1c, stx1d, stx2a–stx2h) plus stx2i and
  stx2k. Sequences sourced from `config/comprehensive_stx_typing.fasta` (commit db5b496^).
- `stxit/stx_db.py`: `List` was missing from the `from typing import ...` line, causing
  a `NameError` when `validate_database()` was called.
- `stxit/stx_db.py`: `_parse_header` now correctly extracts the subtype from plain-name
  FASTA headers (e.g. `>stx2c`) in addition to pipe-separated format.

### Verified
- E. coli O157:H7 EC4115 (NC_011353.1, 5.57 Mb): 242 probe-level BLAST hits collapse
  to 2 clean loci — **stx2c** at 2,696,202 bp (99.9% identity) and **stx2a** at
  3,271,007 bp (99.9% identity). Both within intact PHASTEST-confirmed prophage regions
  (R12 sbcB locus, R15 argW locus).

## [1.0.3] - 2026-04-16

### Added
- Complete STX detection pipeline with BLAST-based subtyping
- Support for all 11 established STX subtypes (stx1a-stx1d, stx2a-stx2h)
- Comprehensive STX reference database with A/B subunit protein accessions
- PHASTEST integration for prophage detection and context analysis
- tRNAscan-SE integration for insertion site analysis
- Codon-aware variant analysis with synonymous/non-synonymous classification
- Indel detection with IS element identification via TnCentral BLAST
- Genome visualization with linear and circular maps
- Multiple output formats (TSV, JSON, PNG, SVG)
- Docker container support with full environment
- VS Code Dev Container configuration
- Comprehensive test suite and CI/CD pipeline
- Auto-detection of database paths (no manual configuration required)

### Database
- STX reference database v1.0.3 with 22 protein sequences
- Metadata mapping for all subtypes with strain references
- Integration with TnCentral for transposon element detection
- Automated database version tracking and validation

### Installation
- Three installation options: Conda, Dev Container, Docker
- One-command setup with `bash install.sh`
- All bioinformatics dependencies managed automatically
- Cross-platform support (Linux, macOS, Windows via containers)

### Documentation
- Comprehensive README with pipeline schema
- Test data examples with expected outputs
- Complete command-line reference
- Installation and troubleshooting guides

## [Unreleased]

### Planned
- Additional STX subtype variants as they are characterized
- Performance optimizations for large genome datasets
- Interactive web dashboard for result visualization
- Integration with additional prophage detection tools
- Phylogenetic analysis of detected STX sequences

---

For more details, see the [GitHub releases page](https://github.com/kramppe/STXit/releases).
