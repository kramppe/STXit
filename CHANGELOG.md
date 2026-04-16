# Changelog

All notable changes to STXit will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
