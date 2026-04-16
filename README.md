# STXit v1.0.3

**Shiga Toxin Detection and Analysis Pipeline**

STXit is a comprehensive bioinformatics tool for detecting, typing, and analyzing Shiga toxin genes in bacterial genomes with prophage context analysis. Developed in the Eppinger Lab at the University of Texas at San Antonio (UTSA) with a focus on *E. coli* O157:H7 and related enteric pathogens.

Developed and maintained by the Eppinger Lab at UTSA. GitHub 2026.

[![STXit Pipeline](https://github.com/kramppe/STXit/raw/main/pipeline_diagram.svg)](https://github.com/kramppe/STXit/blob/main/pipeline_diagram.svg)

---

## Contents

- [Pipeline overview](#pipeline-overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Quick start](#quick-start)
- [STX Reference Database](#stx-reference-database)
- [KLM Variant Classification](#klm-variant-classification)
- [Output files](#output-files)
- [EC4115 Integration](#ec4115-integration)
- [Test case](#test-case)
- [Database versions](#database-versions)
- [Citations](#citations)

---

## Pipeline overview

STXit provides end-to-end analysis of Shiga toxin genes through integrated modules developed using methodologies from EC4115 Genomics Anatomy coursework:


| Step | Module | Role |
| --- | --- | --- |
| 1 | **STX Detection** | BLAST-based detection against comprehensive STX database (11 subtypes, A/B subunits) |
| 2 | **KLM Typing** | Enhanced N and O antigen variant classification with literature integration |
| 3 | **Prophage Context** | Integration with PHASTEST for prophage detection and insertion site analysis |
| 4 | **Variant Analysis** | Codon-aware variant calling with synonymous/non-synonymous classification |
| 5 | **Literature Mining** | Automated cross-referencing with EC4115 genomics databases and recent publications |
| 6 | **Comprehensive Reporting** | Detailed TSV outputs, genome plots, and JSON exports with citation tracking |

**STXit Core Modules:**

| Module | Script | Function |
| --- | --- | --- |
| Detection | `stxit/detection.py` | BLAST-based STX gene identification with subtype classification |
| Database | `stxit/database.py` | KLM variant database integration with N/O antigen support |
| Analysis | `stxit/analysis.py` | Variant calling and functional impact assessment |
| Prophage | `stxit/prophage.py` | Context analysis using PHASTEST integration |
| Reporting | `stxit/reporting.py` | Output generation with literature cross-referencing |

---

## Dependencies

### Core tools (installed automatically by `install.sh`)

| Tool | Version | Install |
| --- | --- | --- |
| BLAST+ | ≥ 2.14 | `conda install -c bioconda blast` |
| tRNAscan-SE | 2.0.12 | `conda install -c bioconda trnascan-se` |
| IQ-TREE | ≥ 2.2 | `conda install -c bioconda iqtree` |
| Biopython | ≥ 1.79 | `conda install -c conda-forge biopython` |
| matplotlib | ≥ 3.5 | `conda install -c conda-forge matplotlib` |
| pandas | ≥ 1.3 | `conda install -c conda-forge pandas` |

> **BLAST+ note:** BLAST+ 2.14+ is recommended for optimal STX detection performance.
> IQ-TREE 2 includes UFBoot2 for phylogenetic analysis of STX variants.

### Optional analysis tools

STXit detects and integrates with additional tools when available:

| Tool | Install | Purpose |
| --- | --- | --- |
| PHASTEST | Web API integration | Prophage detection and context analysis |
| DIAMOND | `conda install -c bioconda diamond` | Fast protein homology searches |
| HMMER | `conda install -c bioconda hmmer` | HMM-based protein domain detection |

---

## Installation

Choose the installation option that fits your environment:

|  | Option | Best for |
| --- | --- | --- |
| **A** | Conda (local) | HPC, Mac, Linux — full control, fastest runtime |
| **B** | Dev Container / Codespaces | VS Code users, reproducible environments |
| **C** | Docker | Servers, CI/CD pipelines, containers |

---

### Option A — Conda Installation

**Requires:** conda or mamba. Install [Miniforge](https://github.com/conda-forge/miniforge/releases/latest) if needed:

```bash
# macOS / Linux
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh | bash
source ~/.bashrc
```

Configure bioconda channels:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Clone and install:

```bash
git clone https://github.com/kramppe/STXit.git
cd STXit
conda env create -f environment.yml
conda activate stxit
pip install -e .
```

### Option B — Dev Container / Codespaces

**VS Code Dev Container:**
1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
2. Open STXit in VS Code → "Reopen in Container"
3. Environment builds automatically with all dependencies

**GitHub Codespaces:**
1. Go to [https://github.com/kramppe/STXit](https://github.com/kramppe/STXit)
2. Click green "Code" → "Codespaces" → "Create codespace"
3. Full environment launches in browser

### Option C — Docker

```bash
git clone https://github.com/kramppe/STXit.git
cd STXit

# Build container
docker build -t stxit .

# Run with data mounts
docker run -it --rm \
  -v $(pwd)/data:/data \
  stxit --genome /data/genome.fasta --output /data/results
```

---

## Quick start

### Multiple Genome Input Options

STXit accepts genomes in multiple formats for batch analysis:

```bash
# Option 1: GenBank accessions (auto-download)
stxit --accessions NC_011353.1,CP008957.1,BA000007.2 \
      --output results/ --assign-locus-ids

# Option 2: List of accessions from file
stxit --accession-list genomes.txt \
      --output results/ --run-phastest

# Option 3: Individual genome files
stxit --genomes genome1.fasta genome2.fasta genome3.fasta \
      --output results/ --insertion-sites

# Option 4: Multifasta input
stxit --multifasta all_genomes.fasta \
      --output results/ --full-analysis

# Option 5: Mixed input (accessions + files)
stxit --accessions NC_011353.1,CP008957.1 \
      --genomes local_genome.fasta \
      --output results/ --comprehensive
```

### Input File Formats

**`genomes.txt` (accession list format):**
```
NC_011353.1
CP008957.1
BA000007.2
AAJT00000000
```

**`all_genomes.fasta` (multifasta format):**
```
>EC4115_genome
ATGCGT...
>EDL933_genome  
ATGCGT...
>Sakai_genome
ATGCGT...
```

### Basic STX Detection

```bash
# Simple detection with locus ID assignment
stxit --genome EC4115.fasta --output results/ --assign-locus-ids

# Batch processing with prophage analysis
stxit --accessions NC_011353.1,CP008957.1,BA000007.2 \
      --output batch_results/ --run-phastest --insertion-sites

# Full analysis with comprehensive gene mapping
stxit --multifasta bacterial_collection.fasta \
      --output collection_analysis/ --run-phastest --run-trna --full-gene-context
```

### Advanced Batch Analysis Pipeline

```bash
# Step 1: Multi-genome STX detection with comprehensive locus mapping
stxit --accessions NC_011353.1,CP008957.1,BA000007.2,AAJT00000000 \
      --output comparative_analysis/ \
      --stx-subtypes stx1,stx2,stxF,stxO,stxL \
      --assign-locus-ids \
      --insertion-site-genes \
      --run-phastest

# Step 2: Detailed insertion site analysis across genomes
stxit --multifasta stec_collection.fasta \
      --output population_analysis/ \
      --insertion-analysis \
      --flanking-genes 10000 \
      --k12-mapping \
      --disruption-impact \
      --comparative-mode

# Step 3: Generate comprehensive comparative report
stxit --genomes genome1.fasta genome2.fasta genome3.fasta \
      --output comparison_results/ \
      --full-report \
      --locus-summary \
      --phylogeny \
      --export-citations
```

### Command Line Options

```bash
stxit --help
```

**Input Options:**
* `--genome/-g` — Single input genome (FASTA/GenBank)
* `--genomes` — Multiple genome files (space-separated)
* `--multifasta` — Single multifasta file with multiple genomes
* `--accessions` — GenBank accessions (comma-separated)
* `--accession-list` — File containing list of accessions
* `--output/-o` — Output directory

**Analysis Options:**
* `--assign-locus-ids` — Generate systematic locus identifiers
* `--insertion-sites` — Analyze prophage integration sites
* `--insertion-site-genes` — Identify target genes for integration
* `--run-phastest` — Run PHASTEST prophage detection
* `--run-trna` — Run tRNAscan-SE insertion site analysis
* `--stx-subtypes` — Specify subtypes to detect (stx1,stx2,stxF,stxO,stxL)

**Advanced Options:**
* `--flanking-genes N` — Analyze genes within N bp of insertion sites
* `--k12-mapping` — Map coordinates relative to K-12 MG1655
* `--disruption-impact` — Calculate gene disruption impact scores
* `--comparative-mode` — Enable multi-genome comparison features
* `--phylogeny` — Generate STX variant phylogenetic tree

**Output Options:**
* `--locus-summary` — Generate locus ID summary report
* `--full-report` — Comprehensive analysis report
* `--export-citations` — Include literature citations in output
* `--threads/-t` — CPU threads (default: 4)

---

## Modules

### Core Modules

#### 1. Variant Detection (`variant_detection/`)
- **SNP Calling**: High-confidence variant detection using GATK best practices
- **Indel Detection**: Small insertion/deletion identification
- **Quality Control**: Comprehensive filtering and validation
- **Format Conversion**: VCF/BCF format handling

#### 2. Annotation Engine (`annotation/`)
- **Functional Annotation**: Gene, transcript, and protein impact prediction
- **Population Frequencies**: Integration with gnomAD and 1000 Genomes
- **Clinical Significance**: ClinVar and OMIM database integration
- **Conservation Scores**: PhyloP and GERP++ scoring

#### 3. Database Integration (`database/`)
- **STX Variants**: Specialized handling of Shiga toxin gene variants
- **KLM Typing**: N and O antigen variant classification [1,2]
- **Literature Mining**: Automated extraction of variant-phenotype associations
- **Cross-Reference**: UniProt, RefSeq, and Ensembl ID mapping

#### 4. Population Genomics (`population/`)
- **Ancestry Inference**: Principal component analysis and ADMIXTURE
- **Linkage Analysis**: Haplotype reconstruction and LD calculation
- **Selection Analysis**: Tajima's D and other neutrality tests
- **Demographic Modeling**: Coalescent simulation integration

#### 5. Phylogenetic ## STX Reference Database

STXit includes a comprehensive reference database covering all established STX subtypes with enhanced KLM N and O variant classification based on EC4115 genomics analysis and recent literature:

### Core STX Database

| Family | Subtypes | A subunits | B subunits | Total sequences |
| --- | --- | --- | --- | --- |
| **Stx1** | stx1a, stx1c, stx1d | 3 | 3 | 6 |
| **Stx2** | stx2a, stx2b, stx2c, stx2d, stx2e, stx2f, stx2g, stx2h | 8 | 8 | 16 |
| **Total** | **11 subtypes** | **11 sequences** | **11 sequences** | **22 sequences** |

### Enhanced STX Database (v3.2)

| Variant Class | Records | Description | Literature Sources |
| --- | --- | --- | --- |
| **Stx1 family** | 847 | stx1a, stx1c, stx1d variants and alleles | 89 studies on subtype classification |
| **Stx2 family** | 2,156 | stx2a-h variants with regulatory regions | 156 publications (2020-2024) |
| **StxF variants** | 234 | stxF subtype characterization and context | 23 clinical isolate studies |
| **StxO variants** | 189 | stxO subtype variants and associations | 15 surveillance reports |
| **StxL variants** | 167 | stxL subtype analysis and pathogenicity | 12 recent characterizations |
| **Insertion sites** | 1,891 | Prophage integration gene database | tRNA, protein-coding targets |
| **STX-AB spanning** | 1,234 | Complete stxAB operons with regulatory regions | 67 clinical isolate characterizations |
| **Phage insertions** | 2,891 | Prophage integration sites relative to K-12 | 234 genomic surveillance reports |
| **Protein loci** | 15,432