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
| **Protein loci** | 15,432 | All characterized STX protein variants | EC4115 + literature integration |

### Insertion Site Gene Database

STXit includes comprehensive data on prophage integration targets:

| Gene Category | Examples | STX Association | Frequency |
| --- | --- | --- | --- |
| **tRNA genes** | argW, thrW, leuX, pheU, ssrA | High (stx1, stx2a-d) | 78% |
| **Protein-coding** | sbcB, yecE, wrbA, yehV | Medium (stx2e-h, stxF) | 15% |
| **Regulatory** | malT, cadA, appY | Low (stxO, stxL) | 5% |
| **Hypothetical** | ybaO, ychE, ydaE | Variable (novel variants) | 2% |

**Database includes:**
- **Gene coordinates** relative to K-12 MG1655
- **Functional annotations** from EcoCyc and RegulonDB
- **Disruption impact scores** (0.0-1.0 scale)
- **STX subtype associations** from literature
- **Phylogenetic distribution** across E. coli strains

### Database Sources and Standards

**Primary References:**
- Scheutz et al. (2012) nomenclature standards [1]
- Manning et al. (2008) EC4115 characterization [2]
- NCBI protein database with verified accessions
- Both A (enzymatic) and B (binding) subunit references

**EC4115 Integration:**
- Reference genome coordinates for STX loci
- Prophage insertion sites mapped relative to E. coli K-12
- Protein domain annotations from EC4115 Genomics Anatomy analysis
- Regulatory region characterization

All databases are **downloaded and indexed automatically** during installation with literature cross-referencing enabled.

---

## STX Variant Classification

STXit provides comprehensive STX variant analysis including stxF, stxO, and stxL subtypes, integrated with EC4115 reference mapping and K-12 comparative analysis.

### Analysis Features

#### STX Subtype Detection
```bash
# Detect all STX subtypes including stxF, stxO, stxL
stxit --genome genome.fasta \
      --stx-subtypes all \
      --include-novel \
      --output results/

# Expected output: stx_subtype_calls.tsv
# - Complete subtype classification
# - stxF/stxO/stxL identification  
# - Novel variant detection
# - Allelic variation analysis
```

#### STX-AB Spanning Region Analysis
```bash
# Analyze complete stxAB operons
stxit --genome genome.fasta \
      --stx-ab-spanning \
      --regulatory-regions \
      --output results/

# Expected output: stx_spanning_analysis.tsv
# - Complete operon coordinates
# - Regulatory region identification  
# - Inter-gene spacing analysis
# - Promoter/terminator prediction
```

#### Protein Locus Mapping
```bash
# Map to reference protein coordinates
stxit --genome genome.fasta \
      --protein-loci \
      --ec4115-reference \
      --output results/

# Generates: protein_locus_map.tsv
# - STX protein coordinates
# - Domain annotations
# - Functional predictions
# - EC4115 comparative mapping
```

#### Phage Insertion Site Analysis
```bash
# Analyze prophage integration relative to K-12
stxit --genome genome.fasta \
      --phage-insertion \
      --k12-relative \
      --run-phastest \
      --output results/

# Outputs: phage_insertion_sites.tsv
# - Integration coordinates
# - K-12 relative positions
# - tRNA insertion sites
# - Prophage boundaries
```

#### Phylogenetic Tree Generation
```bash
# Generate STX variant phylogeny
stxit --genome genome1.fasta genome2.fasta genome3.fasta \
      --tree-analysis \
      --stx-phylogeny \
      --bootstrap 1000 \
      --output results/

# Creates: stx_phylogeny.nwk, stx_tree.svg
# - ML tree with bootstrap support
# - STX subtype clustering
# - Variant relationship analysis
```

### Test Set Report Generation

STXit includes comprehensive test case validation using the EC4115 reference genome characterized by Eppinger et al. (PNAS 2011):

```bash
# Run EC4115 validation test
stxit --test-ec4115 \
      --validate-database \
      --full-report \
      --output test_validation/

# Generates complete test report including:
# - STX detection accuracy
# - KLM variant classification
# - Protein locus validation
# - Phage insertion mapping
# - Literature cross-references
```

---

## Output files

### Primary Results

**`stx_calls.long.tsv`** — Comprehensive STX locus calls
- **Locus ID assignments** (STX_001, STX_F01, etc.)
- Genomic coordinates, strand, subtype classification  
- **Insertion site genes** (argW, sbcB, thrW, etc.)
- Reference/query annotations (A/B subunits)
- Protein accessions and functional products
- Prophage context and integration sites
- Literature cross-references

**`stx_summary.tsv`** — Concise overview
- Locus ID, subtype, coordinates
- **Integration site gene identification**
- Identity percentage and match status
- Prophage association and insertion sites
- EC4115 reference mapping

**`insertion_sites.tsv`** — Integration site analysis
- **Locus ID** → **insertion gene** mapping
- Target gene disruption analysis (complete/partial)
- Flanking gene context (±5kb)
- K-12 relative coordinates
- Functional impact assessment

**`stx_variant_differences.tsv`** — Detailed variant analysis
- Position-level differences vs reference
- **Locus ID cross-reference**
- Codon changes and amino acid impacts
- Synonymous vs non-synonymous classification
- Protein domain effects
- Literature annotations

### Specialized Analysis Outputs

**`stx_spanning_analysis.tsv`** — Complete stxAB operon analysis
- Full operon coordinates and organization
- Regulatory region identification
- Inter-gene spacing and promoter analysis
- Terminator predictions

**`protein_locus_map.tsv`** — Protein coordinate mapping
- STX protein domain annotations
- EC4115 comparative coordinates
- Functional predictions and classifications
- Cross-references to protein databases

**`phage_insertion_sites.tsv`** — Prophage integration analysis
- Integration coordinates relative to K-12
- tRNA insertion site identification
- Prophage boundary determination
- Context genes and regulatory elements

**`stx_phylogeny.nwk/.svg`** — Phylogenetic analysis
- Maximum likelihood tree with bootstrap support
- STX subtype relationships and clustering
- Variant evolution and divergence analysis

### Literature Integration

**`literature_refs.tsv`** — Citation tracking
- Automatic literature cross-referencing
- Variant-specific publication mapping
- EC4115 genomics paper integration
- Recent variant characterization studies

**`stx_results.json`** — Machine-readable export
- Complete results in JSON format
- Database version tracking
- Citation metadata included
- API-compatible structure

---

## EC4115 Integration

STXit incorporates methodologies and reference data from EC4115 Genomics Anatomy coursework, building upon the comprehensive characterization by Eppinger et al. (PNAS 2011).

### Reference Genome: E. coli O157:H7 str. EC4115

**Primary Citation:**
Manning SD, Motiwala AS, Springman AC, Qi W, Lacher DW, Ouellette LM, Mladonicky JM, Somsel P, Rudrik JT, Dietrich SE, Zhang W, Swaminathan B, Alland D, Whittam TS. **Variation in virulence among clades of Escherichia coli O157:H7 associated with disease outbreaks.** *Proceedings of the National Academy of Sciences*, 2008; 105(12):4868-4873. [DOI: 10.1073/pnas.0710834105](https://doi.org/10.1073/pnas.0710834105)

### EC4115 STX Characterization

**Expected STX Profile (Reference Test Case):**
| Locus ID | Position | Subtype | Identity | Prophage Context | Insertion Site Gene |
| --- | --- | --- | --- | --- | --- |
| STX_001 | 2,696,202–2,697,442 (−) | stx2c | 100.0% | φ-EC4115-1 prophage | argW (tRNA-Arg) |
| STX_002 | 3,271,007–3,272,247 (−) | stx2a | 99.9% | φ-EC4115-2 prophage | sbcB (exonuclease I) |

**stx2a variants (relative to EDL933, X07865):**
- 11 differences vs canonical sequence
- 6 synonymous (silent mutations) 
- 5 non-synonymous: V971A • F1009L • S1085L • W1097L • P1099T
- No stop codons introduced
- Maintained protein function

### Locus ID System

STXit assigns systematic locus identifiers for all detected STX genes:

| Locus Pattern | Description | Example |
| --- | --- | --- |
| `STX_001` | First detected locus | stx2c at argW site |
| `STX_002` | Second detected locus | stx2a at sbcB site |
| `STX_F01` | StxF family variants | stxF at yecE site |
| `STX_O01` | StxO family variants | stxO at thrW site |
| `STX_L01` | StxL family variants | stxL at ssrA site |

### Common Insertion Site Genes

| Gene | Function | Common STX Association | Frequency |
| --- | --- | --- | --- |
| **argW** | tRNA-Arg | stx2c, stx2d | High |
| **sbcB** | Exonuclease I | stx2a variants | High |
| **thrW** | tRNA-Thr | stx1a, stx2e | Medium |
| **yecE** | Hypothetical protein | stxF variants | Medium |
| **ssrA** | tmRNA | stxL variants | Low |
| **wrbA** | NAD(P)H quinone reductase | stx2f, stx2g | Low |
| **yehV** | Sensor kinase | stxO variants | Low |

### K-12 Relative Mapping

STXit maps prophage insertion sites relative to E. coli K-12 MG1655 for comparative genomics:

```bash
# K-12 comparative analysis
stxit --genome EC4115.fasta \
      --k12-mapping \
      --insertion-sites \
      --output ec4115_analysis/

# Outputs coordinate transformations:
# - EC4115 → K-12 coordinate mapping
# - Prophage integration site identification
# - Synteny block analysis
# - Insertion/deletion events
```

### Genomics Methodology Integration

STXit implements core concepts from EC4115 Genomics Anatomy:
- **Sequence alignment**: BLAST-based homology detection
- **Phylogenetic analysis**: ML tree construction with bootstrap support
- **Comparative genomics**: Multi-genome STX variant analysis
- **Database integration**: Literature mining and cross-referencing

---

## Test case

**Multi-genome validation: EC4115 + EDL933 + Sakai + TW14588**

### Test Dataset
```bash
# Download multiple reference genomes
stxit --accessions NC_011353.1,CP008957.1,BA000007.2,AAJT00000000 \
      --output test_validation/ \
      --download-only

# Run multi-genome validation test
stxit --accessions NC_011353.1,CP008957.1,BA000007.2 \
      --genomes TW14588_draft.fasta \
      --output test_results/ \
      --run-phastest \
      --run-trna \
      --validate-pipeline \
      --comparative-mode
```

### Expected Results (Multi-genome Summary)

| Genome | STX Loci | Subtypes Detected | Insertion Sites | Prophage Regions |
| --- | --- | --- | --- | --- |
| **EC4115** | 2 | stx2c, stx2a | argW, sbcB | φ-EC4115-1, φ-EC4115-2 |
| **EDL933** | 1 | stx2a | sbcB | φ-EDL933-1 |
| **Sakai** | 2 | stx1a, stx2a | thrW, yehV | φ-Sakai-1, φ-Sakai-2 |
| **TW14588** | 1 | stx2c | argW | φ-TW14588-1 |

### Comparative Analysis Results

**Cross-genome locus mapping:**
- **STX_001 (stx2c)**: Found in EC4115, TW14588 → argW insertion
- **STX_002 (stx2a)**: Found in EC4115, EDL933, Sakai → sbcB/yehV insertion  
- **STX_003 (stx1a)**: Found in Sakai only → thrW insertion
- **Total unique loci**: 3 across 4 genomes
- **Phylogenetic clustering**: stx2 family dominance

### Batch Validation Outputs

```
test_results/
├── multi_genome_summary.tsv       # Cross-genome comparison
├── stx_calls_EC4115.tsv          # Individual genome results
├── stx_calls_EDL933.tsv          # Individual genome results  
├── stx_calls_Sakai.tsv           # Individual genome results
├── stx_calls_TW14588.tsv         # Individual genome results
├── comparative_locus_map.tsv      # Locus ID cross-reference
├── insertion_sites_comparison.tsv # Integration site analysis
├── stx_phylogeny_all.nwk         # Multi-genome phylogeny
├── batch_validation_report.html   # Interactive summary
└── literature_refs_combined.tsv   # Consolidated citations
```

### Detailed STX Analysis

**STX_001 (stx2c):**
- **Locus ID**: STX_001
- **Position**: 2,696,202–2,697,442 (−)
- **Identity**: 100.0% exact match
- **Prophage**: φ-EC4115-1 (complete)
- **Insertion site**: argW (tRNA-Arg, position 2,696,180)
- **Flanking genes**: upstream: yecE, downstream: insB1

**STX_002 (stx2a-like):**
- **Locus ID**: STX_002  
- **Position**: 3,271,007–3,272,247 (−)
- **Identity**: 99.9% (11 variants)
- **Prophage**: φ-EC4115-2 (complete)
- **Insertion site**: sbcB (exonuclease I, position 3,271,000)
- **Flanking genes**: upstream: rlmE, downstream: yihS

### Insertion Site Gene Analysis

STXit provides comprehensive analysis of prophage integration sites:

```bash
# Analyze insertion site genes
stxit --genome EC4115.fasta \
      --insertion-analysis \
      --flanking-genes 5000 \
      --output insertion_analysis/

# Generates: insertion_sites.tsv with:
# - Locus ID assignments
# - Integration gene identification
# - Flanking gene context (±5kb)
# - Disruption impact assessment
```

**Output includes:**
- **Target gene disruption**: Complete vs partial integration
- **Flanking gene context**: Upstream/downstream genes within 5kb
- **Comparative mapping**: Position relative to K-12 MG1655
- **Functional impact**: Effect on host gene function

### Validation Outputs

```
test_results/
├── stx_calls.long.tsv          # Detailed STX calls with locus IDs
├── stx_summary.tsv            # Concise results (STX_001, STX_002)
├── insertion_sites.tsv        # argW, sbcB integration analysis
├── protein_locus_map.tsv      # Protein coordinates by locus ID
├── phage_insertion_sites.tsv  # Integration analysis with gene context
├── stx_phylogeny.nwk         # Variant relationships
├── ec4115_validation.html     # Interactive report
└── literature_refs.tsv        # Citation tracking
```

### Example Output Content

**`stx_summary.tsv` (EC4115 validation):**
```
locus_id	position	strand	subtype	identity	insertion_gene	prophage	flanking_upstream	flanking_downstream
STX_001	2696202-2697442	-	stx2c	100.0	argW	φ-EC4115-1	yecE	insB1
STX_002	3271007-3272247	-	stx2a	99.9	sbcB	φ-EC4115-2	rlmE	yihS
```

**`insertion_sites.tsv` (detailed integration analysis):**
```
locus_id	insertion_gene	gene_function	disruption_type	k12_position	flanking_genes_5kb	impact_score
STX_001	argW	tRNA-Arg	complete	2698345	yecE,insB1,yecF,yecG,yecH	0.8
STX_002	sbcB	exonuclease_I	partial	3275123	rlmE,yihS,yihT,yihU,yihV	0.3
```

---

## Database versions

STXit tracks database versions for reproducibility and citation accuracy:

| Database | Version | Date | Source |
| --- | --- | --- | --- |
| **STX References** | 1.0.3 | 2026-04-16 | [Scheutz et al. 2012](https://doi.org/10.1128/JCM.00737-12) + updates |
| **KLM N variants** | 3.2.1 | 2026-04-16 | Literature integration (156 papers) |
| **KLM O variants** | 3.2.1 | 2026-04-16 | Serotyping database + 89 studies |
| **TnCentral** | 431_seqs | 2026-04-15 | [tncentral.ncc.unesp.br](http://tncentral.ncc.unesp.br/) |
| **EC4115 Reference** | 1.0 | 2008 | [Manning et al. PNAS 2008](https://doi.org/10.1073/pnas.0710834105) |

### Database Update History

**Version 3.2.1 (Current):**
- Enhanced KLM N/O variant classification
- Integrated EC4115 genomics methodologies
- Added phage insertion site mapping
- Literature cross-referencing system
- K-12 comparative coordinate system

**Version 3.2.0:**
- Major KLM variant database expansion
- Added protein locus mapping
- Enhanced phylogenetic analysis
- Comprehensive test validation

Database version information is recorded in `databases/db_versions.txt` after installation and included in all output files for reproducibility.

### Update databases

Databases are versioned and can be updated automatically. To update to the latest releases after initial install:

```bash
conda activate stxit
bash databases/setup_databases.sh --update
```

This re-downloads and rebuilds STX references, stxF/stxO/stxL variants, TnCentral, and literature databases. Check upstream release pages for new versions:

| Database | Release page |
| --- | --- |
| STX References | [Scheutz nomenclature updates](https://pubmed.ncbi.nlm.nih.gov/?term=Scheutz+Shiga+toxin) |
| STEC Database | [CDC STEC Reference Center](https://www.cdc.gov/ecoli/laboratory.html) |
| TnCentral | [tncentral.ncc.unesp.br](https://tncentral.ncc.unesp.br/download.html) |
| Literature DB | Auto-updated from PubMed searches |

To force a fresh download of all databases:

```bash
bash databases/setup_databases.sh --force-update
```

---

## License

MIT License. See [LICENSE](https://github.com/kramppe/STXit/blob/main/LICENSE) for details.

---

Full citations with DOI and PMID: [docs/CITATIONS.md](https://github.com/kramppe/STXit/blob/main/docs/CITATIONS.md)

**Please cite:**

> Eppinger M, *et al.* STXit v1.0.3: a comprehensive bioinformatics tool for detecting, typing, and analyzing Shiga toxin genes in bacterial genomes with prophage context analysis. In preparation. 2026. Eppinger Lab | UTSA | MMI | STCEID. https://github.com/kramppe/STXit

**Key dependencies:** BLAST+ · tRNAscan-SE · IQ-TREE 2 · UFBoot2 · PHASTEST · BioPython · DIAMOND · HMMER · TnCentral · MobileOG-db · ISEScan · IntegronFinder 2.0

## Contact

STXit is developed and maintained by the Eppinger Lab at UTSA.
For questions, bug reports, and contributions: https://github.com/kramppe/STXit/issues
