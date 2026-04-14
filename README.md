# STXit v1.0

**STXit** is a standalone subtype-calling, prophage-context, and reporting workflow for assembled bacterial genomes focused on Shiga toxin (`stx`) loci, subtype reference matching, prophage localization, and insertion-site interpretation. The tool supports single FASTA, multi-FASTA, annotated GenBank, unannotated draft assemblies, and accession-based retrieval. STXit can also update and version its subtype reference database and known insertion-site catalog with timestamped outputs.

Developed and maintained by the Eppinger Lab at the University of Texas at San Antonio (UTSA).  
GitHub 2026.

## Contents

- Pipeline overview
- Dependencies
- Installation
- Input files
- Query modes
- Quick start
- Output files
- STX subtype logic
- Prophage and insertion-site logic
- Variant analysis
- Tree building
- Database update mode
- Citations
- Contact

## Pipeline overview

| Step | Module / Script | Role |
|------|------------------|------|
| 0 | `stx_db.py` | Load, validate, and optionally update the STX subtype database |
| 1 | `stx_typing.py` | Detect `stx` loci and assign nearest subtype references |
| 2 | `annotation.py` | Derive reference-centered and optional query-side gene annotation |
| 3 | `phage_prediction.py` | Parse PHASTEST results or submit genomes for prophage prediction |
| 4 | `trna_scan.py` | Detect tRNA loci or ingest precomputed tRNAscan-SE results |
| 5 | `phage_insertion.py` / `neighborhood.py` | Infer prophage insertion context, left/right flanking loci, and backbone comparison |
| 6 | `plotting.py` | Draw closed-genome phage/STX context figures |
| 7 | `reporting.py` | Produce long tables, summary tables, overlap tables, and metadata reports |
| 8 | `tree_utils.py` | Build small gene trees for non-exact STX variants and optional context trees |

## Dependencies

| Tool | Version | Install |
|------|---------|---------|
| Python | ≥ 3.10 | — |
| Biopython | ≥ 1.80 | `pip install biopython` |
| BLAST+ | ≥ 2.12 | `conda install -c bioconda blast` |
| IQ-TREE 2 | ≥ 2.0 | `conda install -c bioconda iqtree` |
| FastTree 2 | ≥ 2.1 | `conda install -c bioconda fasttree` |
| tRNAscan-SE | ≥ 2.0 | `conda install -c bioconda trnascan-se` |

Optional:
- `ncbi-datasets-cli` for accession-based retrieval
- `entrez-direct` for nucleotide record retrieval
- PHASTEST submission or parsed output ingestion

## Installation

### Option A — manual Conda environment

```bash
conda create -n stxit_v1 python=3.10
conda activate stxit_v1
conda install -c bioconda blast iqtree fasttree entrez-direct trnascan-se
conda install -c conda-forge ncbi-datasets-cli
pip install biopython pandas openpyxl matplotlib
```

### Option B — `environment.yml`

```bash
conda env create -f environment.yml
conda activate stxit_v1
```

### Install STXit

```bash
git clone https://github.com/kramppe/STXit.git
cd STXit
pip install .
```

## Input files

### 1. STX subtype reference FASTA

The subtype database FASTA stores canonical or curated subtype references. FASTA headers should include accession whenever possible.

Recommended header format:

```text
>reference_id|accession|stx_family|stx_subtype|variant_name|strain|citation_key
```

### 2. STX reference metadata

Each sequence in the subtype FASTA should have a metadata entry. The subtype reference database is literature-linked. Each FASTA entry should map to a metadata row and a citation row documenting the subtype designation, accession, and source publication used for curation.

```tsv
reference_id	stx_family	stx_subtype	variant_name	accession	strain	length	is_primary_reference	citation_key
stx2a_ref_1	stx2	stx2a	canonical	X07865	EDL933	1234	true	Scheutz2012
stx2k_ref_1	stx2	stx2k	canonical	MN200195	12GZSW01	1234	true	Hughes2020_Stx2k
```

### 3. STX citation table

```tsv
citation_key	authors	year	title	journal	doi	pmid	pmcid	note
Scheutz2012	Scheutz F et al.	2012	Multicenter evaluation of a sequence-based protocol for subtyping Shiga toxins and standardizing Stx nomenclature	J Clin Microbiol	10.1128/JCM.00860-12	22760050	PMC3421821	Standardized subtype nomenclature
```

### 4. Known insertion-site catalog

STXit supports a configurable catalog of known chromosomal insertion sites. These sites may be defined from K-12, EC4115, EDL933, Sakai, or other non-STX backbone references, and each row may include the source accession used for curation. More known insertion sites may be added over time and the catalog can be updated.

```tsv
site_name	feature_type	match_gene	match_locus_tag	backbone_accession	backbone_label	category	note
wrbA-like	CDS	wrbA		U00096.3	Ecoli_K12_MG1655	backbone_named_site	Common Stx phage insertion site
argW-associated	tRNA	argW		U00096.3	Ecoli_K12_MG1655	trna_site	Common Stx phage insertion site
sbcB-like	CDS	sbcB		U00096.3	Ecoli_K12_MG1655	backbone_named_site	Known Stx phage insertion site
```

### 5. Optional custom phage category file

```tsv
sample	region_id	custom_name	category	note
EC4115	phage_1	Stx2a_phi_A	stx2a_prophage	manual review
```

### 6. Optional backbone comparison reference

By default, STXit compares prophage insertion context to an E. coli K-12 backbone. This may be replaced with another non-STX *Escherichia coli* or *Shigella* reference genome.

### 7. Optional subject metadata for database mode

```tsv
subject_id	label	genbank_file	group	note
NZ_CP008957.1	EDL933	annotations/EDL933.gbk	outbreak1	complete genome
```

## Query modes

### Genome mode

```bash
python stxit.py \
  --query genome1.gbk genome2.fasta genome3.fasta \
  --outdir results
```

### Accession mode

```bash
python stxit.py \
  --query-accessions NC_002695 NZ_CP008957 NC_011353 \
  --outdir results
```

### Database update mode

```bash
python stxit.py \
  --update-db \
  --db-fasta config/comprehensive_stx_typing.fasta \
  --db-metadata config/stx_references.tsv \
  --db-citations config/stx_citations.tsv \
  --known-sites config/known_insertion_sites.tsv \
  --timestamp-db \
  --db-output-dir config/db_updates
```

## Quick start

### Example 1 — subtype + prophage context on genomes

```bash
python stxit.py \
  --query isolate1.gbk isolate2.fasta isolate3.gbk \
  --run-trnascan \
  --phastest-mode parse \
  --outdir run_stx
```

### Example 2 — default K-12 backbone comparison

```bash
python stxit.py \
  --query isolate1.gbk isolate2.gbk \
  --compare-backbone default_k12 \
  --outdir run_backbone
```

### Example 3 — custom backbone comparison

```bash
python stxit.py \
  --query isolate1.gbk \
  --backbone-comparison-genbank K12_MG1655.gbk \
  --outdir run_custom_backbone
```

## Output files

- `stx_calls.long.tsv`
- `stx_summary.tsv`
- `STX_stats.tsv`
- `phage_regions.tsv`
- `stx_phage_overlap.tsv`
- `phage_insertion_context.tsv`
- `phage_insertion_coords.tsv`
- `trna_sites.tsv`
- `run_metadata.tsv`

For non-exact subtype matches:
- `stx_variant_differences.tsv`
- `stx_alignment_overview.tsv`
- `stx_variant_tree.fasta`
- `stx_variant_tree.nwk`
- `stx_variant_tree.svg`

For closed genomes:
- `sample_phage_map.svg`
- `sample_phage_map.png`

Database provenance outputs:
- `stx_references_YYYY-MM-DD.fasta`
- `stx_references_YYYY-MM-DD.tsv`
- `stx_citations_YYYY-MM-DD.tsv`
- `known_insertion_sites_YYYY-MM-DD.tsv`
- `db_update_log_YYYY-MM-DD.txt`

## STX subtype logic

- If a detected STX hit is an exact 100% identity and 100% coverage match to a subtype reference, report the subtype and best hit without generating a variant-difference table by default.
- If the hit is not a 100% exact match, identify the nearest subtype reference and report SNP-level differences, alignment overview, and optional gene-level tree placement.
- If the sequence is divergent enough to be unresolved, flag it as a candidate novel or unresolved STX variant.

## Prophage and insertion-site logic

STXit combines subtype results with prophage-context analysis.

Per prophage region, report:
- prophage boundaries
- prophage completeness/class
- STX overlap
- left and right chromosomal flanking loci
- nearby tRNA loci
- attachment-site fields when available
- integrase presence
- insertion-site interpretation relative to the default K-12 backbone or a user-provided alternative backbone
- named backbone insertion-site matches such as `wrbA-like`, `argW-associated`, `sbcB-like`, and additional curated sites when recognized

This is not limited to integrase; the key goal is to identify the true left and right chromosomal loci bordering the prophage insertion site, determine which normal backbone feature or interval is disrupted, and preserve named site calls when known.

## Variant analysis

For non-exact subtype hits, STXit should report:
- nearest reference
- nearest reference accession
- percent identity
- percent coverage
- SNP-level differences relative to nearest reference
- compact alignment summary
- optional phylogenetic placement among subtype references

If the query is an exact match, these extra outputs are not required by default.

## Tree building

STXit may build a small subtype-reference tree when a non-exact STX variant is detected.

Priority:
1. IQ-TREE 2
2. FastTree 2
3. Biopython NJ fallback

The goal is to show phylogenetic placement of the variant relative to curated subtype references.

## Database update mode

STXit supports updating the STX subtype database and known insertion-site catalog using a provided FASTA and metadata/citation tables.

Recommended behavior:
- validate sequence headers
- validate subtype labels
- validate accession-linked metadata
- write timestamped FASTA
- write timestamped metadata table
- write timestamped citation table
- validate and optionally update `known_insertion_sites.tsv`
- preserve accession-linked backbone site provenance
- write a changelog

## Citations

Please cite the specific subtype reference source, backbone reference, prophage tool, and major software dependencies used in your workflow.

### STX subtype nomenclature and references
- Scheutz F, Teel LD, Beutin L, et al. Multicenter evaluation of a sequence-based protocol for subtyping Shiga toxins and standardizing Stx nomenclature. *J Clin Microbiol*. 2012;50(9):2951–2963. PMID: 22760050. PMCID: PMC3421821. DOI: 10.1128/JCM.00860-12.
- Persson S, Olsen KEP, Ethelberg S, Scheutz F. Subtyping method for *Escherichia coli* Shiga toxin 2 variants and correlations to clinical manifestations. *J Clin Microbiol*. 2007;45(6):2020–2024. PMID: 17446326. PMCID: PMC1933035. DOI: 10.1128/JCM.02591-06.
- Gill A, Dussault F, McMahon T, et al. Characterization of atypical Shiga toxin gene sequences and description of Stx2j, a new subtype. *J Clin Microbiol*. 2022;60(3):e02229-21. PMID: 35225693. PMCID: PMC8925903. DOI: 10.1128/JCM.02229-21.
- Hughes AC, Zhang Y, Bai X, et al. Structural and functional characterization of Stx2k, a new subtype of Shiga toxin 2. *Microorganisms*. 2020;8(1):4. PMID: 31892121. PMCID: PMC7022315. DOI: 10.3390/microorganisms8010004.
- Lindsey RL, Prasad A, Feldgarden M, et al. Identification and characterization of ten *Escherichia coli* strains encoding novel Shiga toxin 2 subtypes, Stx2n as well as Stx2j, Stx2m, and Stx2o, in the United States. *Microorganisms*. 2023;11(10):2561. PMID: 37894219. PMCID: PMC10608928. DOI: 10.3390/microorganisms11102561.
- Bai X, Scheutz F, Dahlgren HM, et al. Genomic characterization of *Escherichia coli* O8 strains producing Shiga toxin 2l subtype. *Microorganisms*. 2022;10(6):1188. PMID: 35744526. PMCID: PMC9227347. DOI: 10.3390/microorganisms10061188.

### Prophage and tRNA tools
- Wishart DS, Han S, Saha S, Oler E, Peters H, Grant JR, Stothard P, Gautam V. PHASTEST: faster than PHASTER, better than PHAST. *Nucleic Acids Res.* 2023;51(W1):W443–W450. DOI: 10.1093/nar/gkad382.
- Chan PP, Lin BY, Mak AJ, Lowe TM. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res*. 2021;49(16):9077–9096. PMID: 34417604. PMCID: PMC8450103. DOI: 10.1093/nar/gkab688.

### Software dependencies
- Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–1423. PMID: 19304878. PMCID: PMC2682512. DOI: 10.1093/bioinformatics/btp163.
- Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. PMID: 20003500. PMCID: PMC2803857. DOI: 10.1186/1471-2105-10-421.
- Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530–1534. PMID: 32011700. PMCID: PMC7182206. DOI: 10.1093/molbev/msaa015.
- Price MN, Dehal PS, Arkin AP. FastTree 2—approximately maximum-likelihood trees for large alignments. *PLoS One*. 2010;5(3):e9490. PMID: 20224823. PMCID: PMC2835736. DOI: 10.1371/journal.pone.0009490.

### Tool citation
Please cite:

> Eppinger M, *et al.* STXit v1.0: a standalone subtype-calling, prophage-context, and reporting workflow for Shiga toxin loci in assembled bacterial genomes. GitHub 2026. Eppinger Lab | UTSA | MMI | STCEID. https://github.com/kramppe/STXit

## Contact

STXit is developed and maintained by the Eppinger Lab at the University of Texas at San Antonio (UTSA).  
For questions, bug reports, and contributions, please use the repository issue tracker.
