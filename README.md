# STXit

**STXit** is a standalone subtype-calling, prophage-context, and reporting workflow for assembled bacterial genomes focused on Shiga toxin (`stx`) loci.

This merged snapshot includes:
- STX subtype calling from FASTA / multi-FASTA / GenBank input
- exact vs variant subtype logic
- codon-aware variant consequence reporting
- subtype-reference tree generation
- manual or automatic PHASTEST workflows
- manual or automatic tRNAscan-SE workflows
- prophage overlap and insertion-context reporting
- named insertion-site matching
- annotation-aware backbone comparison using a supplied GenBank file
- closed-genome plotting with PNG + SVG outputs

## Main outputs
- `stx_calls.long.tsv`
- `stx_summary.tsv`
- `STX_stats.tsv`
- `phage_regions.tsv`
- `stx_phage_overlap.tsv`
- `phage_insertion_context.tsv`
- `phage_insertion_coords.tsv`
- `stx_variant_differences.tsv`
- `stx_alignment_overview.tsv`
- `stx_variant_tree.index.tsv`
- `run_metadata.tsv`

## Installation

```bash
conda create -n stxit python=3.10
conda activate stxit
conda install -c bioconda blast trnascan-se
pip install -e .
```

Optional for tree building:
- IQ-TREE 2
- FastTree

## Example

Manual PHASTEST/tRNAscan input:

```bash
python stxit.py \
  --query EC4115.gbk \
  --db-fasta config/comprehensive_stx_typing.fasta \
  --phastest-tsv phastest_regions.tsv \
  --trnascan-tsv trnascan.tsv \
  --known-sites config/known_insertion_sites.tsv \
  --backbone-comparison-genbank K12_MG1655.gbk \
  --compare-backbone K12_MG1655 \
  --build-variant-trees \
  --plot-closed-genome-maps \
  --outdir results_ec4115
```

Automatic PHASTEST/tRNAscan mode:

```bash
python stxit.py \
  --query EC4115.gbk \
  --db-fasta config/comprehensive_stx_typing.fasta \
  --run-phastest \
  --phastest-email you@example.org \
  --run-trnascan \
  --known-sites config/known_insertion_sites.tsv \
  --build-variant-trees \
  --plot-closed-genome-maps \
  --outdir results_ec4115
```
