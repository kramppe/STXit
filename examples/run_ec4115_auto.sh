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
