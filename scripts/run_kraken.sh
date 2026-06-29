#!/bin/bash
# Taxonomic classification of cleaned viral reads with Kraken2 + k2_viral DB.
# Two confidence settings: 0.1 (primary) and 0.5 (robustness check).
set -euo pipefail
ROOT=/insomnia001/depts/pmg/users/na2933/pivotal-zephyr
source /insomnia001/depts/pmg/users/na2933/miniforge/etc/profile.d/conda.sh
conda activate kraken2
cd "$ROOT"
DB="$ROOT/kraken/db"
THREADS=${THREADS:-4}
mkdir -p kraken

for clean in data/clean/*.clean.fasta; do
  pool=$(basename "$clean" .clean.fasta)
  # primary pass: confidence 0.1, --use-names so per-read output carries the
  # taxon name (used for embedding labels + coverage read selection)
  kraken2 --db "$DB" --threads "$THREADS" --confidence 0.1 --use-names \
    --report "kraken/${pool}.report" --output "kraken/${pool}.out" "$clean"
  # robustness pass: confidence 0.5, report only
  kraken2 --db "$DB" --threads "$THREADS" --confidence 0.5 \
    --report "kraken/${pool}.c50.report" --output /dev/null "$clean"
  echo "KRAKEN $pool done"
done
echo "KRAKEN_DONE"
