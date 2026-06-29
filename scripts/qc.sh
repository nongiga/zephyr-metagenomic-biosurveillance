#!/bin/bash
# QC for Zephyr respiratory viral reads.
#
# The published reads are FASTA (no per-base quality), already basecalled,
# human-depleted and virus-filtered by Zephyr upstream. So adapter/quality
# trimming is not applicable here; we apply a light length filter + exact-dup
# removal and record before/after stats. ONT long reads => no pairing.
set -euo pipefail
ROOT=/insomnia001/depts/pmg/users/na2933/pivotal-zephyr
source /insomnia001/depts/pmg/users/na2933/miniforge/etc/profile.d/conda.sh
conda activate qc
cd "$ROOT"
mkdir -p data/clean results
MINLEN=${MINLEN:-100}

for gz in data/*.respiratory.fasta.gz; do
  pool=$(basename "$gz" .respiratory.fasta.gz)
  seqkit seq -g -m "$MINLEN" "$gz" 2>/dev/null \
    | seqkit rmdup -s 2>/dev/null \
    > "data/clean/${pool}.clean.fasta"
  echo "QC $pool -> $(grep -c '^>' data/clean/${pool}.clean.fasta) reads kept"
done

# before/after stats (tab-separated, with N50 etc.)
seqkit stats -T -a data/*.respiratory.fasta.gz   > results/raw_stats.tsv
seqkit stats -T -a data/clean/*.clean.fasta      > results/clean_stats.tsv
echo "QC_DONE"
