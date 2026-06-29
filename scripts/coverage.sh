#!/bin/bash
# Genome coverage of common respiratory viruses: map each pool's cleaned reads
# to each reference in refs/ (minimap2 map-ont), record per-position depth for
# pairs with >0 mapped reads.
set -euo pipefail
ROOT=/insomnia001/depts/pmg/users/na2933/pivotal-zephyr
export PATH=/insomnia001/depts/pmg/users/na2933/metabat2_baseline/metabat2_env/bin:$PATH
cd "$ROOT"
mkdir -p coverage
: > coverage/mapped_counts.tsv

for ref in refs/*.fa; do
  [ -e "$ref" ] || { echo "no refs found"; exit 1; }
  v=$(basename "$ref" .fa)
  for clean in data/clean/*.clean.fasta; do
    pool=$(basename "$clean" .clean.fasta)
    bam="coverage/${v}__${pool}.bam"
    minimap2 -ax map-ont -t 4 "$ref" "$clean" 2>/dev/null \
      | samtools sort -o "$bam" - 2>/dev/null
    samtools index "$bam"
    mapped=$(samtools view -c -F 0x904 "$bam")
    printf "%s\t%s\t%s\n" "$v" "$pool" "$mapped" >> coverage/mapped_counts.tsv
    if [ "$mapped" -gt 0 ]; then
      samtools depth -a "$bam" > "coverage/${v}__${pool}.depth"
    fi
    rm -f "$bam" "$bam.bai"
  done
done
echo "COVERAGE_DONE"
