# Zephyr metagenomic biosurveillance — Pivotal work test (Parts 2 & 3)

Taxonomic classification and embedding-based clustering of vertebrate-virus reads from
SecureBio's **Zephyr** program (pooled nasal swabs, Oxford Nanopore long-read
metatranscriptomic sequencing). See [`memo.md`](memo.md) for the write-up.

## What this does

- **Part 2 — taxonomy & coverage.** Classify reads from 10 swab pools with **Kraken2**
  (`k2_viral` DB), report which viruses are present and how confident we are, and compute
  **genome coverage** (breadth + depth) for common respiratory viruses by mapping reads to
  RefSeq genomes with minimap2.
- **Part 3, Track B — embedding clustering.** Embed cleaned reads with **NTv3-650M**
  (a nucleotide language model) via the ctx-gene harness, then **PCA + KMeans** and compare
  the clusters to the Kraken taxonomy (cluster purity, homogeneity/completeness/V-measure).

## Data

10 of 315 pools from `https://data.securebio.org/zephyr/` `respiratory-reads/`, selected at
random with **seed 42** (no specific pools were given). The exact list is in
[`results/chosen_pools.txt`](results/chosen_pools.txt). Pools span **1–30,417 reads**.

## Pipeline (`scripts/`)

| Step | Script | Env / resource |
|------|--------|----------------|
| QC (length filter + dedup, stats) | `qc.sh` | `qc` (seqkit) |
| Kraken2 classification (conf 0.1 + 0.5) | `run_kraken.sh` | `kraken2` |
| Species×pool matrix + confidence + heatmap | `aggregate_kraken.py` | `ctx-gene` |
| Sample ≤1000 reads/pool + labels | `prep_embed.py` | stdlib |
| k-mer composition embedding (reported) | `kmer_embed.py` | local, stdlib+numpy |
| NTv3 LM embedding (8M; optional swap) | `embed_reads.py` | `ctx-gene`, GPU |
| minimap2 best-hit per-read labels | `besthit.py` | local, mappy |
| Embedding taxonomy-organization eval | `cluster.py` | local |
| Coverage (minimap2 + samtools depth) | `coverage.sh` | minimap2/samtools |
| Coverage summary + plot | `coverage_summarize.py` | `ctx-gene` |

All compute runs on the Columbia **insomnia** cluster via `srun -A pmg -p short`
(`SLURM_CONF=/etc/slurm/slurm.conf`; GPU via `--gres=gpu:1`). Kraken2 and a small `qc` env
were created with conda/bioconda; NTv3-650M is run offline from the cluster HF cache.

## Key choices

- **k2_viral** DB (RefSeq viral): correct scope since reads are pre-filtered to vertebrate
  viruses; small and fast.
- Reads are **FASTA** (no quality) → QC is length-filter + exact-dup removal; adapter/quality
  trimming was done upstream by Zephyr basecalling.
- Part 3B embedding is a **canonical k-mer composition spectrum** (k=4, 136-dim), per the task's
  "avoid a huge model" guidance; the **NTv3-650M** nucleotide-LM harness (`embed_reads.py`) is also
  implemented as a drop-in alternative (blocked here only by transient cluster filesystem slowness).
- Taxonomic organization of the embedding measured **directly** (clustering-free): same/diff-taxon
  **AUROC** + **kNN agreement** vs chance, against **minimap2 best-hit** labels, with a **pool**
  batch-effect control; KMeans homogeneity/completeness/AMI only as descriptors (not ARI/purity).

## Outputs (`results/`, `figures/`)

`chosen_pools.txt`, `raw_stats.tsv`/`clean_stats.tsv`, `species_by_pool.tsv`,
`classified_summary.tsv`, `species_confidence.tsv`, `coverage.tsv`,
`cluster_vs_taxonomy.tsv`; figures: composition heatmap, best-covered genome depth,
contingency heatmap, PCA scatter.
