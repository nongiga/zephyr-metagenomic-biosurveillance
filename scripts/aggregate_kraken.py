#!/usr/bin/env python3
"""Aggregate per-pool Kraken2 reports into a species x pool matrix + confidence
summaries, and render a composition heatmap.

Kraken2 report columns (tab-separated):
  0 pct  1 clade_reads  2 taxon_reads  3 rank  4 taxid  5 name (indented)
"""
import glob
import os
import re
import sys
import pandas as pd

ROOT = "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr"
KRAKEN = os.path.join(ROOT, "kraken")
RESULTS = os.path.join(ROOT, "results")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)


def parse_report(path):
    """Return (species_clade_counts {name: clade_reads}, total, unclassified)."""
    species, total, unclassified = {}, 0, 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            clade = int(f[1]); rank = f[3]; name = f[5].strip()
            if rank == "U":
                unclassified = clade
            if rank == "R" and f[4] == "1":
                total += clade  # classified clade under root
            if rank == "S":  # species level only (exclude S1/S2 strains)
                species[name] = species.get(name, 0) + clade
    total = total + unclassified  # grand total reads in the pool
    return species, total, unclassified


def pool_name(path):
    return re.sub(r"\.report$", "", os.path.basename(path))


def main():
    reports = sorted(glob.glob(os.path.join(KRAKEN, "*.report")))
    reports = [r for r in reports if not r.endswith(".c50.report")]
    if not reports:
        sys.exit("no kraken reports found")

    mat, classified_rows = {}, []
    for r in reports:
        pool = pool_name(r)
        sp, total, uncl = parse_report(r)
        mat[pool] = sp
        classified = total - uncl
        classified_rows.append({
            "pool": pool, "total_reads": total,
            "classified_reads": classified,
            "pct_classified": round(100 * classified / total, 2) if total else 0.0,
        })

    # species x pool matrix (clade read counts, conf 0.1)
    df = pd.DataFrame(mat).fillna(0).astype(int)
    df = df.loc[df.sum(axis=1).sort_values(ascending=False).index]
    df.to_csv(os.path.join(RESULTS, "species_by_pool.tsv"), sep="\t")

    classified = pd.DataFrame(classified_rows).sort_values("pool")
    classified.to_csv(os.path.join(RESULTS, "classified_summary.tsv"),
                      sep="\t", index=False)

    # confidence summary: read support, #pools, and conf-0.5 survival
    c50 = {}
    for r in glob.glob(os.path.join(KRAKEN, "*.c50.report")):
        sp, _, _ = parse_report(r)
        for k, v in sp.items():
            c50[k] = c50.get(k, 0) + v
    conf = pd.DataFrame({
        "species": df.index,
        "reads_conf0.1": df.sum(axis=1).values,
        "n_pools": (df > 0).sum(axis=1).values,
        "reads_conf0.5": [c50.get(s, 0) for s in df.index],
    })
    conf["frac_surviving_conf0.5"] = (
        conf["reads_conf0.5"] / conf["reads_conf0.1"].replace(0, pd.NA)
    ).round(3)
    conf.to_csv(os.path.join(RESULTS, "species_confidence.tsv"),
                sep="\t", index=False)

    print("== classified summary ==")
    print(classified.to_string(index=False))
    print("\n== top 15 species (reads, #pools, conf0.5 survival) ==")
    print(conf.head(15).to_string(index=False))

    # composition heatmap (top 25 species)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        top = df.head(25)
        plt.figure(figsize=(max(8, 0.7 * df.shape[1]), 10))
        sns.heatmap(np.log10(top + 1), cmap="viridis",
                    cbar_kws={"label": "log10(reads + 1)"})
        plt.title("Viral species composition per pool (Kraken2, conf 0.1)")
        plt.ylabel("species"); plt.xlabel("pool")
        plt.tight_layout()
        plt.savefig(os.path.join(FIGS, "composition_heatmap.png"), dpi=150)
        print("\nwrote figures/composition_heatmap.png")
    except Exception as e:  # plotting is non-critical
        print(f"[warn] heatmap skipped: {e}")


if __name__ == "__main__":
    main()
