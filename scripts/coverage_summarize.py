#!/usr/bin/env python3
"""Summarize per-position depth files into breadth + mean depth per (virus, pool),
and plot depth-along-genome for the best-covered case."""
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr"
COV = os.path.join(ROOT, "coverage")
RESULTS = os.path.join(ROOT, "results")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)


def main():
    rows, depth_cache = [], {}
    for f in sorted(glob.glob(os.path.join(COV, "*.depth"))):
        base = os.path.basename(f)[:-6]
        virus, pool = base.split("__", 1)
        d = np.loadtxt(f, usecols=2, dtype=int, ndmin=1)
        if d.size == 0:
            continue
        reflen = d.size
        breadth = round(100 * np.count_nonzero(d) / reflen, 2)
        mean_depth = round(float(d.mean()), 2)
        rows.append({"virus": virus, "pool": pool, "reflen": reflen,
                     "breadth_pct": breadth, "mean_depth": mean_depth,
                     "n_pos_covered": int(np.count_nonzero(d))})
        depth_cache[base] = d

    if not rows:
        print("no coverage (no reads mapped to any reference)")
        return
    df = pd.DataFrame(rows).sort_values(["breadth_pct", "mean_depth"],
                                        ascending=False)
    df.to_csv(os.path.join(RESULTS, "coverage.tsv"), sep="\t", index=False)
    print(df.to_string(index=False))

    # depth-along-genome for the best-covered (virus, pool)
    best = df.iloc[0]
    key = f"{best['virus']}__{best['pool']}"
    d = depth_cache[key]
    plt.figure(figsize=(11, 3.5))
    plt.fill_between(np.arange(len(d)), d, step="mid", color="steelblue")
    plt.xlabel("genome position (bp)"); plt.ylabel("depth")
    plt.title(f"{best['virus']} in {best['pool']} — "
              f"breadth {best['breadth_pct']}%, mean depth {best['mean_depth']}x")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "coverage_best.png"), dpi=150)
    print(f"wrote figures/coverage_best.png ({key})")


if __name__ == "__main__":
    main()
