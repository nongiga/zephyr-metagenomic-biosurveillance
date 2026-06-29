#!/usr/bin/env python3
"""Evaluate whether the read embedding is organized by taxonomy.

We have labels (minimap2 best-hit per read), so we measure the EMBEDDING's
label-structure directly rather than scoring an arbitrary clustering:

  PRIMARY (clustering-free, label-aware):
    - same- vs different-taxon AUROC over pairwise embedding distances
      (threshold-free generalization of pair-counting indices; chance = 0.5)
    - kNN leave-one-out label agreement (vs the chance baseline sum p_c^2)
  CONFOUND CONTROL:
    - the same two measures with POOL as the grouping (batch effect check)
  CLUSTERING DESCRIPTORS (KMeans, reported as descriptors only — no ARI/purity):
    - homogeneity / completeness / AMI

Reference labels = minimap2 best-hit (embed/besthit.tsv); metrics computed on the
mapped subset. Unmapped reads are kept only for the PCA figure (discovery view).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import (roc_auc_score, homogeneity_completeness_v_measure,
                             adjusted_mutual_info_score)

ROOT = os.environ.get("PZ_ROOT", "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr")
EMBED = os.path.join(ROOT, "embed")
RESULTS = os.path.join(ROOT, "results")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)
K = 10


def auroc_same_vs_diff(X, labels):
    """P(a same-label pair is closer than a different-label pair)."""
    d = pdist(X)                       # condensed pairwise distances
    li = np.array(labels)
    n = len(li)
    iu = np.triu_indices(n, 1)
    same = (li[iu[0]] == li[iu[1]]).astype(int)
    return roc_auc_score(same, -d)     # closer (smaller d) -> more likely same


def knn_agreement(X, labels, k=K):
    """Mean fraction of each point's k nearest neighbours sharing its label."""
    li = np.array(labels)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, idx = nn.kneighbors(X)
    agree = np.mean([np.mean(li[row[1:]] == li[row[0]]) for row in idx])
    p = pd.Series(li).value_counts(normalize=True).values
    return agree, float((p ** 2).sum())   # observed, chance baseline


def main():
    X = np.load(os.path.join(EMBED, "embeddings.npy"))
    ids = [l.strip() for l in open(os.path.join(EMBED, "ids.txt")) if l.strip()]
    meta = pd.read_csv(os.path.join(EMBED, "labels.tsv"), sep="\t").set_index("uid")
    bh = pd.read_csv(os.path.join(EMBED, "besthit.tsv"), sep="\t").set_index("uid")
    pool = meta.loc[ids, "pool"].values
    taxon = bh.loc[ids, "besthit"].values
    mapped = taxon != "unmapped"
    print(f"{X.shape[0]} reads, {X.shape[1]}-d; {mapped.sum()} mapped "
          f"({len(set(taxon[mapped]))} taxa)")

    # standardize + PCA (denoise); all distance/cluster work in this space
    Xs = StandardScaler().fit_transform(X)
    Xp = PCA(n_components=min(50, *Xs.shape), random_state=42).fit_transform(Xs)
    Xm, tax_m, pool_m = Xp[mapped], taxon[mapped], pool[mapped]

    # PRIMARY: taxonomy structure on the mapped reads
    tax_auroc = auroc_same_vs_diff(Xm, tax_m)
    tax_knn, tax_chance = knn_agreement(Xm, tax_m)
    # CONFOUND: pool structure on the same reads
    pool_auroc = auroc_same_vs_diff(Xm, pool_m)
    pool_knn, pool_chance = knn_agreement(Xm, pool_m)

    # CLUSTERING DESCRIPTORS (KMeans vs taxonomy, mapped reads)
    k = len(set(tax_m))
    km = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(Xm)
    h, c, _ = homogeneity_completeness_v_measure(tax_m, km)
    ami = adjusted_mutual_info_score(tax_m, km)

    rows = [
        ("taxon_AUROC_same_vs_diff", round(tax_auroc, 3)),
        ("taxon_kNN_agreement",      round(tax_knn, 3)),
        ("taxon_kNN_chance",         round(tax_chance, 3)),
        ("pool_AUROC_same_vs_diff",  round(pool_auroc, 3)),
        ("pool_kNN_agreement",       round(pool_knn, 3)),
        ("pool_kNN_chance",          round(pool_chance, 3)),
        ("KMeans_homogeneity",       round(h, 3)),
        ("KMeans_completeness",      round(c, 3)),
        ("KMeans_AMI",               round(ami, 3)),
        ("n_reads_total",            int(X.shape[0])),
        ("n_reads_mapped",           int(mapped.sum())),
        ("n_taxa",                   k),
    ]
    out = pd.DataFrame(rows, columns=["metric", "value"])
    out.to_csv(os.path.join(RESULTS, "embedding_eval.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))

    # PCA figure: colour by best-hit taxon (unmapped greyed) — incl. discovery view
    plt.figure(figsize=(9, 7))
    top = pd.Series(taxon[mapped]).value_counts().index.tolist()
    palette = dict(zip(top, sns.color_palette("tab10", len(top))))
    m = ~mapped
    plt.scatter(Xp[m, 0], Xp[m, 1], s=6, alpha=0.3, color="0.8", label="unmapped")
    for t in top:
        mm = taxon == t
        plt.scatter(Xp[mm, 0], Xp[mm, 1], s=14, alpha=0.85,
                    color=palette[t], label=t)
    plt.xlabel("PC1"); plt.ylabel("PC2")
    plt.title(f"k-mer read embeddings (PCA) by minimap2 best-hit\n"
              f"taxon AUROC={tax_auroc:.2f} (pool {pool_auroc:.2f}), "
              f"kNN={tax_knn:.2f} vs chance {tax_chance:.2f}")
    plt.legend(markerscale=2, fontsize=7, loc="best")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS, "pca_clusters.png"), dpi=150)
    print("wrote figures/pca_clusters.png, results/embedding_eval.tsv")


if __name__ == "__main__":
    main()
