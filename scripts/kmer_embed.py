#!/usr/bin/env python3
"""Fast, dependency-free fallback embedding: canonical k-mer composition per read.

Each read -> vector of canonical k-mer frequencies (strand-independent, since ONT
cDNA strand orientation is arbitrary). This is the classic composition signal used
in metagenomic binning, and a sensible baseline when a heavy sequence LM is
unavailable. Output matches embed_reads.py: embeddings.npy + ids.txt.
"""
import os
import numpy as np

ROOT = os.environ.get("PZ_ROOT", "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr")
EMBED = os.path.join(ROOT, "embed")
K = 4
COMP = str.maketrans("ACGT", "TGCA")


def canonical_index():
    """Map each k-mer to a canonical-pair column index."""
    bases = "ACGT"
    kmers = [""]
    for _ in range(K):
        kmers = [p + b for p in kmers for b in bases]
    idx, col = {}, {}
    nxt = 0
    for km in kmers:
        rc = km.translate(COMP)[::-1]
        key = min(km, rc)
        if key not in col:
            col[key] = nxt
            nxt += 1
        idx[km] = col[key]
    return idx, nxt


def read_fasta(path):
    uid, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if uid is not None:
                    yield uid, "".join(buf)
                uid = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if uid is not None:
        yield uid, "".join(buf)


def main():
    idx, ncol = canonical_index()
    X, ids = [], []
    for uid, seq in read_fasta(os.path.join(EMBED, "reads.fasta")):
        seq = seq.upper()
        v = np.zeros(ncol, dtype=np.float32)
        n = 0
        for i in range(len(seq) - K + 1):
            j = idx.get(seq[i:i + K])
            if j is not None:  # skip k-mers with N
                v[j] += 1
                n += 1
        if n:
            v /= n
        X.append(v); ids.append(uid)
    X = np.vstack(X)
    np.save(os.path.join(EMBED, "embeddings.npy"), X)
    with open(os.path.join(EMBED, "ids.txt"), "w") as fh:
        fh.write("\n".join(ids) + "\n")
    print(f"kmer embeddings {X.shape} (k={K}, canonical) for {len(ids)} reads")


if __name__ == "__main__":
    main()
