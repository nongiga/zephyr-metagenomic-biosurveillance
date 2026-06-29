#!/usr/bin/env python3
"""QC for Zephyr respiratory viral reads (stdlib only).

Reads are FASTA (no per-base quality), already basecalled / human-depleted /
virus-filtered upstream by Zephyr, so adapter/quality trimming is N/A. We apply
a light length filter + exact-duplicate removal and record before/after stats.
"""
import glob
import gzip
import os

ROOT = "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr"
DATA = os.path.join(ROOT, "data")
CLEAN = os.path.join(DATA, "clean")
RESULTS = os.path.join(ROOT, "results")
MINLEN = 100
os.makedirs(CLEAN, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)


def read_fasta_gz(path):
    op = gzip.open if path.endswith(".gz") else open
    rid, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(buf)
                rid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if rid is not None:
        yield rid, "".join(buf)


def stats(lengths):
    if not lengths:
        return dict(n=0, sum=0, min=0, mean=0.0, max=0, n50=0)
    s = sorted(lengths, reverse=True)
    tot = sum(s)
    half, run, n50 = tot / 2, 0, s[-1]
    for L in s:
        run += L
        if run >= half:
            n50 = L
            break
    return dict(n=len(s), sum=tot, min=min(s), mean=round(tot / len(s), 1),
                max=max(s), n50=n50)


def main():
    raw_rows, clean_rows = [], []
    for gz in sorted(glob.glob(os.path.join(DATA, "*.respiratory.fasta.gz"))):
        pool = os.path.basename(gz).replace(".respiratory.fasta.gz", "")
        raw_len, kept_len, seen = [], [], set()
        out = open(os.path.join(CLEAN, f"{pool}.clean.fasta"), "w")
        for rid, seq in read_fasta_gz(gz):
            raw_len.append(len(seq))
            if len(seq) < MINLEN or seq in seen:
                continue
            seen.add(seq)
            out.write(f">{rid}\n{seq}\n")
            kept_len.append(len(seq))
        out.close()
        r, c = stats(raw_len), stats(kept_len)
        raw_rows.append((pool, r)); clean_rows.append((pool, c))
        print(f"{pool}: raw {r['n']} -> clean {c['n']} reads "
              f"(len {c['min']}-{c['max']}, mean {c['mean']}, N50 {c['n50']})")

    hdr = "pool\tnum_seqs\tsum_len\tmin_len\tmean_len\tmax_len\tN50\n"
    for fname, rows in [("raw_stats.tsv", raw_rows), ("clean_stats.tsv", clean_rows)]:
        with open(os.path.join(RESULTS, fname), "w") as fh:
            fh.write(hdr)
            for pool, s in rows:
                fh.write(f"{pool}\t{s['n']}\t{s['sum']}\t{s['min']}\t"
                         f"{s['mean']}\t{s['max']}\t{s['n50']}\n")
    print("QC_DONE")


if __name__ == "__main__":
    main()
