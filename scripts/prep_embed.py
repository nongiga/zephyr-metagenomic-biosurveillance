#!/usr/bin/env python3
"""Build the embedding read set: up to 1000 reads/pool (seed 42) sampled from ALL
cleaned reads, each tagged with its Kraken2 label (species name, or "unclassified").

Sampling from all reads — not just classified — lets the embedding speak to the
large unclassified fraction (reads from viruses divergent from / absent in the
RefSeq-viral DB). Writes embed/reads.fasta and embed/labels.tsv. Stdlib only."""
import glob
import os
import random
import re

ROOT = "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr"
CLEAN = os.path.join(ROOT, "data", "clean")
KRAKEN = os.path.join(ROOT, "kraken")
EMBED = os.path.join(ROOT, "embed")
os.makedirs(EMBED, exist_ok=True)
CAP = 1000
random.seed(42)

NAME_RE = re.compile(r"^(.*) \(taxid (\d+)\)$")


def read_fasta(path):
    seqs, rid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rid is not None:
                    seqs[rid] = "".join(buf)
                rid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if rid is not None:
        seqs[rid] = "".join(buf)
    return seqs


def label_map(out_path):
    """read_id -> (taxon, taxid, classified_bool) from a kraken --use-names .out"""
    m = {}
    with open(out_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            status, rid = f[0], f[1]
            mm = NAME_RE.match(f[2])
            if status == "C" and mm and mm.group(2) != "0":
                m[rid] = (mm.group(1), mm.group(2), True)
            else:
                m[rid] = ("unclassified", "0", False)
    return m


def main():
    fa = open(os.path.join(EMBED, "reads.fasta"), "w")
    lab = open(os.path.join(EMBED, "labels.tsv"), "w")
    lab.write("uid\tpool\tread_id\ttaxon\ttaxid\tclassified\tlength\n")
    n_total, n_cls = 0, 0
    for out in sorted(glob.glob(os.path.join(KRAKEN, "*.out"))):
        pool = os.path.basename(out)[:-4]
        labels = label_map(out)
        seqs = read_fasta(os.path.join(CLEAN, f"{pool}.clean.fasta"))
        rids = [r for r in seqs if r in labels]
        if len(rids) > CAP:
            rids = random.sample(rids, CAP)
        for rid in rids:
            taxon, taxid, classified = labels[rid]
            seq = seqs[rid]
            uid = f"{pool}|{rid}"
            fa.write(f">{uid}\n{seq}\n")
            lab.write(f"{uid}\t{pool}\t{rid}\t{taxon}\t{taxid}\t"
                      f"{int(classified)}\t{len(seq)}\n")
            n_total += 1
            n_cls += int(classified)
        print(f"{pool}: {len(rids)} reads sampled "
              f"({sum(labels[r][2] for r in rids)} classified)")
    fa.close(); lab.close()
    print(f"TOTAL embedded reads: {n_total} ({n_cls} classified)")


if __name__ == "__main__":
    main()
