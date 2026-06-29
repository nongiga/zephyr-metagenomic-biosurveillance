#!/usr/bin/env python3
"""Per-read reference label = minimap2 best hit (via mappy), against the curated
respiratory-virus panel. More trustworthy than Kraken for divergent ONT reads
(it labels reads Kraken left unclassified) and it resolves the RV-A/B/C
cross-mapping by taking the single best-scoring target. Reads with no alignment
-> 'unmapped'. Writes embed/besthit.tsv."""
import glob
import os
import mappy

ROOT = os.environ.get("PZ_ROOT", "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr")
REFS = os.path.join(ROOT, "refs")
EMBED = os.path.join(ROOT, "embed")


def build_named_refs(path):
    """Concatenate refs with header = virus name (file stem)."""
    out = os.path.join(EMBED, "refs_named.fa")
    with open(out, "w") as o:
        for fa in sorted(glob.glob(os.path.join(REFS, "*.fa"))):
            name = os.path.basename(fa)[:-3]
            seq = "".join(l.strip() for l in open(fa) if not l.startswith(">"))
            o.write(f">{name}\n{seq}\n")
    return out


def read_fasta(path):
    uid, buf = None, []
    for line in open(path):
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
    ref = build_named_refs(REFS)
    aln = mappy.Aligner(ref, preset="map-ont")
    if not aln:
        raise SystemExit("failed to build mappy index")
    n_map = 0
    with open(os.path.join(EMBED, "besthit.tsv"), "w") as out:
        out.write("uid\tbesthit\tmatch_len\n")
        for uid, seq in read_fasta(os.path.join(EMBED, "reads.fasta")):
            best = None
            for h in aln.map(seq):  # primary + secondary; keep best by matches
                if best is None or h.mlen > best.mlen:
                    best = h
            if best is not None:
                out.write(f"{uid}\t{best.ctg}\t{best.mlen}\n")
                n_map += 1
            else:
                out.write(f"{uid}\tunmapped\t0\n")
    print(f"best-hit labels written; {n_map} mapped of reads")


if __name__ == "__main__":
    main()
