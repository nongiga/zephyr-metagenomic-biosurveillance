#!/usr/bin/env python3
"""Embed cleaned viral reads with NTv3-650M, reusing the ctx-gene harness.

Whole-read embedding = NTv3Evaluator.embed_gene(seq, 0, len(seq)) with layer=None
(mean-pooled full-resolution last hidden state). Runs offline against the cached
InstaDeepAI/NTv3_650M_pre. GPU required.
"""
import os
import sys
import time

os.environ.setdefault("HF_HOME", "/insomnia001/home/na2933/.cache/huggingface")
os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
sys.path.insert(0, "/insomnia001/depts/pmg/users/na2933/ctx-gene")

import numpy as np
from evaluators.ntv3 import NTv3Evaluator

ROOT = "/insomnia001/depts/pmg/users/na2933/pivotal-zephyr"
EMBED = os.path.join(ROOT, "embed")
MODEL_NAME = os.environ.get("NTV3_MODEL", "InstaDeepAI/NTv3_8M_pre")
MAX_LEN = 20000  # NTv3 context ceiling; rare longer reads are truncated


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
    reads = list(read_fasta(os.path.join(EMBED, "reads.fasta")))
    print(f"loaded {len(reads)} reads; model = {MODEL_NAME}")
    ev = NTv3Evaluator(model_name=MODEL_NAME)
    print("model loaded")

    embs, ids = [], []
    n_trunc = 0
    t0 = time.time()
    for i, (uid, seq) in enumerate(reads):
        seq = seq.upper()
        if len(seq) > MAX_LEN:
            seq = seq[:MAX_LEN]
            n_trunc += 1
        emb = ev.embed_gene(seq, 0, len(seq)).float().numpy()
        embs.append(emb)
        ids.append(uid)
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(reads)}  ({time.time()-t0:.0f}s)")

    X = np.vstack(embs).astype(np.float32)
    np.save(os.path.join(EMBED, "embeddings.npy"), X)
    with open(os.path.join(EMBED, "ids.txt"), "w") as fh:
        fh.write("\n".join(ids) + "\n")
    print(f"saved embeddings {X.shape}; truncated {n_trunc} long reads; "
          f"{time.time()-t0:.0f}s total")


if __name__ == "__main__":
    main()
