#!/bin/bash
#SBATCH -A pmg
#SBATCH -p short
#SBATCH --gres=gpu:1
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 60
#SBATCH --exclude=ins082,ins084,ins085
#SBATCH -o /insomnia001/depts/pmg/users/na2933/pivotal-zephyr/embed/embed.log
#SBATCH -e /insomnia001/depts/pmg/users/na2933/pivotal-zephyr/embed/embed.log
export PYTHONUNBUFFERED=1
export NTV3_MODEL=InstaDeepAI/NTv3_8M_pre   # small NTv3 variant: fast load + forward
ENV=/insomnia001/depts/pmg/users/na2933/miniforge/envs/ctx-gene
SP=$ENV/lib/python3.11/site-packages
MODEL=/insomnia001/home/na2933/.cache/huggingface/hub/models--InstaDeepAI--NTv3_8M_pre

# Cold GPFS small-file reads are latency-bound (~30 KB/s single-threaded), so
# prefetch import-critical libs + the model into page cache with 64-way
# parallelism before importing. This collapses ~10 min of serial stalls.
echo "[$(date +%T)] warming page cache (full site-packages + model, parallel)..."
( find "$SP" -type f 2>/dev/null; find "$MODEL" -type f 2>/dev/null ) \
  | xargs -P 96 -n 1 -I{} dd if={} of=/dev/null bs=4M 2>/dev/null
echo "[$(date +%T)] cache warmed; importing + embedding"

source /insomnia001/depts/pmg/users/na2933/miniforge/etc/profile.d/conda.sh
conda activate ctx-gene
cd /insomnia001/depts/pmg/users/na2933/pivotal-zephyr
python -u scripts/embed_reads.py
echo "EMBED_SBATCH_DONE"
