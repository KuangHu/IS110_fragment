#!/bin/bash
#SBATCH --job-name=is_sizes
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/is_sizes_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/is_sizes_%j.err

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"
python3 /tmp/compute_is_sizes.py
