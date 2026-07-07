#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_fluxmap
#SBATCH -p short
#SBATCH -t 0-02:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
export OMP_NUM_THREADS=8 XLA_FLAGS="--xla_force_host_platform_device_count=8"
PP=/ginsburg/astro/users/tjk2147/spf_work/spf_pp
cd "$PP/python"
echo "=== DESC real sqrt(g) fluxmap: precise_QA example equilibrium (real, converged) ==="
"$HOME/desc_venv/bin/python" -c "import quasr_fluxmap as q; q.build(0, example='precise_QA', n_rho=12, n_theta=48, n_zeta=144)"
echo DONE_FLUXMAP
