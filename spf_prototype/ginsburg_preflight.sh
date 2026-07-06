#!/bin/bash
# Ginsburg LOGIN-NODE preflight checks for the conformal stellarator SPF run.
# 100% READ-ONLY: it inspects storage, conda/modules, cross sections, SLURM, and
# build tools so we can set the real paths in sweep/sweep.sbatch (see sweep/docs/README.md)
# BEFORE cloning + creating the (multi-GB) conda env. Nothing is created or changed.
#
# Usage on the login node:
#     bash ginsburg_preflight.sh 2>&1 | tee ginsburg_preflight.out
# then paste ginsburg_preflight.out back.

echo "===== A. STORAGE / WHERE TO WORK ====="
echo "USER=$USER  HOME=$HOME  PWD=$PWD"
checkquota 2>/dev/null || mmlsquota 2>/dev/null || quota -s 2>/dev/null || echo "(no quota tool found - try 'gpfsquota')"
ls -ld /burg/astro /burg/astro/users/$USER /burg/home/$USER 2>/dev/null
df -h /burg/home/$USER /burg/astro 2>/dev/null
ls -ld /burg/scratch /scratch /burg/astro/scratch 2>/dev/null || echo "(no obvious scratch mount)"

echo; echo "===== B. CONDA / MODULES ====="
which conda; conda --version; echo "base prefix: ${CONDA_PREFIX:-<none>}"
module avail anaconda 2>&1 | grep -i anaconda | head
conda config --show channels 2>/dev/null

echo; echo "===== C. CROSS SECTIONS (a shared copy saves a multi-GB download) ====="
echo "OPENMC_CROSS_SECTIONS=${OPENMC_CROSS_SECTIONS:-<unset>}"
find /burg -maxdepth 4 -iname "cross_sections.xml" 2>/dev/null | head
ls -d /burg/opt/*endf* /burg/opt/*openmc* /burg/opt/*nndc* 2>/dev/null || echo "(no shared XS under /burg/opt)"

echo; echo "===== D. SLURM (account / partition / walltime) ====="
sinfo -s 2>/dev/null | head
sacctmgr -n show assoc user=$USER format=account,partition,qos%40 2>/dev/null
scontrol show partition 2>/dev/null | grep -E "PartitionName|Default=|MaxTime" | head -20

echo; echo "===== E. BUILD TOOLS + GIT ====="
which git cmake g++ gcc 2>/dev/null
cmake --version 2>/dev/null | head -1
g++ --version 2>/dev/null | head -1

echo; echo "===== F. (optional) does DAGMC-enabled OpenMC solve on conda-forge here? ====="
echo "Run this MANUALLY on the login node (has internet); it is a dry-run, installs nothing:"
echo "  conda create -n _solvetest --dry-run -c conda-forge \"openmc=0.15.*=dagmc*\" 2>&1 | tail -20"

echo; echo "===== preflight done ====="
