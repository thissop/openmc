# Source this INSIDE the Lima 'spf' VM to get a working OpenMC+DAGMC environment.
#   limactl shell spf
#   source /Users/tkiker/Documents/GitHub/openmc/spf_prototype/local_worker/activate_worker.sh
#
# Built 2026-07-26 on the M2 (arm64) Lima VM 'spf'. See local_worker/README.md.
source ~/miniforge3/etc/profile.d/conda.sh
conda activate spf
export PATH=$HOME/opt/openmc/bin:$PATH
export LD_LIBRARY_PATH=$HOME/opt/openmc/lib:$HOME/opt/dagmc/lib:$CONDA_PREFIX/lib
export OPENMC_CROSS_SECTIONS=$HOME/xs/endfb80/cross_sections.xml
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
echo "[spf worker] openmc $(openmc --version 2>/dev/null | head -1 | awk '{print $3}') | DAGMC $(openmc --version 2>/dev/null | grep -i 'DAGMC support' | awk '{print $3}') | XS=$OPENMC_CROSS_SECTIONS"
