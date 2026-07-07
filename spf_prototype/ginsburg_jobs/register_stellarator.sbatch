#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_register_stell
#SBATCH -p short
#SBATCH -t 0-01:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
source /burg/opt/anaconda3-2023.09/etc/profile.d/conda.sh
conda activate spf-stellarator
export PATH="$CONDA_PREFIX/bin:$PATH"
SPF=/ginsburg/astro/users/tjk2147/spf_work; SRC=$SPF/openmc_src; NS=$SPF/spf_pp/native_spf

echo "=== copy StellaratorSource into the tree ==="
cp "$NS/stellarator_source.h"   "$SRC/include/openmc/stellarator_source.h"
cp "$NS/stellarator_source.cpp" "$SRC/src/stellarator_source.cpp"
cp "$NS/stellarator_source.py"  "$SRC/openmc/stellarator_source.py"

echo "=== patch the 3 registration hooks (idempotent, anchor-checked) ==="
python - <<'PY'
SRC="/ginsburg/astro/users/tjk2147/spf_work/openmc_src"
# 1. CMakeLists: new TU next to src/source.cpp
cm=f"{SRC}/CMakeLists.txt"; t=open(cm).read()
if "src/stellarator_source.cpp" not in t:
    assert "  src/source.cpp\n" in t, "CMake anchor missing"
    open(cm,"w").write(t.replace("  src/source.cpp\n","  src/source.cpp\n  src/stellarator_source.cpp\n",1))
    print("CMake patched")
else: print("CMake already patched")
# 2. src/source.cpp: include + Source::create branch
sc=f"{SRC}/src/source.cpp"; t=open(sc).read(); o=t
if "stellarator_source.h" not in t:
    t=t.replace('#include "openmc/source.h"\n','#include "openmc/source.h"\n#include "openmc/stellarator_source.h"\n',1)
if 'source_type == "stellarator"' not in t:
    a='      return make_unique<TokamakSource>(node);\n'
    assert t.count(a)==1, f"cpp anchor count={t.count(a)}"
    t=t.replace(a, a+'    } else if (source_type == "stellarator") {\n      return make_unique<StellaratorSource>(node);\n',1)
if t!=o: open(sc,"w").write(t); print("source.cpp patched")
else: print("source.cpp already patched")
# 3. openmc/source.py: from_xml_element dispatch (local import avoids circularity)
sp=f"{SRC}/openmc/source.py"; t=open(sp).read()
if "source_type == 'stellarator'" not in t:
    a="                return TokamakSource.from_xml_element(elem)\n"
    assert t.count(a)==1, f"py anchor count={t.count(a)}"
    br=("            elif source_type == 'stellarator':\n"
        "                from openmc.stellarator_source import StellaratorSource\n"
        "                return StellaratorSource.from_xml_element(elem)\n")
    open(sp,"w").write(t.replace(a, a+br,1)); print("source.py patched")
else: print("source.py already patched")
PY

echo "=== rebuild (compiles stellarator_source.cpp + relink) ==="
cd "$SRC/build" && cmake .. >/dev/null 2>&1; make -j8 2>&1 | tail -12

echo "=== is StellaratorSource in the built lib? ==="
nm "$SRC/build/lib/libopenmc.so" 2>/dev/null | grep -c StellaratorSource

echo "=== numpy-mirror unit test (cluster) ==="
export PYTHONPATH="$SRC:${PYTHONPATH:-}"; export LD_PRELOAD="$SRC/build/lib/libopenmc.so"
cd "$NS" && python -m pytest test_stellarator_source_cpp.py -q 2>&1 | tail -8
echo DONE_REGISTER
