"""Precompute magnet (and breeder/first_wall) OpenMC cell IDs for a DAGMC h5m via pymoab,
so coil_run_v3 can CellFilter the coils without openmc.lib c.fill (NotImplementedError on DAGMC).
cell_id = MOAB volume GLOBAL_ID + 1 (DAGMC auto_geom_ids convention). Usage: python step1_cells.py H5M OUT.npz"""
import sys
import numpy as np
from pymoab import core, types

H5M = sys.argv[1]
OUT = sys.argv[2]
mb = core.Core(); mb.load_file(H5M)
name_tag = mb.tag_get_handle(types.NAME_TAG_NAME)
cat_tag  = mb.tag_get_handle(types.CATEGORY_TAG_NAME)
gid_tag  = mb.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
root = mb.get_root_set()
groups = mb.get_entities_by_type_and_tag(root, types.MBENTITYSET,
            np.array([cat_tag]), np.array(["Group"]))

def cells_for(keyword):
    ids = []
    for g in groups:
        try: nm = mb.tag_get_data(name_tag, g, flat=True)[0]
        except Exception: nm = ""
        nm = nm.strip("\x00") if isinstance(nm, str) else str(nm)
        if keyword in nm.lower():
            for v in mb.get_entities_by_type(g, types.MBENTITYSET):
                try: ids.append(int(mb.tag_get_data(gid_tag, v, flat=True)[0]) + 1)
                except Exception: pass
    return sorted(set(ids))

mag = cells_for("magnet"); brd = cells_for("breeder"); fw = cells_for("first_wall")
print(f"{H5M}", flush=True)
print(f"  magnets   n={len(mag):3d} ids={mag}", flush=True)
print(f"  breeder   n={len(brd):3d} ids={brd}", flush=True)
print(f"  first_wall n={len(fw):3d} ids={fw}", flush=True)
np.savez(OUT, magnet_cells=np.array(mag), breeder_cells=np.array(brd), first_wall_cells=np.array(fw), h5m=H5M)
print(f"WROTE {OUT}", flush=True)
