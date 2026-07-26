#!/usr/bin/env python
"""Map DAGMC coil volumes (openmc cell IDs 8-27) to toroidal/poloidal angles.

Uses MOAB directly on dagmc_baseline.h5m. DAGMC volume GLOBAL_ID == openmc cell id
(auto_geom_ids, no collisions: in-vessel volumes 3-7, coils 8-27). We compute each
coil volume centroid from its triangle vertices, then phi=atan2(y,x) (toroidal) and
a poloidal angle about the local magnetic axis ring.
"""
import numpy as np
from pymoab import core, types, rng

H5M = "/burg-archive/home/tjk2147/pstl_test/dagmc_baseline.h5m"

mb = core.Core()
mb.load_file(H5M)

# category tag -> find Volumes and Groups
cat = mb.tag_get_handle(types.CATEGORY_TAG_NAME)
gid = mb.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
name_tag = mb.tag_get_handle(types.NAME_TAG_NAME)

root = mb.get_root_set()

# --- group (material) membership: which global-ids are 'mat:magnets' ---
groups = mb.get_entities_by_type_and_tag(
    root, types.MBENTITYSET, np.array([cat]), np.array(["Group"]))
magnet_gids = set()
for g in groups:
    try:
        nm = mb.tag_get_data(name_tag, g, flat=True)[0]
    except Exception:
        nm = ""
    if isinstance(nm, bytes):
        nm = nm.decode(errors="ignore")
    if "magnet" in str(nm).lower():
        # volumes in this group
        vols = mb.get_entities_by_type_and_tag(
            g, types.MBENTITYSET, np.array([cat]), np.array(["Volume"]))
        for v in vols:
            magnet_gids.add(int(mb.tag_get_data(gid, v, flat=True)[0]))
print("magnet group global-ids:", sorted(magnet_gids))

# --- per-volume centroid from triangle vertices ---
vols = mb.get_entities_by_type_and_tag(
    root, types.MBENTITYSET, np.array([cat]), np.array(["Volume"]))

records = []
for v in vols:
    vid = int(mb.tag_get_data(gid, v, flat=True)[0])
    if vid not in magnet_gids:
        continue
    # volume -> child surfaces -> triangles -> vertices
    surfs = mb.get_child_meshsets(v)
    all_coords = []
    for s in surfs:
        tris = mb.get_entities_by_type(s, types.MBTRI)
        if len(tris) == 0:
            continue
        verts = mb.get_connectivity(tris)
        all_coords.append(mb.get_coords(verts).reshape(-1, 3))
    if not all_coords:
        continue
    coords = np.vstack(all_coords)
    c = coords.mean(axis=0)  # cm
    # openmc cell id = MOAB GLOBAL_ID + 1 (auto_geom_ids offset)
    records.append((vid + 1, c[0], c[1], c[2]))

records.sort()
vids = np.array([r[0] for r in records])
cx = np.array([r[1] for r in records])
cy = np.array([r[2] for r in records])
cz = np.array([r[3] for r in records])
phi = np.degrees(np.arctan2(cy, cx)) % 360.0
R = np.hypot(cx, cy)

print("\n cellID    X(cm)     Y(cm)     Z(cm)     R(cm)   phi(deg)")
for i in range(len(vids)):
    print(f"  {vids[i]:4d}  {cx[i]:8.1f} {cy[i]:8.1f} {cz[i]:8.1f} {R[i]:8.1f}  {phi[i]:7.2f}")

np.savez("/burg-archive/home/tjk2147/pstl_test/ww_pert/coil_centroids.npz",
         cell_ids=vids, cx=cx, cy=cy, cz=cz, R=R, phi=phi)
print("\nSAVED coil_centroids.npz")
