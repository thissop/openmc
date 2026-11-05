#!/usr/bin/env python
"""Expand each matched-mass design family (4x4 sector delta, cm) onto build_delta.py's (13,19) toroidal/poloidal
grid, save one npz per family with key 'delta' (shape (13,19)) ready for `build_delta.py <fam.npz> <name>`.

Sector mapping is the EXACT inverse of design_families.py's aggregation: for grid point (a,b), delta = fam[ti[a],pj[b]]
where ti,pj digitize (tor13,pol19) into the 4x4 sectors. Hard-gates build_delta's fixed-envelope invariants
(0 <= delta <= 35 -> breeder in [15,50], envelope conserved) so a bad field fails HERE, not after a cluster build.

Families: uniform / clearance / response_aware (from design_families.npz) + minimax (minimax_candidate.npz,
12 built sectors; the 4 unbuilt corner sectors t0p0,t0p1,t3p0,t3p1 are set to 0 -- flagged, they are not in the
12x20 operator). Each family carries the SAME target mass M*, so peak-coil differences are a pure placement effect."""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SH = "/Users/tkiker/lima_share"
OUT = f"{HERE}/family_grids"; os.makedirs(OUT, exist_ok=True)

cd = np.load(f"{HERE}/clearance_domain.npz")
tor13, pol19 = cd["tor13"], cd["pol19"]                       # (13,), (19,) grid node angles (deg)
ti = np.clip(np.digitize(tor13, np.linspace(0, 90, 5)) - 1, 0, 3)   # 13 -> sector row 0..3  (matches design_families)
pj = np.clip(np.digitize(pol19, np.linspace(0, 360, 5)) - 1, 0, 3)  # 19 -> sector col 0..3

def expand(fam44):
    g = np.zeros((13, 19))
    for a in range(13):
        for b in range(19):
            g[a, b] = fam44[ti[a], pj[b]]
    return g

fams = {}
df = np.load(f"{HERE}/design_families.npz")
for k, name in [("fam_uniform", "uniform"), ("fam_clearance", "clearance"), ("fam_response_aware", "response")]:
    fams[name] = df[k].reshape(4, 4)

# minimax: 12 built-sector deltas -> (4,4); unbuilt corners = 0
mm = np.load(f"{HERE}/minimax_candidate_matched.npz", allow_pickle=True)  # re-matched to the common M* (see minimax_rematch.py)
patch_order = list(mm["patch_order"]); mmd = mm["delta"]
mfield = np.zeros((4, 4))
for val, p in zip(mmd, patch_order):
    p = p if isinstance(p, str) else p.decode()
    mfield[int(p[1]), int(p[3])] = val
fams["minimax"] = mfield

# MC-calibrated per-sector shield volume per cm (for a mass audit at expansion time, same coefs as design_families)
import json
mc = json.load(open(f"{SH}/masses_local.json")); z = mc["zerodelta"]
dV20 = np.load(f"{SH}/surface_extract.npz", allow_pickle=True)["patch_dV_shield"].copy()
built = {"t0p2":(0,2),"t0p3":(0,3),"t1p0":(1,0),"t1p1":(1,1),"t1p2":(1,2),"t1p3":(1,3),
         "t2p0":(2,0),"t2p1":(2,1),"t2p2":(2,2),"t2p3":(2,3),"t3p2":(3,2),"t3p3":(3,3)}
ratios = []
for kk, (a, b) in built.items():
    r_ = (mc[kk]["V_shield"] - z["V_shield"]) / dV20[a, b]; ratios.append(r_); dV20[a, b] = mc[kk]["V_shield"] - z["V_shield"]
r = float(np.mean(ratios))
for kk, (a, b) in [("t0p0",(0,0)),("t0p1",(0,1)),("t3p0",(3,0)),("t3p1",(3,1))]:
    dV20[a, b] *= r
dV_per_cm = dV20 / 20.0
def mass_t(fam44): return float(np.sum(dV_per_cm * fam44) * (15.6 - 1.94) / 1000.0) / 1e3  # tonnes

print(f"{'family':12} {'mass_t':>8} {'delta[min,max]':>16} {'grid[min,max]':>16}  invariants")
for name, fam in fams.items():
    g = expand(fam)
    assert g.shape == (13, 19)
    assert np.all(g >= -1e-9), f"{name}: delta<0 -> build_delta rejects"
    assert np.all(g <= 35 + 1e-9), f"{name}: delta>35 -> breeder<15 floor violated"
    g = np.clip(g, 0.0, 35.0)
    np.savez(f"{OUT}/fam_{name}.npz", delta=g, fam44=fam, tor13=tor13, pol19=pol19)
    print(f"{name:12} {mass_t(fam):8.0f} [{fam.min():5.1f},{fam.max():5.1f}]   [{g.min():5.1f},{g.max():5.1f}]   OK")
print(f"\nWROTE {OUT}/fam_{{uniform,clearance,response,minimax}}.npz  (key 'delta', shape (13,19))")
print("Cluster: for f in uniform clearance response minimax; do python build_delta.py family_grids/fam_$f.npz fam_$f; done")
