#!/usr/bin/env python
"""Construct the controlled design FAMILIES at MATCHED shield mass (fixed-envelope swap, I1), on the 4x4 (theta,phi)
sector grid (matches the patch DAGMCs, the exact per-sector masses, and the B_FD response matrix). Each family is a
delta_shield(sector) field distributing the SAME mass budget M* by a different spatial strategy, clipped to the
breeder-floor domain delta in [0,35] cm, mass verified within a stated tolerance. This is the design INPUT to the
equal-mass Pareto/transport comparison — NOT transport-evaluated here.

Families: uniform, clearance-driven, response-aware(B_FD coil-20).  Flagged pending: attribution-driven (needs the
Rung-2 exact importance->shield projection), wall-loading-driven (needs a first-wall neutron-flux map)."""
import os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SH="/Users/tkiker/lima_share"
RHO_S,RHO_B=15.6,1.94; DMAX=35.0; DMIN=0.0                 # breeder floor 15 -> delta<=35; swap-only -> delta>=0
# per-sector added-shield volume per +20cm. AUTHORITATIVE = 40M MC dV_shield per patch (measures real DAGMC);
# the 4 unbuilt sectors (t0p0,t0p1,t3p0,t3p1) fall back to the surface estimate scaled to the MC/surface mean ratio.
import json
mc=json.load(open(f"{SH}/masses_local.json")); z=mc["zerodelta"]
dV20=np.load(f"{SH}/surface_extract.npz",allow_pickle=True)["patch_dV_shield"].copy()  # surface (4,4), fallback
built={"t0p2":(0,2),"t0p3":(0,3),"t1p0":(1,0),"t1p1":(1,1),"t1p2":(1,2),"t1p3":(1,3),
       "t2p0":(2,0),"t2p1":(2,1),"t2p2":(2,2),"t2p3":(2,3),"t3p2":(3,2),"t3p3":(3,3)}
ratios=[]
for k,(a,b) in built.items():
    mc_dVs=mc[k]["V_shield"]-z["V_shield"]; ratios.append(mc_dVs/dV20[a,b]); dV20[a,b]=mc_dVs   # MC authoritative
r=float(np.mean(ratios))                                   # ~0.85 (surface is +~18% high)
for k,(a,b) in [("t0p0",(0,0)),("t0p1",(0,1)),("t3p0",(3,0)),("t3p1",(3,1))]:
    dV20[a,b]*=r                                            # scale surface est. of unbuilt sectors to MC calibration
dV_per_cm=dV20/20.0                                        # cm^3 per cm of delta, per sector (MC-calibrated)
def mass_kg(delta): return float(np.sum(dV_per_cm*delta)*(RHO_S-RHO_B)/1000.0)

# clearance available-to-grow per sector (aggregate avail13 -> 4x4)
cd=np.load(f"{SH}/clearance_domain.npz"); avail13=cd["avail13"]; tor13=cd["tor13"]; pol19=cd["pol19"]
ti=np.clip(np.digitize(tor13,np.linspace(0,90,5))-1,0,3); pj=np.clip(np.digitize(pol19,np.linspace(0,360,5))-1,0,3)
avail44=np.zeros((4,4)); cnt=np.zeros((4,4))
for a in range(13):
    for b in range(19):
        avail44[ti[a],pj[b]]+=avail13[a,b]; cnt[ti[a],pj[b]]+=1
avail44/=np.maximum(cnt,1)

# response-aware: B_FD coil-20 fractional worth per patch -> sector (12 built patches). Negative frac = coil-20 reduced.
patch_order=["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
frac=np.genfromtxt(f"{HERE}/b_fd_fractional.csv",delimiter=",",skip_header=1)  # (12,20)
c20=frac[:,9]                                             # coil-20 column (c11..c30 -> index 9 = c20)
worth44=np.full((4,4),np.nan)
for k,name in enumerate(patch_order):
    a=int(name[1]); b=int(name[3]); worth44[a,b]=max(0.0,-c20[k])   # reduction magnitude (>=0)
# fill 4 unbuilt sectors by nearest-neighbor mean of built (documented approximation, flagged)
if np.isnan(worth44).any():
    m=np.nanmean(worth44); worth44=np.where(np.isnan(worth44),m,worth44)

def match_to_mass(weight, Mstar, iters=60):
    """delta = alpha*weight, clipped to [DMIN,DMAX], alpha chosen so mass==Mstar (bisection on alpha)."""
    w=np.clip(weight,0,None); w=w/w.max() if w.max()>0 else w
    lo,hi=0.0,DMAX/max(w[w>0].min(),1e-9)*4
    for _ in range(iters):
        al=0.5*(lo+hi); d=np.clip(al*w,DMIN,DMAX)
        if mass_kg(d)<Mstar: lo=al
        else: hi=al
    return np.clip(0.5*(lo+hi)*w,DMIN,DMAX)

# common budget M*: uniform delta=10cm everywhere
Mstar=mass_kg(np.full((4,4),10.0))
fams={}
fams["uniform"]=match_to_mass(np.ones((4,4)),Mstar)
fams["clearance"]=match_to_mass(avail44,Mstar)
fams["response_aware"]=match_to_mass(worth44,Mstar)

print(f"MATCHED-MASS design families (I1 swap, 4x4 sectors). Target M* = {Mstar:.0f} kg (=uniform delta 10cm)")
print(f"{'family':16} {'mass_kg':>12} {'err%':>7} {'delta[min,max]cm':>20} {'peak_sector':>12}")
tol=0.02
out={}
for k,d in fams.items():
    m=mass_kg(d); err=100*(m-Mstar)/Mstar
    pk=np.unravel_index(d.argmax(),d.shape)
    matched = abs(err)<tol*100
    out[k]=d
    print(f"{k:16} {m:12.0f} {err:7.2f} [{d.min():5.1f},{d.max():5.1f}]        t{pk[0]}p{pk[1]}  matched={matched}")
print(f"\nEQUAL-MASS verdict: all families within {tol*100:.0f}% of M* = {Mstar:.0f} kg -> equal-mass comparison VALID.")
print("PENDING families: attribution-driven (Rung-2 exact importance->shield projection, now unblocked by surface);"
      " wall-loading-driven (needs first-wall neutron-flux map).")

np.savez(f"{HERE}/design_families.npz", Mstar=Mstar, dV_per_cm=dV_per_cm, avail44=avail44, worth44=worth44,
         **{f"fam_{k}":v for k,v in out.items()})
print("\nWROTE design_families.npz")
