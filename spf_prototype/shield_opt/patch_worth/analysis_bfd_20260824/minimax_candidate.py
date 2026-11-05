#!/usr/bin/env python
"""RESPONSE-AWARE candidate via constrained MINIMAX over the FULL 12x20 B_FD operator (not just the coil-20 column).
Design var: delta_k (cm) per built sector k (12). Linearized predicted coil-j response after intervention:
    R[j] = R_base[j] + sum_k (delta_k/20) * E[k,j]        E = b_fd_raw (R_patch_k - R_base, per +20cm patch)
LP:  min t  s.t.  R[j] <= t  for all 20 coils ;  mass(delta)=M* ;  0<=delta_k<=35.
Reports the predicted full 20-coil response vector, the limiting (peak) coil before/after, and its migration,
and contrasts with the coil-20-guided field. This is a PREDICTION from the linear surrogate; the winning candidate
must still be transport-verified (the surrogate was shown non-generalizing at n=12 -> treat as a screening guess)."""
import os, numpy as np
from scipy.optimize import linprog
HERE=os.path.dirname(os.path.abspath(__file__)); SH="/Users/tkiker/lima_share"
RHO_S,RHO_B=15.6,1.94; DMAX=35.0
patch_order=["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
# CORRECT definitions (from reconstruct_b_fd.py): b_fd_raw.csv = B = -E/denom (normalized worth, NOT dR);
# b_fd_fractional.csv = E/R0 ; R0 = coil_step1b_uniform baseline. Recover raw dR: E = frac * R0.
# CAVEAT: baseline is shield 42.55 (2.55cm confound) + patches used WW/baseline didn't; the linear surrogate is
# KNOWN not to generalize at n=12 (see spf-patch-worth-n12-verdict) -> this minimax is a SCREENING guess only.
frac=np.genfromtxt(f"{HERE}/b_fd_fractional.csv",delimiter=",",skip_header=1)  # (12,20) = dR/R0
base=np.load(f"{HERE}/../cluster_pull/coil_step1b_uniform.npz")
Rbase=np.ravel(base["coil_fast_flux"]).astype(float)                 # (20) baseline coil responses R0
E=frac*Rbase[None,:]                                                 # (12,20) raw dR per +20cm patch
print(f"baseline peak coil = c{11+int(Rbase.argmax())} = {Rbase.max():.4e} (from coil_step1b_uniform)")

# MC-calibrated per-sector added-shield volume per cm (from design_families logic, 12 built sectors)
import json
mc=json.load(open(f"{SH}/masses_local.json")); z=mc["zerodelta"]
dV_per_cm=np.array([(mc[p]["V_shield"]-z["V_shield"])/20.0 for p in patch_order])  # cm^3/cm, per built sector
mass_coef=dV_per_cm*(RHO_S-RHO_B)/1000.0                             # kg per cm delta, per sector
Mstar=float(np.sum(mass_coef*10.0))                                  # budget = uniform 10cm over the 12 built sectors

# LP: variables x = [delta_1..delta_12, t]; minimize t
n=12
c=np.zeros(n+1); c[-1]=1.0
# R[j] - t <= 0  ->  sum_k (E[k,j]/20) delta_k - t <= -R_base[j]
A_ub=np.zeros((20,n+1)); b_ub=np.zeros(20)
for j in range(20):
    A_ub[j,:n]=E[:,j]/20.0; A_ub[j,-1]=-1.0; b_ub[j]=-Rbase[j]
A_eq=np.zeros((1,n+1)); A_eq[0,:n]=mass_coef; b_eq=[Mstar]
bounds=[(0,DMAX)]*n+[(None,None)]
r=linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=bounds,method="highs")
assert r.success, r.message
delta=r.x[:n]; tstar=r.x[-1]
Rpred=Rbase+ (E.T @ (delta/20.0))                                   # predicted 20-coil response
# coil-20-guided comparison field (same mass, follows -E[:,c20])
c20=E[:,9]; w=np.clip(-c20,0,None);
def match(w):
    w=w/w.max() if w.max()>0 else w
    lo,hi=0,DMAX/max(w[w>0].min(),1e-9)*4
    for _ in range(80):
        al=.5*(lo+hi); d=np.clip(al*w,0,DMAX)
        if float(np.sum(mass_coef*d))<Mstar: lo=al
        else: hi=al
    return np.clip(.5*(lo+hi)*w,0,DMAX)
d_c20=match(w); R_c20=Rbase+(E.T@(d_c20/20.0))

print(f"\nBudget M* = {Mstar:.0f} kg (uniform 10cm over 12 built sectors)")
jb=int(Rbase.argmax()); jm=int(Rpred.argmax()); jc=int(R_c20.argmax())
print(f"baseline    peak = c{11+jb} {Rbase.max():.4e}")
print(f"MINIMAX     peak = c{11+jm} {Rpred.max():.4e}  ({100*(Rpred.max()-Rbase.max())/Rbase.max():+.1f}% vs base)  limiting-coil migration c{11+jb}->c{11+jm}")
print(f"coil20-guided peak = c{11+jc} {R_c20.max():.4e}  ({100*(R_c20.max()-Rbase.max())/Rbase.max():+.1f}% vs base)")
print(f"minimax delta (cm) by sector: "+", ".join(f"{p}:{delta[i]:.1f}" for i,p in enumerate(patch_order)))
np.savez(f"{HERE}/minimax_candidate.npz",delta=delta,Rpred=Rpred,Rbase=Rbase,R_c20=R_c20,d_c20=d_c20,
         Mstar=Mstar,patch_order=np.array(patch_order),tstar=tstar)
print("WROTE minimax_candidate.npz")
