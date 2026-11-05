#!/usr/bin/env python
"""LOCAL integrity closeout: VERIFY (not recompute) that the mass-audit and nested-basis artifacts contain every
required field and that their executable proofs still pass. Reads existing artifacts; re-runs only the cheap proofs
(nesting subset+lstsq; mass spread / swap / not-equal-mass verdict from stored MC volumes). Prints PASS/FAIL per gate."""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SH="/Users/tkiker/lima_share"
ok=True
def check(name, cond, detail=""):
    global ok; ok = ok and cond
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  — {detail}" if detail else ""))

print("=== A. NESTED FOURIER BASIS ===")
man=json.load(open(f"{HERE}/nested_basis_manifest.json"))
lv=man["levels"]; names=["F0","F1","F2","F3","F4"]
check("5 levels F0..F4 present", all(n in lv for n in names), ",".join(lv))
dofs=[lv[n]["n_dof"] for n in names]
check("DOF monotone increasing 5/13/32/59/124", dofs==[5,13,32,59,124], str(dofs))
# required per-level fields
req=["M","N","n_dof","lambda_pol_min_cm","lambda_tor_min_cm","sub_coil_pol","resolves_8coil_tor"]
check("per-level fields complete", all(all(k in lv[n] for k in req) for n in names))
check("nesting_proof recorded", len(man.get("nesting_proof",[]))==4)
check("numeric_nesting_resid < 1e-10", man.get("numeric_nesting_resid",1) < 1e-10, f"{man.get('numeric_nesting_resid'):.1e}")
check("feasibility cap present (m_nyq=9,n_nyq=6)", man["feasibility"]["m_nyquist_pol"]==9 and man["feasibility"]["n_nyquist_tor"]==6)
check("F4 sub-coil poloidal True, 8-coil-tor False", lv["F4"]["sub_coil_pol"] and not lv["F4"]["resolves_8coil_tor"])
# EXECUTABLE PROOF: reload mode arrays, re-verify literal subset + lstsq containment
d=np.load(f"{HERE}/nested_basis_modes.npz")
sets={L:{k.split("__",1)[1] for k in d.files if k.startswith(L+"__")} for L in names}
subset_ok=all(sets[a].issubset(sets[b]) and len(sets[a])<len(sets[b]) for a,b in zip(names[:-1],names[1:]))
check("EXEC nesting literal subset F0⊂F1⊂F2⊂F3⊂F4", subset_ok)
B=np.stack([d[f"F2__{nm}"].ravel() for nm in sorted(sets["F2"])],1)
tgt=d["F1__c_2_2"].ravel(); coef,*_=np.linalg.lstsq(B,tgt,rcond=None); resid=float(np.max(np.abs(B@coef-tgt)))
check("EXEC F1 mode reproduced in F2 basis (lstsq resid<1e-10)", resid<1e-10, f"{resid:.1e}")

print("\n=== B. MASS AUDIT (authoritative 40M MC) ===")
mc=json.load(open(f"{SH}/masses_local.json"))
patches=["t0p2","t0p3","t1p0","t1p1","t1p2","t1p3","t2p0","t2p1","t2p2","t2p3","t3p2","t3p3"]
reqf=["hash","shield_cells","breeder_cells","V_shield","V_shield_err","V_breeder","V_breeder_err","bbox","field_period_fraction"]
check("uniform+zerodelta+12 patches present", all(k in mc for k in ["uniform","zerodelta"]+patches), f"{len(mc)} configs")
check("per-config required fields (hash,vols,errs,bbox,fp)", all(all(f in mc[c] for f in reqf) for c in mc))
check("all DAGMC hashes distinct", len({mc[c]["hash"] for c in mc})==len(mc))
check("comparability: identical bbox + field_period=1 across configs",
      all(np.allclose(mc[c]["bbox"],mc["uniform"]["bbox"],atol=1.0) and mc[c]["field_period_fraction"]==1.0 for c in mc))
# EXECUTABLE PROOF: recompute added mass, spread, swap, not-equal-mass from stored volumes
z=mc["zerodelta"]; RHO_S,RHO_B=15.6,1.94
net={}; swap_ok=True
for p in patches:
    dVs=mc[p]["V_shield"]-z["V_shield"]; dVb=mc[p]["V_breeder"]-z["V_breeder"]
    net[p]=(dVs*RHO_S+dVb*RHO_B)/1000.0
    if abs(dVs+dVb)/abs(dVs) > 0.10: swap_ok=False
lo,hi=min(net.values()),max(net.values())
check("EXEC added-mass spread ~2.45x (163-400 t)", 2.3<hi/lo<2.6 and 150e3<lo<175e3 and 385e3<hi<415e3,
      f"{lo/1e3:.0f}-{hi/1e3:.0f} t, {hi/lo:.2f}x")
check("EXEC fixed-envelope swap dVs≈-dVb (<10%)", swap_ok)
Muni=(mc["uniform"]["V_shield"]-z["V_shield"])*RHO_S/1000+(mc["uniform"]["V_breeder"]-z["V_breeder"])*RHO_B/1000
check("EXEC every patch lighter than uniform control", all(net[p]<Muni for p in patches), f"uniform={Muni/1e3:.0f} t")
check("EXEC t3p3 refutes equal-mass vs uniform (>3σ lighter)", (Muni-net['t3p3'])/ (np.sqrt(2)*15.9e3) > 3,
      f"Δ={(Muni-net['t3p3'])/1e3:.0f} t")

print("\n=== C. SURFACE + FILAMENT-CLEARANCE artifacts ===")
se=np.load(f"{HERE}/surface_extract.npz",allow_pickle=True)
check("surface_extract has loci+areas+layer_names", all(k in se.files for k in ["loci","A_shield","layer_names","V_shield_exact"]))
cd=np.load(f"{HERE}/clearance_domain.npz",allow_pickle=True)
check("clearance_domain has Cfield+avail+domain bounds", all(k in cd.files for k in ["Cfield","avail","I1_delta_min","I2_grow_max"]))
check("filament-clearance min ~29.7cm recorded", 25<float(cd["Cfield"].min())<35, f"{float(cd['Cfield'].min()):.1f} cm")

print("\n"+("ALL CLOSEOUT GATES PASS" if ok else "*** SOME GATES FAILED ***"))
