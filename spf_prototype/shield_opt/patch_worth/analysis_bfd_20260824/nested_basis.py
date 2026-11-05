#!/usr/bin/env python
"""CORRECTED basis: ONE genuinely nested stellarator-symmetric 2-D Fourier sequence F0 ⊂ F1 ⊂ F2 ⊂ F3 (⊂ F4)
for the shield-thickness field t(theta_pol, phi_tor) on the ParaStell (NT=13 toroidal x NPOL=19 poloidal) grid, nfp=4.

Replaces the earlier INVALID 20/12/35 A/B/C sequence, which mixed basis FAMILIES (indicator sectors vs Fourier) and
was NOT nested. Here every level is the same family (2-D Fourier) truncated at increasing order, so nesting is exact.

Stellarator symmetry t(theta,phi)=t(-theta,-phi)  ->  keep cosines  cos(m*theta - n*nfp*phi)  (VMEC/ParaStell convention,
n ranges +/- for m>=1 to carry both helicities; m=0 keeps n>=0 only). Index set for truncation (M,N):
    S(M,N) = {(0,n): 0<=n<=N}  U  {(m,n): 1<=m<=M, -N<=n<=N}
Because M_K, N_K increase monotonically, S(M_K,N_K) ⊂ S(M_{K+1},N_{K+1})  ==>  span(F_K) ⊂ span(F_{K+1}). PROVEN below
by asserting each level's mode-name set is a strict subset of the next (literal containment, not numerical).

Physical scales (exact device constants): R0=1367 cm, a_minor=217 cm, nfp=4.
    poloidal:  lambda_pol(m) = 2*pi*a / m
    toroidal:  lambda_tor(n) = (2*pi*R0/nfp) / n     (n oscillations per field period)
Coil footprint (given): ~227 cm poloidal; 8 coils/period -> toroidal coil pitch (2*pi*R0/nfp)/8 = 268 cm.
GRID NYQUIST (hard feasibility cap of the 13x19 surface): poloidal m_max = (NPOL-1)//2 = 9; toroidal n_max = (NT-1)//2 = 6.
"""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
R0, A, NFP = 1367.0, 217.0, 4
NT, NPOL = 13, 19
tor=np.linspace(0,90,NT); pol=np.linspace(0,360,NPOL)
TH=np.radians(pol)[None,:]*np.ones((NT,1)); PH=np.radians(tor)[:,None]*np.ones((1,NPOL))
M_NYQ, N_NYQ = (NPOL-1)//2, (NT-1)//2                       # 9, 6
COIL_POL, COIL_TOR = 227.0, (2*np.pi*R0/NFP)/8              # 227, 268 cm
lam_pol=lambda m:(2*np.pi*A/m) if m else np.inf
lam_tor=lambda n:((2*np.pi*R0/NFP)/abs(n)) if n else np.inf

def index_set(M,N):
    idx=[(0,n) for n in range(0,N+1)]
    for m in range(1,M+1):
        for n in range(-N,N+1): idx.append((m,n))
    return idx
def mode_array(m,n): return np.cos(m*TH - n*NFP*PH)
def modes(M,N): return {f"c_{m}_{n}":mode_array(m,n) for (m,n) in index_set(M,N)}

# nested truncation ladder (M,N) monotone increasing -> nesting guaranteed
LADDER=[("F0",1,1),("F1",2,2),("F2",4,3),("F3",6,4),("F4",M_NYQ,N_NYQ)]  # F4 = grid Nyquist cap
levels={}; mode_sets={}; all_arrays={}
for name,M,N in LADDER:
    idx=index_set(M,N); md=modes(M,N); mode_sets[name]=set(md)
    m_hi=M; n_hi=N
    levels[name]={"M":M,"N":N,"n_dof":len(idx),
        "lambda_pol_min_cm":round(lam_pol(m_hi),0),"lambda_tor_min_cm":round(lam_tor(n_hi),0),
        "finest_vs_coil_pol":f"{lam_pol(m_hi)/COIL_POL:.2f}x coil", "finest_vs_coil_tor":f"{lam_tor(n_hi)/COIL_TOR:.2f}x coil",
        "sub_coil_pol": bool(lam_pol(m_hi)<COIL_POL), "resolves_8coil_tor": bool(n_hi>=8)}
    all_arrays.update({f"{name}__{k}":v for k,v in md.items()})

# PROOF OF NESTING: each level's mode-name set is a strict subset of the next
nesting_ok=True; proof=[]
order=[n for n,_,_ in LADDER]
for a,b in zip(order[:-1],order[1:]):
    sub=mode_sets[a].issubset(mode_sets[b]); strict=sub and len(mode_sets[a])<len(mode_sets[b])
    proof.append(f"{a} subset {b}: {sub} (|{a}|={len(mode_sets[a])} < |{b}|={len(mode_sets[b])} strict={strict})")
    nesting_ok &= strict
assert nesting_ok, "NESTING PROOF FAILED"

# numeric confirmation: least-squares fit an F_K mode with the F_{K+1} basis -> exact (residual ~ 0)
def basis_matrix(level,names): return np.stack([all_arrays[f"{level}__{nm}"].ravel() for nm in names],1)
Bnext=basis_matrix("F2",sorted(mode_sets["F2"])); tgt=all_arrays["F1__c_2_2"].ravel()
coef,res,rk,sv=np.linalg.lstsq(Bnext,tgt,rcond=None); recon=Bnext@coef
num_resid=float(np.max(np.abs(recon-tgt)))
assert num_resid<1e-10, f"numeric nesting check failed resid={num_resid}"

feas={"grid":[NT,NPOL],"m_nyquist_pol":M_NYQ,"n_nyquist_tor":N_NYQ,
 "coil_pol_cm":COIL_POL,"coil_tor_pitch_cm":round(COIL_TOR,0),
 "sub_coil_needs_m_gt":6,"finest_resolvable_lambda_pol_cm":round(lam_pol(M_NYQ),0),
 "toroidal_8coil_needs_n":8,"grid_max_n":N_NYQ,
 "FEASIBILITY_LIMIT":("On the current 13x19 surface the poloidal grid resolves only up to m=9 (lambda_pol~151cm, "
   "marginally sub-coil), and the toroidal grid resolves only up to n=6 (lambda_tor~358cm) -- it CANNOT resolve the "
   "8-coil-per-period toroidal structure, which needs n~8 (lambda_tor~268cm). A genuinely sub-coil-in-BOTH design "
   "study therefore requires re-meshing the ParaStell surface finer (>=17 toroidal, >=21 poloidal); this is the exact "
   "geometric feasibility cap, not a modeling choice.")}

robustness={"separate_families_NOT_in_Fourier_hierarchy":
   ["sector_indicators (16 coarse patches) -- interpretable, discontinuous; use as a basis-FAMILY robustness check",
    "B-splines (smooth compact support) -- a LATER separate family for the same convergence test, not a continuation of F0..F4"]}

np.savez(f"{HERE}/nested_basis_modes.npz", **all_arrays, tor=tor, pol=pol, R0=R0, a=A, nfp=NFP,
         ladder=np.array([[n,M,N] for n,M,N in LADDER],dtype=object))
json.dump({"levels":levels,"nesting_proof":proof,"numeric_nesting_resid":num_resid,
           "feasibility":feas,"robustness_families":robustness},
          open(f"{HERE}/nested_basis_manifest.json","w"),indent=2)

print("NESTED STELLARATOR-SYMMETRIC FOURIER LADDER (single family, proven nested):")
for name,M,N in LADDER:
    d=levels[name]
    print(f"  {name}: M={M} N={N}  DOF={d['n_dof']:3d}  lam_pol_min={d['lambda_pol_min_cm']:.0f}cm "
          f"({d['finest_vs_coil_pol']})  lam_tor_min={d['lambda_tor_min_cm']:.0f}cm  sub_coil_pol={d['sub_coil_pol']}")
print("\nNESTING PROOF:"); [print("  "+p) for p in proof]
print(f"  numeric lstsq residual (F1 mode in F2 basis): {num_resid:.1e}")
print(f"\nGRID NYQUIST CAP: m<= {M_NYQ} pol, n<= {N_NYQ} tor. "
      f"Sub-coil poloidal reachable (m7-9); 8-coil toroidal NOT (needs n=8 > {N_NYQ}).")
print("wrote nested_basis_modes.npz + nested_basis_manifest.json")
