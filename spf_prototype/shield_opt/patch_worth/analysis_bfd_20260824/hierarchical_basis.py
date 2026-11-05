#!/usr/bin/env python
"""RUNG 5 (step d): hierarchical, field-period-compatible spatial basis for the shield-thickness design field
t(theta_pol, phi_tor) on the ParaStell (13 toroidal x 19 poloidal) grid, nfp=4. Records physical wavelength and
DOF at every level so the spatial-resolution CONVERGENCE study can add levels and watch the controllable subspace.

Geometry constants (exact, from repo notes): R0 = 1367 cm, a_minor = 217 cm, nfp = 4.
Stellarator symmetry: t(theta,phi) = t(-theta,-phi) -> retain COSINE (even) modes cos(m*theta - n*nfp*phi).
No proxy: this is the design-variable basis (dimensionless mode shapes); the CLEARANCE bounds that turn it into a
feasible domain come from the exact ParaStell/DAGMC surface (cluster job), applied later.

Physical wavelengths:
  poloidal mode m -> lambda_pol = 2*pi*a / m               (a=217 cm)
  toroidal mode n (n oscillations per field period)        -> arc/period = 2*pi*R0/nfp; lambda_tor = (2*pi*R0/nfp)/n
Level A: coarse sectors (uniform, inboard/outboard, up/down, the existing 4x4=16 patch indicator basis).
Level B: low-order stellarator-symmetric Fourier (m<=3, n<=2).
Level C: higher Fourier (m<=6, n<=4) approximating overlapping smooth B-splines (~30-60 effective DOF).
"""
import os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R0, A_MINOR, NFP = 1367.0, 217.0, 4
NT, NPOL = 13, 19
tor = np.linspace(0, 90, NT); pol = np.linspace(0, 360, NPOL)
TH = np.radians(pol)[None, :] * np.ones((NT, 1))       # poloidal angle rad, shape (NT,NPOL)
PH = np.radians(tor)[:, None] * np.ones((1, NPOL))     # toroidal angle rad (0..pi/2 over one period)

def lam_pol(m): return (2*np.pi*A_MINOR/m) if m else np.inf
def lam_tor(n): return ((2*np.pi*R0/NFP)/n) if n else np.inf

levels = {}

# ---- Level A: coarse sectors (interpretable indicator/large-scale modes) ----
A = {}
A["uniform"] = np.ones((NT, NPOL))
A["inboard"] = 0.5*(1 + np.cos(TH - np.pi))            # peak poloidal 180 (inboard)
A["outboard"] = 0.5*(1 + np.cos(TH))                   # peak poloidal 0/360 (outboard)
A["updown"] = 0.5*(1 + np.cos(2*TH))                   # m=2 up/down
# the existing 4x4=16 patch indicator sectors (coarse B_FD basis)
for it in range(4):
    for jp in range(4):
        ind = np.zeros((NT, NPOL))
        te = np.linspace(0, 90, 5); pe = np.linspace(0, 360, 5)
        ti = (np.clip(np.digitize(tor, te)-1, 0, 3) == it)
        pj = (np.clip(np.digitize(pol, pe)-1, 0, 3) == jp)
        ind[np.ix_(np.where(ti)[0], np.where(pj)[0])] = 1.0
        A[f"sector_t{it}p{jp}"] = ind
levels["A_coarse"] = {"modes": list(A), "n_dof": len(A),
                      "min_lambda_pol_cm": round(2*np.pi*A_MINOR/2, 0),  # up/down m=2 -> ~682
                      "min_lambda_tor_cm": round((2*np.pi*R0/NFP)/1, 0), # sector ~one lobe per period
                      "note": "interpretable large-scale + the 16 coarse patch sectors (the coarse B_FD basis)"}

# ---- Level B: low-order stellarator-symmetric Fourier cos(m*th - n*nfp*ph) ----
def fourier_modes(m_max, n_max):
    modes = {}
    for m in range(0, m_max+1):
        for n in range(0, n_max+1):
            if m == 0 and n == 0:
                modes["c0_0"] = np.ones((NT, NPOL)); continue
            modes[f"c{m}_{n}"] = np.cos(m*TH - n*NFP*PH)
            # stellarator symmetry keeps cosines; sines break it -> excluded (documented)
    return modes
B = fourier_modes(3, 2)
levels["B_lowfourier"] = {"modes": list(B), "n_dof": len(B), "m_max": 3, "n_max": 2,
    "lambda_pol_cm_by_m": {m: round(lam_pol(m), 0) for m in range(1, 4)},
    "lambda_tor_cm_by_n": {n: round(lam_tor(n), 0) for n in range(1, 3)},
    "note": "stellarator-symmetric cosines only; finest here lambda_pol~455cm (m=3), lambda_tor~1074cm (n=2)"}

# ---- Level C: higher-order (m<=6, n<=4) -> ~30-60 effective smooth DOF ----
C = fourier_modes(6, 4)
levels["C_highfourier"] = {"modes": list(C), "n_dof": len(C), "m_max": 6, "n_max": 4,
    "lambda_pol_cm_by_m": {m: round(lam_pol(m), 0) for m in range(1, 7)},
    "lambda_tor_cm_by_n": {n: round(lam_tor(n), 0) for n in range(1, 5)},
    "note": "finest lambda_pol~227cm (m=6, ~1 coil footprint), lambda_tor~537cm (n=4); ~= overlapping smooth B-splines"}

# save the actual mode arrays (for building design families + fitting B_FD onto each level)
np.savez(f"{HERE}/hierarchical_basis_modes.npz",
         **{f"A__{k}": v for k, v in A.items()},
         **{f"B__{k}": v for k, v in B.items()},
         **{f"C__{k}": v for k, v in C.items()},
         tor=tor, pol=pol, R0=R0, a=A_MINOR, nfp=NFP)
json.dump({"grid": [NT, NPOL], "R0_cm": R0, "a_cm": A_MINOR, "nfp": NFP, "levels": levels},
          open(f"{HERE}/hierarchical_basis_manifest.json", "w"), indent=2)

print("HIERARCHICAL BASIS (field-period-compatible, stellarator-symmetric):")
for L, d in levels.items():
    print(f"  {L}: {d['n_dof']} DOF")
print(f"\nphysical resolution ladder (poloidal lambda / coil footprint ~{2*np.pi*A_MINOR/6:.0f}-500 cm):")
print(f"  A coarse sectors: ~device-scale")
print(f"  B m<=3,n<=2: lambda_pol {lam_pol(3):.0f}-{lam_pol(1):.0f} cm, lambda_tor {lam_tor(2):.0f}-{lam_tor(1):.0f} cm")
print(f"  C m<=6,n<=4: lambda_pol {lam_pol(6):.0f}-{lam_pol(1):.0f} cm (finest ~1 coil), lambda_tor {lam_tor(4):.0f}+ cm")
print("wrote hierarchical_basis_modes.npz + manifest.json")
