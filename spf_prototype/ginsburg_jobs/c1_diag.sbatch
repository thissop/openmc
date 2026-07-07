#!/bin/bash
#SBATCH -A astro
#SBATCH -J spf_c1diag
#SBATCH -p short
#SBATCH -t 0-00:20
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=6G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
module load anaconda/3-2023.09 2>/dev/null
cd /ginsburg/astro/users/tjk2147/spf_work/spf_pp/python
"$HOME/desc_venv/bin/python" - <<'PY'
import numpy as np, quasr_equilibrium_field as qef
for ID in (803097, 886079):
    nfp, Rm, Zm, meta = qef._quasr_boundary_modes(ID)
    # axisymmetric (n=0) part: R(theta)=sum rc cos(m th), Z=sum zs sin(m th)
    r0 = {m: rc for (m, n), (rc, rs) in Rm.items() if n == 0 and abs(rc) > 1e-9}
    z0 = {m: zs for (m, n), (zc, zs) in Zm.items() if n == 0 and abs(zs) > 1e-9}
    th = np.linspace(0, 2*np.pi, 400, endpoint=False)
    R = sum(rc*np.cos(m*th) for m, rc in r0.items())
    Z = sum(zs*np.sin(m*th) for m, zs in z0.items())
    area = 0.5*abs(np.sum(R*np.gradient(Z, th) - Z*np.gradient(R, th)))
    # signed poloidal jacobian sign (does theta wind the 'right' way for DESC?)
    Rt, Zt = np.gradient(R, th), np.gradient(Z, th)
    Rr = R - R.mean()  # crude radial dir
    jac = Rr*Zt - (Z-Z.mean())*Rt   # proportional to the local dtheta orientation
    print(f"ID {ID} nfp {nfp}: n=0 R modes {sorted(r0)} Z modes {sorted(z0)}")
    print(f"  axisym cross-section: R[{R.min():.3f},{R.max():.3f}] Z[{Z.min():.3f},{Z.max():.3f}] "
          f"area {area:.4f}  (m=1 R={r0.get(1,0):.3f} Z={z0.get(1,0):.3f})")
    print(f"  theta-orientation sign: mean(jac) {np.mean(jac):+.3e}  (consistent sign? "
          f"{np.all(jac>0) or np.all(jac<0)})")
    # try building the DESC surface + its own jacobian check
    try:
        from desc.geometry import FourierRZToroidalSurface
        surf = qef._desc_surface(nfp, Rm, Zm)
        from desc.grid import LinearGrid
        g = LinearGrid(rho=1.0, theta=64, zeta=64, NFP=nfp)
        d = surf.compute(["sqrt(g)"], grid=g) if hasattr(surf, "compute") else None
        print(f"  DESC surface built OK; (surface has compute: {d is not None})")
    except Exception as e:
        print(f"  DESC surface issue: {type(e).__name__}: {e}")
PY
echo DONE_C1DIAG
