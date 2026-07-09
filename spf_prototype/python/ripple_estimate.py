"""TF-ripple NWL back-of-envelope: toroidal ripple vs poloidal in/out, and the
exp(-N*a/R0) line-of-sight low-pass suppression. See docs/notes/RIPPLE_AND_PREDICTOR.md.

Filamentary ring source on the magnetic axis circle (R0, Z=0), isotropic emission,
toroidal emissivity S(phi)=1+eps_s*cos(N*phi); wall = torus of minor radius a about
that circle. First-flight, free-streaming. Any O(delta) N-periodic source channel
(strength, radial shift, SPF B-hat tilt) enters the same way.
"""
import numpy as np

R0, a, N, delta = 620.0, 200.0, 18, 0.005
eps_s = delta

Ns = 4000
phi_s = np.linspace(0, 2 * np.pi, Ns, endpoint=False)
dphi = 2 * np.pi / Ns
S = 1.0 + eps_s * np.cos(N * phi_s)
xs, ys, zs = R0 * np.cos(phi_s), R0 * np.sin(phi_s), np.zeros(Ns)


def nwl(phi_w, theta):
    Rw = R0 + a * np.cos(theta)
    xw, yw, zw = Rw * np.cos(phi_w), Rw * np.sin(phi_w), a * np.sin(theta)
    cx, cy = R0 * np.cos(phi_w), R0 * np.sin(phi_w)
    nx, ny, nz = xw - cx, yw - cy, zw
    nn = np.sqrt(nx * nx + ny * ny + nz * nz)
    nx, ny, nz = nx / nn, ny / nn, nz / nn
    dx, dy, dz = xw - xs, yw - ys, zw - zs
    d2 = dx * dx + dy * dy + dz * dz
    cos_inc = np.clip((nx * dx + ny * dy + nz * dz) / np.sqrt(d2), 0, None)
    return np.sum(S * cos_inc / d2 * dphi) / (4 * np.pi)


def modamp(q):
    return (q.max() - q.min()) / (2 * q.mean())


if __name__ == "__main__":
    q_out, q_in = nwl(0.0, 0.0), nwl(0.0, np.pi)
    qpol = np.array([nwl(0.0, t) for t in np.linspace(0, 2 * np.pi, 181)])
    phiw = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    tor_out = modamp(np.array([nwl(p, 0.0) for p in phiw]))
    tor_in = modamp(np.array([nwl(p, np.pi) for p in phiw]))
    print(f"R0={R0} a={a} a/R0={a/R0:.3f} N={N} delta={delta:.2%} eps_s={eps_s:.2%}")
    print(f"POLOIDAL  q_in/q_out={q_in/q_out:.3f}  swing={(qpol.max()-qpol.min())/qpol.mean():.1%}")
    print(f"TOROIDAL  outboard={tor_out:.4%}  inboard={tor_in:.4%}  suppression={eps_s/tor_out:.0f}x")
    print(f"low-pass  exp(-N*a/R0)=exp(-{N*a/R0:.2f})={np.exp(-N*a/R0):.4f}")
    print(f"RATIO poloidal/toroidal = {((qpol.max()-qpol.min())/qpol.mean())/tor_out:.0f}x")
