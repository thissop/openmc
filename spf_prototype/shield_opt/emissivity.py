"""Pluggable fusion-reactivity emissivity weight r(rho) for the culprit/NWL machinery.
Uniform (|sqrt g| only) isolates geometry+polarization; Bosch-Hale folds in a realistic
core-peaked DT reaction-rate profile so users can drop in their own plasma scenario.

r(rho) = n(rho)^2 * <sigma v>_DT(T(rho))   (times the voxel |sqrt g| applied by the caller)
Default profiles match Miralles-Dolz et al. (IEEE TPS 2026): n(s)=4.8e20(1-s^5) m^-3,
T(s)=11.5(1-s) keV, s=rho^2 (normalized toroidal flux). Bosch-Hale (1992) DT <sigma v>.
"""
import numpy as np

# Bosch-Hale 1992 DT reactivity (m^3/s), T in keV, valid ~0.2-100 keV
_BG=34.3827; _MRC2=1.124656e6
_C=[1.17302e-9,1.51361e-2,7.51886e-2,4.60643e-3,1.35000e-2,-1.06750e-4,1.36600e-5]
def sigma_v_DT(T_keV):
    T=np.maximum(np.asarray(T_keV,float),1e-3)
    theta=T/(1-(T*(_C[1]+T*(_C[3]+T*_C[5])))/(1+T*(_C[2]+T*(_C[4]+T*_C[6]))))
    xi=(_BG**2/(4*theta))**(1/3)
    sv=_C[0]*theta*np.sqrt(xi/(_MRC2*T**3))*np.exp(-3*xi)  # cm^3/s
    return sv*1e-6  # -> m^3/s

def profiles_miralles(s):
    n=4.8e20*(1-s**5); T=11.5*(1-s); return n,T

def emissivity_weight(rho, kind="bosch_hale", profile=profiles_miralles):
    """r(rho), normalized to max 1. kind='uniform' -> ones; 'bosch_hale' -> n^2<sigma v>(T)."""
    rho=np.asarray(rho,float)
    if kind=="uniform": return np.ones_like(rho)
    s=rho**2                      # normalized toroidal flux
    n,T=profile(s)
    r=n**2*sigma_v_DT(T)
    return r/ (r.max()+1e-300)

if __name__=="__main__":
    rho=np.linspace(0,1,12)
    for k in ("uniform","bosch_hale"):
        print(k, np.round(emissivity_weight(rho,k),3))
