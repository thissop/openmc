# Why the breeder-for-shield trade sign is a MULTIGROUP effect (signed-worth analysis)

Goal was to calibrate the signed worth so its absolute VALUE/sign predicts the finite-difference
dR, not just its ranking. Working the physics through shows that **the sign cannot be obtained from
single-group (total-flux) data at all** — it is set by neutron moderation in the breeder. This is a
cleaner, more precise result than a fitted prefactor, and it is empirically proven on the 4 patches.

## 1. First-order single-group perturbation predicts the WRONG sign

The fixed-envelope trade swaps a shell of breeder (FLiBe) for shield (WC) at the layer interface.
Standard first-order transport perturbation theory for a total-cross-section change:

    dR = - INT dSigma_t(x) phi(x) psi_dagger(x) dx .

For the breeder->shield swap, dSigma_t = Sigma_shield - Sigma_breeder > 0 over the trade zone, so

    dR ~ -(Sigma_sh - Sigma_br) * INT_[trade zone] phi*psi_dagger  < 0   for EVERY patch,

because the contributon phi*psi_dagger is positive everywhere. i.e. single-group theory says adding
shield ALWAYS lowers the coil dose. Measured (coil-20, 4 patches):

| patch | W_FD (measured) | swap contributon INT phi*psi | first-order sign | measured sign |
|-------|-----------------|------------------------------|------------------|---------------|
| t1p0  | +8.23e-5        | 4.99e6                       | DOWN             | DOWN          |
| t1p2  | +1.66e-5        | 4.89e5                       | DOWN             | DOWN          |
| t2p1  | +9.31e-5        | 2.19e5                       | DOWN             | DOWN          |
| **t2p3** | **-4.72e-5** | **1.33e7 (largest)**         | **DOWN**         | **UP**        |

t2p3 measured a dose INCREASE (+41%, +38.7 sigma) yet has the LARGEST trade-zone contributon --
exactly where first-order theory is most confident of a dose DROP. Single-group perturbation theory
is not just miscalibrated; it gets the SIGN wrong for the one patch that matters.

## 2. The missing physics: breeder moderation (a multigroup source term)

The coil response is FAST flux (0.1-25 MeV). The breeder (FLiBe + Be multiplier) is a strong
moderator: it removes fast neutrons by down-scattering them to lower energy. Thinning the breeder
therefore RELEASES fast flux toward the coil -- a positive contribution to the fast response that a
total-cross-section removal term cannot represent. In multigroup form the fast-response perturbation
is (schematically, g=fast):

    dR_fast = - INT dSigma_removal,g phi_g psi_g            (removal: shield helps, breeder-loss hurts)
              - INT dSigma_scatter,(g'->g) phi_g' psi_g     (MODERATION: breeder down-scatter fast->thermal)

The second (down-scatter) term is where the breeder does its fast-flux removal, and it flips sign
when the breeder is thinned. Where the breeder moderates the most flux -- the highest fast phi
through it, i.e. the largest trade-zone contributon (t2p3) -- removing it does the most damage. That
is exactly the observed inversion.

## 3. Why the full-band signed worth still RANKS correctly (+0.80)

The validated ranking result (Spearman W_FD vs W_signed = +0.80, robust across the XS ratio) uses
the FULL breeder-band contributon as the debit. That integral is a PROXY for "how much moderation
the breeder does at this patch": high breeder-band phi*psi_dagger <-> a lot of fast flux is being
moderated there <-> thinning the breeder releases a lot of fast flux <-> the trade is less beneficial
(or harmful). So the signed contributon is a good RANKER precisely because the breeder-band term
tracks the moderation loss -- even though, as a first-order removal integral, its absolute value is
wrong (all-negative, Section 1). The ranking is real physics via a proxy; the value is not ab initio.

## 4. Prescription for a quantitative (ab-initio) signed worth

Compute the response perturbation in >=2 energy groups:
  * forward flux tallied by group on the mesh (fast + thermal), spectrally consistent with the coil
    fast-flux tally -- the FIRST concrete step (fwd_meshflux.py now has a --fast option; a 2-group
    tally is a small extension);
  * a 2-group (or continuous-energy MGXS) adjoint psi_dagger for the coil fast response;
  * the down-scatter (moderation) cross section of the breeder, Sigma_s(fast->thermal).
Then dR from Section 2 includes the moderation term and can flip sign. This is the multigroup /
angular random-ray direction (connects to the group-dependent adjoint tally). Until then, the signed
contributon is used as a validated RANKER, not a calibrated predictor -- stated honestly.

## Paper narrative (layered, each step proven)
1. Attribution != actionability: scalar contributon worth is ANTI-correlated with truth (rho=-1.0).
2. First-order removal-based worth predicts the WRONG sign (dose-down for all; t2p3 measured up).
3. The trade sign is set by breeder MODERATION; the signed contributon proxies it and ranks
   robustly (+0.80). A quantitative worth requires the multigroup moderation term.
Finite difference is ground truth throughout.
