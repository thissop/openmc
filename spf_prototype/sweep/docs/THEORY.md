# THEORY: field-direction coherence as the predictor of the SPF benefit

## The physical picture
Spin-polarized fuel (SPF) makes DT neutrons emit anisotropically about the local
magnetic field direction b_hat(x) = B/|B|. The emission kernel is
w(theta_B) ~ 1 + a2 P2(cos theta_B) + ..., with coefficients set by the
polarization state and theta_B the angle to b_hat. The steering benefit at a given
detector (a wall patch, a coil) is the coherent sum of these biased emissions over
all source points that illuminate that detector.

In a tokamak b_hat is nearly unidirectional across the burning volume (dominantly
toroidal, slowly varying), so the biased emissions from different source points add
coherently and the steering is strong. In a stellarator b_hat ROTATES across the
volume (helical excursion, pitch varying with position), so emissions from
different source points are biased in different directions and partially cancel at
any fixed detector. That cancellation erodes the SPF benefit. The amount of
cancellation is the physics under study.

## Why the independent variable is field-direction coherence, not QS error
The natural temptation is to order configs by "deviation from axisymmetry" or by
quasisymmetry (QS) error. Both are wrong for this problem.

Quasisymmetry is a symmetry of the field STRENGTH |B| in Boozer coordinates, not a
symmetry of the real-space field DIRECTION or of the plasma shape. Two configs with
identical QS error can have very different coherence of b_hat over the source
volume: QS constrains |B|, while the SPF kernel depends on the DIRECTION b_hat.
"Deviation from axisymmetry" is a shape statement, and shape does not map one to one
onto field-direction spread either. The quantity that actually controls the
coherent-versus-cancelling sum is the SOURCE-WEIGHTED DISPERSION OF b_hat itself, so
that is the independent variable we build the study around.

## The order parameter C (script-C in figures)
Define the source-weighted first moment of the field-direction distribution on the
unit sphere and its magnitude:

    C = | integral s(x) b_hat(x) dV | / integral s(x) dV ,   C in [0,1].

s(x) is the fusion source density the neutronics source actually uses (the
(1 - rho^2) birth profile of the compiled source; see coherence_metrics.source_weights).
C is exactly the mean RESULTANT LENGTH R-bar of directional statistics: the length
of the vector average of unit directions weighted by s. It is the magnitude of the
first trigonometric moment of the direction distribution on S^2 (Mardia and Jupp,
Directional Statistics). C = 1 means every source point shares one field direction
(perfectly coherent, tokamak limit, maximal SPF steering). C -> 0 means the
directions are dispersed over the sphere and the steering cancels.

### The frame is decisive: local cylindrical basis
C is only meaningful once the FRAME for b_hat is fixed, and this choice is not
cosmetic. We express b_hat in the LOCAL CYLINDRICAL basis (e_R, e_phi, e_Z) at each
point before taking the moment. Reasoning:

- A purely toroidal field (the tokamak limit) is b_hat = e_phi, which in the local
  cylindrical basis is (0, 1, 0) at EVERY point. Its resultant length is 1, so the
  tokamak limit calibrates to C = 1, matching the physical statement that tokamak
  steering is maximally coherent.
- In the fixed lab (x, y, z) frame the SAME toroidal field winds around the torus:
  e_phi points in a different lab direction at each toroidal angle, and the moment
  integrates to nearly zero. Lab-frame C is therefore near zero for a tokamak AND a
  stellarator, so it cannot discriminate between them. (Measured directly:
  device 59509 has cylindrical C = 0.997 but lab-frame C = 0.006.)

So the cylindrical basis is the frame in which "field direction is the same
everywhere" corresponds to the tokamak, which is the reference the SPF benefit is
measured against. A stellarator's helical and pitch variation makes b_hat depart
from (0, 1, 0) over the source, driving C below 1. This is the single most
important definitional choice in the study; if a field-line (Frenet) frame were
intended instead the metric would change, so the choice is stated here explicitly.

## Beyond the scalar: structure of the spread
C is a magnitude and loses the STRUCTURE of the dispersion. Two configs can share
|first moment| yet spread b_hat differently (a thin planar fan versus an isotropic
smear). We therefore also compute the source-weighted direction tensor

    T = integral s(x) b_hat b_hat^T dV / integral s(x) dV ,   trace(T) = 1,

and its eigenvalues (l1 >= l2 >= l3, summing to 1). A concentrated direction gives
(1, 0, 0); a planar fan gives roughly (a, 1-a, 0); an isotropic smear gives
(1/3, 1/3, 1/3). The anisotropy l1 - l2 and the full spectrum are the structure
predictors, used when the scalar C fails to collapse the classes (see below). We
also report the source-weighted angular standard deviation of b_hat about its mean
direction, and the directional-statistics circular standard deviation sqrt(-2 ln C).

## The hypothesis the code tests (not assumes)
The directional efficiency eta (fraction of the ideal SPF benefit surviving) is
hypothesized to collapse onto a single curve eta = eta_source(C) across symmetry
classes (QA, QH, QI), and the blanket rescales it by one scalar attenuation factor:

    eta(C, blanket) ~ eta_source(C) * A(tau_scatter).

Field geometry lives in eta_source(C); blanket physics lives in A. Both the
universality (one curve across classes) and the separability (factorization into
eta_source times A) are FALSIFIABLE outputs of the sweep, not assumptions. See
FACTORIZATION.md. If the classes do not collapse under scalar C, the direction
tensor structure (not C alone) is the real predictor, and that is itself a result.

## Naming
We call C the field-direction coherence, symbol script-C. The name is provisional:
it is earned only if C predicts eta with low scatter across the config family. If
the tensor structure is required instead, the scalar name is retired in favor of the
structural predictor.

## What is cheap here
C, the tensor, and the angular std need the field DIRECTION and the source weight
only. No neutronics. Given the coils they are seconds of compute per config
(Biot-Savart on a source-volume sample). The whole point is to predict the expensive
neutronics eta from this cheap geometric quantity, so that SPF benefit can be
screened across a device family without transport.
