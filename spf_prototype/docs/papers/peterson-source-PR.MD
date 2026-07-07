<!--
If you are a first-time contributor to OpenMC, please have a look at our
contributing guidelines:
https://github.com/openmc-dev/openmc/blob/develop/CONTRIBUTING.md
-->

# Description

This PR implements a new type of source based on a [common parametric description of a tokamak plasma](https://www.sciencedirect.com/science/article/abs/pii/S0920379612000853). A schematic of the tokamak shaping parameters can be seen in the image below where $R_0$ is the major radius, $a$ is the minor radius, $\delta$ is the triangularity, $\kappa$ is the elongation, $\Delta$ is the Shafranov shift, $\tilde{r}$ is the normalized minor radial coordinate $r/a$, and $\alpha$ is the poloidal angle measured counterclockwise from the outboard midplane.

<img width="1620" height="1662" alt="image" src="https://github.com/user-attachments/assets/cf8cec9f-731f-44e8-8061-715edb4030de" />

The parameterization defines the flux surfaces with the following equations:
```math
\large
\begin{aligned}
  R\left(r,\alpha\right) &= R_0 + r\cos\left(\alpha + \delta\sin\alpha\right) +\Delta\left(1-\left(\dfrac{r}{a}\right)^2\right)\\ 
  Z\left(r,\alpha\right) &= \kappa r\sin{\alpha},
\end{aligned}
```

The total source emission rate, $S_0$ in neutrons / second is given by the integral over the plasma volume of the source emission density $S(R,Z)$ as shown in the equation below.
```math
\large
    S_0 = \int_VS\left(R,Z\right)RdRdZd\phi
```
To sample the plasma conformally and continuously without rejection we want to transform coordinates from $R, Z$ to $r, \alpha$ because the plasma is defined to be on $r \in [0,a] \times \alpha \in [0, 2\pi]$ and the source emission density is a flux function and therefore only a function of $r$. This means the spatial probability distribution we need to sample from is (assuming axisymmetry so that $\phi$ can be sampled uniformly):

```math
\large
p(\tilde{r},\alpha) \propto S(\tilde{r})R(\tilde{r}, \alpha)\mathcal{J}(\tilde{r}, \alpha)
```
where $R(\tilde{r}, \alpha)$ and $\mathcal{J}(\tilde{r}, \alpha)$ are the major radius and jacobian of the transformation respectively. I won't work through the remainder of the math here, but the algorithm to sample the neutron source site proceeds as follows 

1. Sample $\tilde{r}$ from the marginal distribution $p(\tilde{r})$
2. Sample $\alpha$ from the conditional distribution $p(\alpha | \tilde{r})$
3. Sample $\phi$ uniformly over the toroidal extent of the plasma
4. Sample the neutron energy 
5. Sample the angle isotropically
6. Convert $(\tilde{r}, \alpha, \phi)$ coordinates to $(x, y, z)$ and return the source site.


Fixes # (issue)

# Checklist

- [x] I have performed a self-review of my own code
- [x] I have run [clang-format](https://docs.openmc.org/en/latest/devguide/styleguide.html#automatic-formatting) (version 18) on any C++ source files (if applicable)
- [x] I have followed the [style guidelines](https://docs.openmc.org/en/latest/devguide/styleguide.html#python) for Python source files (if applicable)
- [x] I have made corresponding changes to the documentation (if applicable)
- [x] I have added tests that prove my fix is effective or that my feature works (if applicable)
<!--
While tests will automatically be checked by CI, it is good practice to
ensure that they pass locally first. See instructions here:
https://docs.openmc.org/en/latest/devguide/tests.html
-->
