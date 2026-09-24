# Continuous Geological Particle Size Distribution Framework

Type: prototype
Status: resolved
Blocked by: 

## Question

How should geological and synthetic particle size distributions (Log-normal, Weibull, Rosin-Rammler, truncated Gaussian, and custom empirical sieving curves) be formulated, parameterized, and sampled with analytical number-weighted vs volume-weighted conversions?

Specifically:
1. In geological sieve analysis, distributions are measured by weight/volume fraction, whereas particle generation algorithms sample discrete particles by number frequency. How should the API bridge number-frequency vs volume-frequency PSDs analytically?
2. What geological parameters should be exposed as first-class inputs (e.g. median grain diameter $d_{50}$, Folk-Ward sorting coefficient $\sigma_\phi$, skewness, uniformity coefficient $C_u = d_{60}/d_{10}$)?
3. How should truncation bounds $[r_{\min}, r_{\max}]$ be handled to prevent unphysical sub-micron or oversized single grains that violate REV constraints?

## Answer

### 1. Analytical Volume-to-Number Frequency Bridge
- **The Physics**: Laboratory sieve analyses measure volume/mass fraction $f_V(r)$. Sampling grains by number requires the number frequency distribution:
  $$f_N(r) = \frac{f_V(r) / r^3}{\int_{r_{min}}^{r_{max}} \frac{f_V(u)}{u^3} \mathrm{d}u}$$
- **Log-Normal Invariance (Hatch-Choate Equation)**:
  When $f_V(r)$ is log-normal with volume median diameter $d_{50, V}$ and geometric standard deviation $\sigma_g$, $f_N(r)$ is also strictly log-normal with identical variance $\sigma_g$, but with a shifted number median diameter:
  $$\ln d_{50, N} = \ln d_{50, V} - 3 \ln^2 \sigma_g$$
- Sampling from $f_N(r)$ guarantees that the generated synthetic rock's cumulative solid volume matches laboratory sieve curves exactly.

### 2. Supported Geological Distributions
- `LogNormal(d50, sigma_phi)`: Supports Folk-Ward geological sorting scale ($\phi = -\log_2(d)$):
  - Well sorted ($\sigma_\phi < 0.50$)
  - Moderately sorted ($\sigma_\phi \in [0.71, 1.00]$)
  - Poorly sorted ($\sigma_\phi > 1.00$)
- `Weibull(scale, shape)` / `RosinRammler`: Common for crushed gravel, proppants, and mechanical sediments.
- `TruncatedGaussian(mean, std, bounds)`.
- `EmpiricalSieve(diameters, cumulative_percent_passing)`: Cubic spline interpolation of experimental laboratory sieve analysis.

### 3. Truncation and REV Bounds Enforcement
- Bounds $[r_{min}, r_{max}]$ are strictly enforced via truncated inverse-CDF transforms.
- Built-in validation check:
  - Error if $2 r_{max} > \min(L_x, L_y, L_z) / 3$ (ensures domain contains a valid Representative Elementary Volume and satisfies the 3-cell periodic boundary guard condition).
  - Warning if $r_{min} < 0.5 \times \text{voxel\_size}$ (sub-voxel resolution loss during rasterization).
