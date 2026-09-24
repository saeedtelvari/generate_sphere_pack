# Continuum Porous Media: Spectral Gaussian Random Fields & QSGS

Type: prototype
Status: resolved
Blocked by: 01

## Question

How should the level-cut Gaussian Random Field (GRF) and Quartet Structure Generation Set (QSGS) generators be parameterized and implemented using 3D FFT and directional growth for synthetic carbonates, micrites, and anisotropic media?

Specifically:
1. For GRF: what spectral filter kernels (Gaussian: $e^{-(k/k_c)^2}$, exponential, or Von Kármán) produce microstructures that best resemble carbonate vugs, chalk, and spinodal-like porous networks?
2. How is exact analytical target porosity enforced through inverse error function ($\mathrm{erf}^{-1}$) thresholding on the standardized 3D field?
3. For QSGS: how should directional growth probabilities ($P_x, P_y, P_z$) and seed density $c_d$ be mapped to macroscopic anisotropy ratios (e.g. horizontal-to-vertical permeability ratio $k_h/k_v$)?

## Answer

### 1. 3D Spectral Filter Kernels for GRF
The GRF is synthesized in Fourier space via 3D FFT in $O(V \log V)$ time:
$$Z(\mathbf{x}) = \text{Re}\left(\mathcal{F}^{-1}\left\{ \sqrt{S(\mathbf{k})} \cdot \mathcal{F}\{\mathcal{N}(0, 1)\} \right\}\right)$$
where $k = \|\mathbf{k}\|$ is radial wavenumber and $k_c = 2\pi / l_c$ sets the spatial correlation length $l_c$:
- **Gaussian Kernel**: $S(k) \propto \exp\left(-(k / k_c)^2\right)$, generating smooth, bicontinuous spinodal-like pore structures.
- **Exponential Kernel**: $S(k) \propto \left(1 + (k / k_c)^2\right)^{-2}$, generating rougher, multi-scale surfaces characteristic of micritic carbonates.
- **Von Kármán / Matérn Kernel**: $S(k) \propto \left(1 + (k / k_c)^2\right)^{-(H + 1.5)}$, parameterized by Hurst roughness exponent $H \in (0, 1)$ to reproduce self-affine natural fracture and vug boundaries.

### 2. Exact Analytical Thresholding via Inverse Normal CDF
After standardizing the real-space field to zero mean and unit variance ($Z \leftarrow (Z - \mu) / \sigma \sim \mathcal{N}(0, 1)$), the threshold $\alpha$ for target porosity $\phi$ is computed analytically:
$$\alpha = \sqrt{2}\,\text{erf}^{-1}(2\phi - 1)$$
The solid matrix is partitioned as `matrix = (Z >= alpha)`, guaranteeing that voxelized porosity strictly matches the target $\phi$ without iteration.

### 3. QSGS Directional Growth for Anisotropic Media
- **Nucleation**: Randomly seeds solid core sites with probability $c_d \ll 1$.
- **Directional Growth**: Evaluates neighboring voxels in 26 directions with probabilities $(P_x, P_y, P_z)$. Setting $P_x, P_y > P_z$ naturally produces stratified, laminated porous media (mimicking shale or layered sandstone) with directional permeability anisotropy $k_h / k_v > 1$.
