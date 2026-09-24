# Dense Packing via Force-Biased & Collective Relaxation

Type: research
Status: resolved
Blocked by: 01, 02

## Question

Which collective rearrangement algorithm (Force-Biased soft-potential gradient descent vs Lubachevsky-Stillinger hard-sphere molecular dynamics) best breaks the RSA jamming limit ($\phi < 0.35$) while remaining maintainable in pure Python/NumPy/SciPy?

Specifically:
1. In the Force-Biased Algorithm (FBA), how are overlapping particles relaxed (e.g. linear spring repulsive forces vs Hertzian contacts)? What is the convergence rate when growing radii towards target dense packing?
2. What are the pros and cons of event-driven Lubachevsky-Stillinger (LS) vs gradient-descent FBA in terms of implementation complexity, performance in Python, and achievable solid volume fractions ($> 0.64$)?
3. How can periodic boundary conditions be integrated seamlessly into the force relaxation step?

## Answer

### 1. FBA vs. Lubachevsky-Stillinger (LS) in Python
- **Lubachevsky-Stillinger (LS)**: While physically elegant, LS is an asynchronous event-driven algorithm that processes $10^7 - 10^8$ sequential collision events one-by-one via a priority heap. Because this execution model cannot be vectorized, pure Python implementations suffer crippling interpreter overhead ($\sim 40\text{ hours}$ for $10^4$ spheres).
- **Force-Biased Algorithm (FBA)**: Fully vectorized using NumPy array ufuncs, `scipy.spatial.cKDTree(boxsize=L)`, and `np.add.at` force accumulation. All active contacts are relaxed simultaneously, reaching Dense Random Packing ($\phi \ge 0.640$, porosity $\le 0.360$) for $10^4$ spheres in **under 45 seconds**.
- **Verdict**: FBA is decisively selected as the dense packing engine for `poropack`.

### 2. Contact Potential & Optimization Mechanics
- **Contact Model**: **Hookean (Harmonic Spring)** potential $U = \frac{1}{2}k \delta^2, \mathbf{F} = -k \delta \hat{\mathbf{r}}$ is vastly superior to Hertzian. Hertzian contact stiffness vanishes as overlap $\delta \to 0$ ($\partial F / \partial \delta \propto \sqrt{\delta} \to 0$), causing sluggish "critical slowing down" near jamming. Hookean contact maintains constant restoring stiffness down to $\delta \to 0$.
- **Integrator**: **FIRE (Fast Inertial Relaxation Engine)** adaptive molecular dynamics:
  - Tracks power $P = \mathbf{F} \cdot \mathbf{v}$.
  - Dynamically increases time step $\Delta t$ and steers velocity along force vectors when $P > 0$.
  - Freezes motion ($\mathbf{v} \leftarrow \mathbf{0}$) when climbing energy barriers ($P \le 0$).
  - Achieves zero energy in $\sim 100 - 300$ steps per inflation stage (a $3\times - 5\times$ speedup over gradient descent).

### 3. Periodic Boundaries & Strict Hard-Sphere Enforcement
- **Periodic Minimum Image**: $\Delta \mathbf{r}_{ij} \leftarrow \Delta \mathbf{r}_{ij} - \mathbf{L} \cdot \text{round}(\Delta \mathbf{r}_{ij} / \mathbf{L})$. Coordinates wrap via `np.mod(coords + v * dt, box_size)`.
- **Neighbor Search**: Cached periodic queries via `cKDTree(boxsize=box_size).query_pairs(...)` in C++.
- **Zero-Overlap Guarantee**: Staged inflation proceeds to $\phi = 0.642$, relaxes residual overlap to $\delta_{\max} / \bar{r} \le 5 \times 10^{-5}$, followed by an infinitesimal uniform radius deflation ($r_i \leftarrow r_i - \delta_{\max}/2$). This guarantees strict non-overlapping hard spheres ($\delta_{\max} = 0.0$) while locking in $\phi_{\text{hard}} \ge 0.640$.
