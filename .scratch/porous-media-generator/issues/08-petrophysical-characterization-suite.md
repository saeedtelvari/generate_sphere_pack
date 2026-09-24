# Integrated Petrophysical Characterization Suite

Type: grilling
Status: resolved
Blocked by: 01

## Question

What native petrophysical diagnostics should be built into the core library (connected vs isolated porosity, specific surface area via marching cubes/Crofton, 3D two-point correlation $S_2(r)$, chord length distribution, Kozeny-Carman permeability) versus delegating complex pore network extraction to PoreSpy/OpenPNM?

Specifically:
1. What morphological properties must run natively without any heavy external dependencies (total/effective porosity via `scipy.ndimage.label`, specific surface area $S_v$ via marching cubes, $S_2(r)$ via FFT autocorrelation)?
2. For pore network extraction (coordination number, pore body diameter, pore throat diameter & length), should we wrap PoreSpy's SNOW2 algorithm as an optional plugin, or provide a lightweight native watershed partitioner?
3. For transport properties: how should permeability estimators be staged (analytical Kozeny-Carman $k \propto \phi^3 / S_v^2$, Berg equation, vs numerical flow solvers)?

## Answer

### 1. Native Diagnostics (Pure NumPy/SciPy/scikit-image)
The following properties run natively without heavy external dependencies:
- **Porosity Partitioning**:
  - Total porosity: $\phi_t = 1 - N_{solid} / N_{total}$.
  - Effective / Percolating porosity $\phi_{eff}$: 3D connected-component labeling (`scipy.ndimage.label`) identifying void clusters spanning between inlet and outlet boundaries.
  - Isolated porosity: $\phi_{iso} = \phi_t - \phi_{eff}$.
- **Specific Surface Area ($S_v$)**: Evaluates interfacial surface area $A_{interface}$ via marching cubes triangulation (`skimage.measure.marching_cubes`) or the Crofton formula, divided by domain volume: $S_v = A_{interface} / V_{box}$.
- **Two-Point Spatial Correlation ($S_2(r)$)**: Computed via 3D FFT autocorrelation of the binary pore indicator field with radial azimuthal averaging, providing instant verification of Representative Elementary Volume (REV) and characteristic correlation length.
- **Directional Chord Lengths**: Fast ray-sampling across $x, y, z$ to extract pore-body and solid-grain chord distributions.

### 2. Analytical Transport Property Estimators
- **Kozeny-Carman Permeability**:
  $$k_{\text{KC}} = \frac{\phi_{eff}^3}{c \cdot S_v^2} \approx \frac{\phi_{eff}^3}{5 \cdot S_v^2}$$
- **Formation Factor ($F$) & Archie Exponent ($m$)**:
  Estimated via tortuosity $\tau$ approximation ($F \approx \tau / \phi_{eff} \approx \phi_{eff}^{-m}$, $m \approx 1.8 - 2.1$).

### 3. Optional Pore Network Modeling Bridge (PoreSpy / OpenPNM)
- `VoxelGrid.extract_pore_network()` cleanly wraps `porespy.networks.snow2` when `porespy` is installed.
- Returns an OpenPNM-compatible network dictionary containing pore body diameters, throat diameters, throat lengths, and coordination number histograms without hard-coding external dependencies into the core generator.
