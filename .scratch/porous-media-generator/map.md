## Destination

A high-performance, modular Python library for generating, diagenetically transforming, and petrophysically characterizing 3D multi-purpose porous media (granular sphere packs, gravitational sedimentation beds, continuum Gaussian Random Fields, and QSGS) tailored for Digital Rock Physics (DRP) workflows, with clean exports to VTK, TIFF stacks, and reservoir/pore-network simulation tools.

## Notes

- **Domain**: Digital Rock Physics (DRP), Porous Media Generation, Geomechanics/Diagenesis, Petrophysics.
- **Skills to consult**: `wayfinder`, `domain-modeling`, `grilling`, `codebase-design`.
- **Standing preferences**:
  - $O(N)$ spatial indexing for all collision checks; avoid $O(N^2)$ brute-force comparisons.
  - Native Periodic Boundary Conditions (PBC) in $(x,y)$ and $(x,y,z)$ to eliminate artificial boundary skin layers.
  - Decoupled Lagrangian (`GrainPack`) and Eulerian (`VoxelGrid`) architectures.
  - Self-contained, dependency-light rasterization without forcing external C++ dependencies for basic usage, but with optional deep integration with PoreSpy and OpenPNM.
  - Strict preservation of physical conservation laws (analytical vs voxelized porosity, grain size distribution fidelity).

## Decisions so far

<!-- the index: one line per closed ticket, enough to judge relevance, then zoom the link for the detail the ticket holds -->

- [01: Core Lagrangian & Eulerian Data Structures](issues/01-core-representations.md): Decoupled GrainPack and VoxelGrid representations with native local bounding-box periodic rasterizer and PoreSpy/OpenPNM interop.
- [02: Spatial Indexing & Collision Detection Engine](issues/02-spatial-indexing-engine.md): Uniform spatial hash grid selected over cKDTree for O(1) dynamic insertions, 1.5s per 10^5 spheres, and native integer periodic wrapping.
- [03: Continuous Geological PSD Framework](issues/03-continuous-psd-framework.md): Geological Log-normal, Weibull, and sieve curves with analytical Hatch-Choate volume-to-number frequency conversion.
- [04: Ballistic Gravitational Sedimentation (Drop-and-Roll)](issues/04-drop-and-roll-sedimentation.md): Exact 4-stage analytical kinematics without DEM time-stepping, generating realistic loose sand porosity (0.39-0.42) and isostatic coordination (Z~6) in ~3s for 10^5 grains.
- [05: Dense Packing via Force-Biased Relaxation](issues/05-dense-packing-relaxation.md): Vectorized Hookean spring potential with FIRE optimization and adaptive staged inflation, reaching strict hard-sphere dense random packing (phi >= 0.64) in <45s.
- [06: Spectral Gaussian Random Fields & QSGS](issues/06-spectral-grf-generator.md): 3D FFT spectral filtering (Gaussian and Von Kármán kernels) with exact analytical inverse-CDF thresholding and anisotropic QSGS directional growth.
- [07: Post-Depositional Diagenesis Engine](issues/07-diagenetic-transformation-engine.md): Two-stage sequential diagenesis combining Lagrangian vertical uniaxial compaction (strain tensor) with Eulerian Euclidean distance transform cementation.
- [08: Integrated Petrophysical Characterization Suite](issues/08-petrophysical-characterization-suite.md): Native pure NumPy/SciPy diagnostics for connected porosity, specific surface area (Sv), S2(r) autocorrelation, and Kozeny-Carman permeability, with optional PoreSpy/OpenPNM network bridges.
- [09: Export Formats & Simulation Interfaces](issues/09-export-and-simulation-interfaces.md): Exporters for 3D TIFF stacks, VTK (.vti/.vtp), CSV/Parquet, and direct MRST grid/rock structs (.mat).

## Not yet specified

<!-- see "Fog of war": in-scope fog you can't ticket yet; graduates as the frontier advances -->

- **Multi-mineral diagenesis**: Explicitly tagging individual grains as quartz, feldspar, or carbonate, and modeling mineral-specific cementation rates or preferential clay coatings (pore-lining illite/chlorite).
- **Direct Stokes & Navier-Stokes flow solver**: Evaluating whether to include a built-in lattice Boltzmann method (LBM) or finite-difference Stokes flow solver, or rely strictly on PoreSpy/OpenPNM/MRST integration.
- **Two-phase displacement (drainage & imbibition)**: Developing a native invasion percolation with trapping engine to generate primary drainage capillary pressure curves $P_c(S_w)$ and relative permeability curves.
- **Generative AI / 3D Diffusion Models**: Training or wrapping pre-trained 3D score-based diffusion models for conditioned synthetic micro-CT rock generation (conditioned on porosity, permeability, and mineralogy).
- **Anisotropic geostatistical conditioning**: Reconstructing 3D porous volumes conditioned directly on 2D micro-CT thin-section slices via simulated annealing (Yeong-Torquato).

## Out of scope

<!-- see "Out of scope": work ruled beyond the destination; closed, never graduates -->

- Full macroscopic reservoir simulation (handled externally by MRST, CMG, or tNavigator; our responsibility ends at generating and characterizing the rock/grid).
- Direct physical mechanical testing of rock cores under triaxial compression apparatus (finite-element geomechanical solvers are left to specialized tools).
