# Porous Media & Digital Rock Physics Generator

A Python library for generating, diagenetically transforming, and petrophysically characterizing 3D multi-purpose porous media (granular assemblies, continuum random fields, and synthetic digital rocks).

## Language

### Geometry & Spatial Domain

**GrainPack**:
A discrete Lagrangian representation of solid particles containing spatial coordinates $(x,y,z)$, radii $r$, and optional particle attributes (mineralogy, orientation, shape).
_Avoid_: SphereList, ParticleArray, CoordTable

**VoxelGrid**:
An Eulerian 3D discrete grid representing the phase occupancy (pore space vs solid grains/cements) with defined isotropic voxel spacing.
_Avoid_: BinaryImage, Matrix, Image3D

**Periodic Boundary Conditions (PBC)**:
Toroidal domain wrap-around in $(x,y)$ or $(x,y,z)$ ensuring that particles or structures intersecting one boundary seamlessly re-enter the opposing boundary without boundary skin effects.
_Avoid_: Infinite domain, Wrap-around

**Representative Elementary Volume (REV)**:
The minimum physical volume of porous media over which macroscopic petrophysical properties (porosity, specific surface area, permeability) remain statistically invariant with volume scale.
_Avoid_: Sample volume, Unit cell

---

### Generation Paradigms

**Random Sequential Adsorption (RSA)**:
Stochastic insertion of particles one by one into the domain with hard-sphere non-overlapping rejection checks until target porosity or jamming limit is reached.
_Avoid_: Brute-force placement, Random shot

**Sedimentation ("Drop-and-Roll")**:
A gravitational settling algorithm where grains fall along trajectories and roll over contacted surfaces until finding stable mechanical equilibrium (typically 3 contact points).
_Avoid_: Gravity drop, Ballistic landing

**Force-Biased Relaxation (FBA)**:
A collective rearrangement method where overlapping or growing particles interact via repulsive contact potentials and relax dynamically to achieve dense random close packing.
_Avoid_: Molecular dynamics, Particle pusher

**Level-Cut Gaussian Random Field (GRF)**:
A continuum method where 3D white Gaussian noise is filtered in Fourier space using spatial autocorrelation kernels and thresholded to yield bicontinuous porous media.
_Avoid_: Noise thresholding, Perlin rock

---

### Diagenesis & Petrophysics

**Compaction**:
A post-depositional uniaxial mechanical strain deformation that compresses the vertical axis and forces grain interpenetration.
_Avoid_: Squeezing, Flattening

**Cementation**:
Morphological precipitation of secondary minerals (quartz overgrowth, calcite) expanding outward from grain boundaries to fill pore throats.
_Avoid_: Grain growing, Solid dilation

**Effective Porosity**:
The volume fraction of the domain occupied by interconnected, percolating pore space, excluding isolated non-flowing void clusters.
_Avoid_: Connected porosity, Open void fraction

**Two-Point Correlation Function ($S_2(r)$)**:
The probability that two random points separated by Euclidean distance $r$ both fall within the pore phase.
_Avoid_: Autocorrelation, Pair correlation
