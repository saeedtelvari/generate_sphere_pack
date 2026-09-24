# Core Lagrangian & Eulerian Data Structures

Type: prototype
Status: resolved
Blocked by: 

## Question

What are the core data structures for decoupled Lagrangian particles (`GrainPack`) and Eulerian 3D discrete fields (`VoxelGrid`), and how should continuous-to-voxel rasterization with periodic boundary conditions and sub-voxel antialiasing be designed?

Specifically:
1. How should `GrainPack` store coordinates, radii, bounding box dimensions, and periodic flags (e.g. structured NumPy record array or contiguous arrays)?
2. How should `VoxelGrid` handle memory-efficient storage (boolean array vs bitpack vs uint8) and spatial metadata (origin, voxel size)?
3. What is the fastest, self-contained rasterization algorithm to convert `GrainPack` into `VoxelGrid` without requiring heavy external C++ binaries, while supporting true periodic boundary wrap-around?

## Answer

### 1. `GrainPack` Lagrangian Data Structure
- **Storage**: Contiguous NumPy arrays for memory efficiency and vectorized math:
  - `coords`: `ndarray` of shape `(N, 3)`, `float64` (grain centroids $[x, y, z]$).
  - `radii`: `ndarray` of shape `(N,)`, `float64`.
  - `box_size`: `ndarray` of shape `(3,)`, `float64` $[L_x, L_y, L_z]$.
  - `periodic`: `tuple[bool, bool, bool]`, indicating whether boundaries wrap in $x, y, z$.
  - `attributes`: `dict[str, ndarray]`, optional grain properties (mineralogy IDs, contact counts).
- **Core Operations**: Fast affine transformations (compaction strain tensor), analytical non-overlapping porosity, filtering, and export.

### 2. `VoxelGrid` Eulerian Data Structure
- **Storage**: 3D NumPy boolean array (`bool_`, 1 byte/voxel) or `uint8` for multi-phase segmented media (0 = pore, 1 = quartz, 2 = clay, 3 = calcite cement).
  - Shape: `(Nx, Ny, Nz)` where $N_i = \mathrm{round}(L_i / \text{voxel\_size})$.
  - `voxel_size`: float isotropic physical voxel resolution $h$ (e.g. in $\mu\text{m}$ or relative units).
- **Integrity**: Exact numerical porosity computed as `1.0 - np.count_nonzero(matrix) / matrix.size`, eliminating overlapping analytical approximations.

### 3. Native Periodic Sphere Rasterizer
- **Algorithm**: Local Bounding-Box Stamping with Toroidal Wrap:
  1. For each sphere $(x_0, y_0, z_0, r)$, determine the voxel bounding box indices $[i_{min}, i_{max}], [j_{min}, j_{max}], [k_{min}, k_{max}]$ of size $\approx \lceil 2r / h \rceil^3$.
  2. Compute local Euclidean distance grid $(i \cdot h - x_0)^2 + (j \cdot h - y_0)^2 + (k \cdot h - z_0)^2 \le r^2$.
  3. If periodic, map voxel coordinates modulo grid dimensions ($i \pmod{N_x}$, etc.) to paint intersecting spheres seamlessly across boundary faces, edges, and corners.
- **Interoperability**: Provides `.to_porespy()` and `.to_openpnm()` helper methods to pass the boolean grid directly into existing DRP packages without file I/O overhead.
