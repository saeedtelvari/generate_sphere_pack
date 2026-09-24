# AGENTS.md: Developer & Agent Engineering Guide

This document defines the architecture, conventions, testing protocols, and operational workflows for autonomous AI coding agents and human developers contributing to or utilizing `poropack`.

---

## 1. Repository Purpose & Architecture

`poropack` is a Python library for generating, diagenetically transforming, and petrophysically characterizing 3D multi-purpose porous media (granular assemblies, continuum random fields, and synthetic digital rocks) for Digital Rock Physics (DRP) and reservoir simulation.

### Dual Representation Architecture
The codebase decouples discrete Lagrangian particles from continuous Eulerian voxel fields:

```
+---------------------------------------------------------------------------------+
|                                 LAGRANGIAN DOMAIN                               |
|                                                                                 |
|   GrainPack (poropack.core.representations)                                    |
|   - coords: (N, 3) float64 array of grain centroids                             |
|   - radii:  (N,) float64 array of sphere radii                                  |
|   - box_size: (3,) float64 bounding box dimensions [Lx, Ly, Lz]                 |
|   - periodic: (bool, bool, bool) toroidal boundary wrap flags                   |
|   - attributes: dict[str, ndarray] (contacts, mineralogy IDs, etc.)             |
|   - Methods: .analytical_porosity(), .wrap_coordinates(), .to_voxel(voxel_size) |
+---------------------------------------------------------------------------------+
                                      |
                                      | .rasterize(voxel_size) / .to_voxel()
                                      | via rasterize_spheres()
                                      v
+---------------------------------------------------------------------------------+
|                                  EULERIAN DOMAIN                                |
|                                                                                 |
|   VoxelGrid (poropack.core.representations)                                     |
|   - matrix / voxels: 3D bool or uint8 array (0 = pore space, >0 = solid matrix) |
|   - voxel_size: float physical isotropic resolution h                           |
|   - origin: (3,) float64 origin coordinates [x0, y0, z0]                        |
|   - periodic: (bool, bool, bool) boundary indicators                            |
|   - Methods: .porosity(), .solid_fraction(), .pore_mask, .to_porespy()          |
+---------------------------------------------------------------------------------+
```

---

## 2. Directory Layout & Module Index

```
generate_sphere_pack/
├── AGENTS.md                          # This engineering guide
├── CONTEXT.md                         # Authoritative domain vocabulary & nomenclature rules
├── README.md                          # User-facing overview, quickstart & installation
├── pyproject.toml                     # Modern package metadata & pytest configuration
├── requirements.txt                   # Production dependencies
├── genrandsp.py                       # Backward-compatible wrapper for legacy scripts
├── poropack/                          # Core production package
│   ├── __init__.py                    # Top-level API exports & versioning
│   ├── core/                          # Fundamental representations & spatial acceleration
│   │   ├── representations.py         # GrainPack, VoxelGrid, rasterize_spheres
│   │   ├── spatial.py                 # PeriodicSpatialGrid O(1) dynamic spatial hash grid
│   │   └── psd.py                     # Geological continuous & empirical PSDs
│   ├── generators/                    # Physics-grounded porous media generators
│   │   ├── rsa.py                     # Fast Random Sequential Adsorption
│   │   ├── sedimentation.py           # 3D Ballistic Sedimentation ("Drop-and-Roll")
│   │   ├── dense_packing.py           # Force-Biased Algorithm (FBA) + FIRE relaxation
│   │   └── continuum.py               # Spectral Gaussian Random Fields & QSGS
│   ├── diagenesis/                    # Post-depositional geological transformations
│   │   ├── compaction.py              # Uniaxial vertical mechanical strain deformation
│   │   └── cementation.py             # EDT syntaxial cementation & secondary dissolution
│   ├── petrophysics/                  # Digital rock characterization suite
│   │   ├── porosity.py                # Total, effective (percolating), isolated porosity
│   │   ├── morphology.py              # Specific surface area Sv, S2(r), chord distributions
│   │   ├── transport.py               # Kozeny-Carman permeability & formation factor
│   │   └── pnm.py                     # PoreSpy SNOW2 pore network extraction
│   └── io/                            # Simulation & visualization interfaces
│       ├── exporters.py               # TIFF, VTK (.vti/.vtp), MRST .mat, Parquet, CSV
│       └── visualization.py           # Matplotlib orthogonal slices & X-ray projections
├── notebooks/                         # Clean, topic-focused interactive tutorials
│   ├── 01_quickstart_and_psd.ipynb    # Introduction, continuous PSDs, and fast RSA
│   ├── 02_physical_generators.ipynb   # Sedimentation, FBA dense packing, and GRF/QSGS
│   ├── 03_diagenesis_and_petrophysics.ipynb # Compaction, cementation, dissolution, petrophysics
│   └── 04_reservoir_simulation_and_exports.ipynb # MRST integration & simulation exporters
├── tests/                             # Comprehensive automated pytest test suite (21 tests)
└── archive/                           # Preserved historical prototypes and legacy code
    ├── README.md                      # Archive explanation & migration notes
    └── legacy/                        # Legacy notebooks and procedural scripts
```

---

## 3. Strict Domain Vocabulary & Naming Rules

When modifying or generating code, agents must adhere to the conventions defined in [`CONTEXT.md`](CONTEXT.md):

| Domain Concept | Required Term / Class | Terms to Avoid |
|---|---|---|
| Lagrangian Particle Assembly | `GrainPack` | `SphereList`, `ParticleArray`, `CoordTable` |
| Eulerian Discrete Matrix | `VoxelGrid` (underlying `.voxels` / `.matrix`) | `BinaryImage`, `Matrix`, `Image3D` |
| Periodic Wrap-around | `periodic` / `Periodic Boundary Conditions` | `infinite domain`, `torus wrap` |
| Gravitational Settling | `Sedimentation` / `"Drop-and-Roll"` | `gravity drop`, `ballistic landing` |
| Dense Assembly Relaxation | `Force-Biased Algorithm` / `FBA` | `particle pusher`, `molecular dynamics` |
| Continuum Bicontinuous Media | `GaussianRandomField` (GRF) / `QSGS` | `Perlin rock`, `noise thresholding` |
| Post-depositional Strain | `Compaction` / `vertical_strain` | `squeezing`, `flattening` |
| Mineral Precipitation | `Cementation` | `grain growing`, `dilation` |
| Autocorrelation Function | `Two-Point Correlation Function` / $S_2(r)$ | `autocorrelation`, `pair correlation` |
| Flow-capable Void Fraction | `Effective Porosity` / `effective` | `connected porosity`, `open void fraction` |

---

## 4. Operational Environment Rules

1. **Python Environment**:
   - The user develops with conda / isolated Python environments. Ensure dependencies in `requirements.txt` and `pyproject.toml` are satisfied.
   - Run tests using:
     ```powershell
     pytest tests/ -v
     ```
2. **MATLAB / MRST Compatibility**:
   - When generating or testing MATLAB Reservoir Simulation Toolbox (MRST) scripts:
     - Always run `startup` from the working directory first to initialize MRST paths.
     - Include `mrstModule add ...` for all required modules (e.g. `mrstModule add incomp mrst-gui`).
     - `export_mrst_mat` generates compatible Cartesian geometry structs `G` (`cartDims`, `dimensions`, `cells.num`, `faces.num`) and rock structs (`rock.poro`, `rock.perm`).
3. **Non-Interactive Matplotlib Backend**:
   - In automated test scripts or batch execution, always configure non-interactive plotting before importing `pyplot`:
     ```python
     import matplotlib
     matplotlib.use('Agg')
     import matplotlib.pyplot as plt
     ```
   - Never leave naked `plt.show()` calls in batch automation scripts as it blocks execution awaiting GUI interaction.

---

## 5. Standard Code Recipes for Agents

### Recipe A: Generating a Granular Digital Rock
```python
import poropack as pp

# 1. Define sedimentological distribution with Hatch-Choate conversion
psd = pp.LogNormalPSD(d50=180.0, sigma_phi=0.30, bounds=(30.0, 350.0))

# 2. Generate RSA grain assembly with periodic boundaries
rsa = pp.RSAGenerator(
    box_size=(600.0, 600.0, 600.0),
    psd=psd,
    min_throat=2.0,
    periodic=(True, True, True),
    random_state=42,
)
pack = rsa.generate(target_porosity=0.62)

# 3. Rasterize to Eulerian VoxelGrid at 6 um resolution
voxel_grid = pack.rasterize(voxel_size=6.0)
```

### Recipe B: Applying Diagenetic Burial History
```python
# 1. 12% vertical mechanical compaction
compacted = pp.apply_compaction(voxel_grid, vertical_strain=0.12)

# 2. 8% syntaxial quartz overgrowth cementation
cemented = pp.apply_cementation(compacted, cement_fraction=0.08)

# 3. 3% intra-granular core dissolution
dissolved = pp.apply_dissolution(cemented, porosity_increase=0.03, mode="grain_cores")
```

### Recipe C: Petrophysical Characterization & Export
```python
# 1. Porosity breakdown
poro_res = pp.analyze_porosity(dissolved)
print(poro_res.summary())

# 2. Specific surface area and permeability
sv = pp.specific_surface_area(dissolved)
k_mD = pp.kozeny_carman(dissolved) / 9.869233e-16

# 3. Export to MRST and ParaView
pp.export_mrst_mat(dissolved, "sim_grid.mat", permeability_mD=k_mD)
pp.export_vtk_vti(dissolved, "matrix.vti")
pp.export_parquet(pack, "grains.parquet")
```

---

## 6. Verification & Quality Gates

Before committing any modifications:
1. **Pytest Verification**: All 21 tests in `tests/` must pass 100% green:
   ```powershell
   pytest tests/ -v
   ```
2. **Tutorial Notebook Verification**: All notebooks in `notebooks/` must execute without errors.
3. **Backward Compatibility**: Ensure `from genrandsp import generate_spherepack` functions identically for legacy code.
4. **Zero Overlap Guarantee**: Dense packing outputs must strictly verify that particle center distances exceed the sum of radii ($\Delta r_{ij} \ge r_i + r_j$).
