# Export Formats & Simulation Interfaces (VTK, TIFF, MRST)

Type: task
Status: resolved
Blocked by: 01

## Question

What file format exporters (VTK `.vti`/`.vtp`, TIFF stacks, raw binary, CSV, HDF5, and MRST grid structures) must be implemented for seamless integration with external DRP and reservoir flow solvers?

Specifically:
1. For 3D volumetric images: implement direct export to 8-bit/16-bit TIFF image stacks, raw binary (`.raw`), and VTK structured image (`.vti`) for visualization in ParaView, Avizo, Dragonfly, and ImageJ.
2. For particle assemblies (`GrainPack`): implement CSV, Parquet, and VTK PolyData (`.vtp`) containing grain positions, radii, and contact networks.
3. For reservoir simulation: implement an export function generating MATLAB/MRST-compatible Cartesian/unstructured grid matrices (`G.cells`, `G.faces`, `rock.perm`, `rock.poro`) directly readable by MRST.

## Answer

### 1. 3D Volumetric Imaging Exporters (`VoxelGrid`)
- **TIFF Stack**: Exports 8-bit / 16-bit multi-page TIFF stacks readable by ImageJ/Fiji, Dragonfly, Avizo, and GeoDict.
- **Raw Binary (`.raw` + `.mhd` MetaImage)**: Standard format for HPC flow solvers and micro-CT tomography processing.
- **VTK Structured Grid (`.vti`)**: XML-based structured point grid for ParaView visualization, including point and cell data arrays (phase segmentation, distance-to-grain fields).

### 2. Particle Assembly Exporters (`GrainPack`)
- **CSV & Apache Parquet**: Tabular records `['x', 'y', 'z', 'r', 'phase_id', 'coordination']` for fast data analysis with Pandas/Polars.
- **VTK PolyData (`.vtp`)**: Encodes sphere centers as vertices, radii as scalar point attributes, and pairwise contacts as lines (`vtkCellArray`), enabling full 3D contact network visualization in ParaView.

### 3. Reservoir Simulation Bridge (MRST Integration)
- **MATLAB `.mat` Export (`scipy.io.savemat`)**:
  Generates MRST-ready structs matching the reservoir simulator's expected format:
  - `G`: Cartesian grid structure (`G.cartDims`, `G.cells`, `G.faces`).
  - `rock.poro`: Effective porosity per grid block.
  - `rock.perm`: Permeability tensor (mD or $\text{m}^2$) computed from local Kozeny-Carman or upscaled pore network flow.
- Seamlessly loadable in MRST via `load('poro_grid.mat')` and executable out-of-the-box.
