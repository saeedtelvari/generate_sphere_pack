"""
poropack.io.exporters
=====================
Export formats for 3D digital rock visualization and reservoir simulation:
- 3D TIFF image stacks (ImageJ, Dragonfly, Avizo)
- VTK Structured Image (.vti) and PolyData (.vtp) (ParaView)
- Raw binary (.raw)
- MATLAB MRST grid and rock matrix structs (.mat)
- Tabular CSV / Parquet
"""

from __future__ import annotations

import os
from typing import Any, Dict, Optional, Tuple, Union
import numpy as np
from scipy import io as spio

from poropack.core.representations import GrainPack, VoxelGrid


def export_tiff(
    grid: VoxelGrid,
    filepath: str,
    dtype: str = "uint8",
) -> None:
    """
    Export 3D binary/segmented voxel grid as a multi-page TIFF stack.

    Parameters
    ----------
    grid : VoxelGrid
        Porous media grid.
    filepath : str
        Target .tif or .tiff file path.
    dtype : {'uint8', 'uint16'}, default='uint8'
        Data type. For uint8: 0 = pore, 255 = solid.
    """
    try:
        import tifffile
    except ImportError as e:
        raise ImportError("TIFF export requires 'tifffile'. Install it via 'pip install tifffile'.") from e

    # Convention: 0 = pore space, 255 = solid grain
    if dtype == "uint8":
        data = (grid.matrix.astype(np.uint8)) * 255
    elif dtype == "uint16":
        data = (grid.matrix.astype(np.uint16)) * 65535
    else:
        data = grid.matrix.astype(np.uint8) * 255

    # Move Z to first dimension for standard slice stacking: (Nz, Ny, Nx)
    data_transposed = np.moveaxis(data, 2, 0)
    tifffile.imwrite(filepath, data_transposed, photometric="minisblack")


def export_raw(
    grid: VoxelGrid,
    filepath: str,
) -> None:
    """
    Export 3D binary grid as raw binary file (.raw) with companion ASCII header (.mhd).
    """
    raw_path = filepath if filepath.endswith(".raw") else filepath + ".raw"
    mhd_path = os.path.splitext(raw_path)[0] + ".mhd"

    # Write binary bytes
    data_bytes = (grid.matrix.astype(np.uint8) * 255).tobytes()
    with open(raw_path, "wb") as f:
        f.write(data_bytes)

    # Write MetaImage header (.mhd) for ParaView / ImageJ
    Nx, Ny, Nz = grid.shape
    h = grid.voxel_size
    raw_filename = os.path.basename(raw_path)

    mhd_content = (
        "ObjectType = Image\n"
        "NDims = 3\n"
        "BinaryData = True\n"
        "BinaryDataByteOrderMSB = False\n"
        f"DimSize = {Nx} {Ny} {Nz}\n"
        f"ElementSpacing = {h} {h} {h}\n"
        "ElementType = MET_UCHAR\n"
        f"ElementDataFile = {raw_filename}\n"
    )
    with open(mhd_path, "w", encoding="utf-8") as f:
        f.write(mhd_content)


def export_vtk_vti(
    grid: VoxelGrid,
    filepath: str,
    extra_arrays: Optional[Dict[str, np.ndarray]] = None,
) -> None:
    """
    Export VoxelGrid as a VTK XML Structured Image (.vti) for ParaView.
    Uses PyVista or VTK if available, with a clean pure-Python XML fallback.
    """
    if not filepath.endswith(".vti"):
        filepath += ".vti"

    Nx, Ny, Nz = grid.shape
    h = grid.voxel_size

    try:
        import pyvista as pv
        # Construct PyVista UniformGrid (ImageData)
        pv_grid = pv.ImageData(
            dimensions=(Nx + 1, Ny + 1, Nz + 1),
            spacing=(h, h, h),
            origin=tuple(grid.origin),
        )
        # Assign cell data (voxels)
        # PyVista/VTK cell order is Fortran order
        pv_grid.cell_data["Phase"] = grid.matrix.flatten(order="F").astype(np.uint8)
        if extra_arrays:
            for name, arr in extra_arrays.items():
                pv_grid.cell_data[name] = arr.flatten(order="F")
        pv_grid.save(filepath)

    except ImportError:
        # Pure-Python XML VTI fallback
        ox, oy, oz = grid.origin
        matrix_flat = (grid.matrix.astype(np.uint8)).flatten(order="F")
        data_str = " ".join(map(str, matrix_flat))

        vti_xml = (
            '<?xml version="1.0"?>\n'
            '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n'
            f'  <ImageData WholeExtent="0 {Nx} 0 {Ny} 0 {Nz}" Origin="{ox} {oy} {oz}" Spacing="{h} {h} {h}">\n'
            f'    <Piece Extent="0 {Nx} 0 {Ny} 0 {Nz}">\n'
            '      <CellData Scalars="Phase">\n'
            '        <DataArray type="UInt8" Name="Phase" format="ascii">\n'
            f"          {data_str}\n"
            "        </DataArray>\n"
            "      </CellData>\n"
            "    </Piece>\n"
            "  </ImageData>\n"
            "</VTKFile>\n"
        )
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(vti_xml)


def export_vtk_vtp(
    pack: GrainPack,
    filepath: str,
) -> None:
    """
    Export GrainPack as a VTK XML PolyData (.vtp) for 3D sphere glyph and contact network visualization.
    """
    if not filepath.endswith(".vtp"):
        filepath += ".vtp"

    try:
        import pyvista as pv
        pv_cloud = pv.PolyData(pack.coords)
        pv_cloud.point_data["Radius"] = pack.radii
        for k, v in pack.attributes.items():
            if isinstance(v, (list, np.ndarray)) and len(v) == len(pack):
                pv_cloud.point_data[k] = np.asarray(v)
        pv_cloud.save(filepath)

    except ImportError:
        # Fallback CSV
        pack.to_csv(filepath.replace(".vtp", ".csv"))


def export_mrst_mat(
    grid: VoxelGrid,
    filepath: str,
    perm_calc: bool = True,
    permeability_mD: Optional[float] = None,
) -> None:
    """
    Export VoxelGrid as a MATLAB struct (.mat) formatted for the MATLAB Reservoir Simulation Toolbox (MRST).

    Generates:
    - G: Cartesian grid structure (G.cartDims, G.cells.num, G.faces.num)
    - rock.poro: Effective/local porosity values per cell
    - rock.perm: Permeability tensor / values (in milliDarcy) matching MRST conventions
    """
    if not filepath.endswith(".mat"):
        filepath += ".mat"

    from poropack.petrophysics.porosity import analyze_porosity
    from poropack.petrophysics.transport import kozeny_carman

    poro_res = analyze_porosity(grid)
    phi_val = poro_res.effective if poro_res.effective > 0 else poro_res.total
    if permeability_mD is not None:
        k_md = float(permeability_mD)
    else:
        k_md = kozeny_carman(grid, unit="millidarcy") if perm_calc else 100.0

    Nx, Ny, Nz = grid.shape
    h = grid.voxel_size
    n_cells = Nx * Ny * Nz

    # Create voxel-level or coarse rock structs
    # Porosity field: 0.0 for solid voxels, 1.0 for pore voxels
    poro_array = (~grid.matrix).astype(np.float64).flatten(order="F")

    # MRST Cartesian grid struct
    g_struct = {
        "cartDims": np.array([Nx, Ny, Nz], dtype=np.float64),
        "cells": {
            "num": float(n_cells),
        },
        "faces": {
            "num": float((Nx + 1) * Ny * Nz + Nx * (Ny + 1) * Nz + Nx * Ny * (Nz + 1)),
        },
        "voxel_size": float(h),
        "dimensions": grid.physical_size,
    }

    rock_struct = {
        "poro": poro_array[:, None],  # (N_cells, 1)
        "perm": np.full((n_cells, 1), k_md, dtype=np.float64),  # scalar perm per cell
        "mean_poro": float(phi_val),
        "mean_perm_md": float(k_md),
    }

    mat_dict = {
        "G": g_struct,
        "rock": rock_struct,
        "description": "poropack generated porous media for MRST simulation",
    }

    spio.savemat(filepath, mat_dict, do_compression=True)


def export_csv(
    pack: GrainPack,
    filepath: str,
) -> None:
    """
    Export GrainPack coordinates and radii as a tabular CSV file.
    """
    pack.to_csv(filepath)


def export_parquet(
    pack: GrainPack,
    filepath: str,
) -> None:
    """
    Export GrainPack coordinates and radii as a columnar Apache Parquet file.
    """
    if not filepath.endswith(".parquet"):
        filepath += ".parquet"
    pack.to_dataframe().to_parquet(filepath, index=False)
