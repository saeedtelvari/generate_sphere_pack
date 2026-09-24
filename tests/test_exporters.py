"""Unit tests for I/O exporters and visualization utilities."""

import os
import tempfile
import numpy as np
import pytest
from poropack.core.representations import GrainPack, VoxelGrid
from poropack.io.exporters import (
    export_csv,
    export_mrst_mat,
    export_parquet,
    export_raw,
    export_tiff,
    export_vtk_vti,
    export_vtk_vtp,
)
from poropack.io.visualization import (
    plot_orthogonal_slices,
    plot_slice,
    plot_xray,
)


def test_exporters():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test GrainPack and VoxelGrid
        coords = np.array([[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]])
        radii = np.array([4.0, 5.0])
        box_size = np.array([30.0, 30.0, 30.0])
        pack = GrainPack(coords=coords, radii=radii, box_size=box_size)
        grid = pack.to_voxel(voxel_size=1.0)

        # 1. TIFF
        tiff_path = os.path.join(tmpdir, "test.tif")
        export_tiff(grid, tiff_path)
        assert os.path.exists(tiff_path)
        assert os.path.getsize(tiff_path) > 0

        # 2. Raw + MHD
        raw_path = os.path.join(tmpdir, "test.raw")
        export_raw(grid, raw_path)
        assert os.path.exists(raw_path)
        assert os.path.exists(os.path.join(tmpdir, "test.mhd"))

        # 3. VTK VTI
        vti_path = os.path.join(tmpdir, "test.vti")
        export_vtk_vti(grid, vti_path)
        assert os.path.exists(vti_path)

        # 4. VTK VTP
        vtp_path = os.path.join(tmpdir, "test.vtp")
        export_vtk_vtp(pack, vtp_path)
        assert os.path.exists(vtp_path)

        # 5. CSV & Parquet
        csv_path = os.path.join(tmpdir, "test.csv")
        export_csv(pack, csv_path)
        assert os.path.exists(csv_path)

        parquet_path = os.path.join(tmpdir, "test.parquet")
        export_parquet(pack, parquet_path)
        assert os.path.exists(parquet_path)
        assert os.path.getsize(parquet_path) > 0

        # 6. MRST .mat
        mat_path = os.path.join(tmpdir, "test_mrst.mat")
        export_mrst_mat(grid, mat_path)
        assert os.path.exists(mat_path)


def test_visualizations():
    matrix = np.zeros((20, 20, 20), dtype=bool)
    matrix[5:15, 5:15, 5:15] = True
    grid = VoxelGrid(matrix=matrix, voxel_size=1.0)

    fig1, ax1 = plot_xray(grid, axis=2, show=False)
    assert fig1 is not None

    fig2, ax2 = plot_slice(grid, slice_idx=10, axis=2, show=False)
    assert fig2 is not None

    fig3, axes3 = plot_orthogonal_slices(grid, show=False)
    assert len(axes3) == 3
