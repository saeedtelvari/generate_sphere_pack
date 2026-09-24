"""Unit tests for poropack.core representations and native rasterizer."""

import numpy as np
import pytest
from poropack.core.representations import GrainPack, VoxelGrid, rasterize_spheres


def test_grain_pack_creation_and_properties():
    coords = np.array([[10.0, 10.0, 10.0], [30.0, 30.0, 30.0]])
    radii = np.array([5.0, 8.0])
    box_size = np.array([50.0, 50.0, 50.0])

    pack = GrainPack(coords=coords, radii=radii, box_size=box_size)
    assert len(pack) == 2
    assert pack.volume == 125000.0

    expected_solid = 4.0 / 3.0 * np.pi * (5.0 ** 3 + 8.0 ** 3)
    assert np.isclose(pack.solid_volume, expected_solid)
    assert np.isclose(pack.analytical_porosity(), 1.0 - (expected_solid / 125000.0))

    df = pack.to_dataframe()
    assert list(df.columns) == ["X", "Y", "Z", "R"]
    assert len(df) == 2

    pack_copy = pack.copy()
    assert np.array_equal(pack_copy.coords, pack.coords)
    assert pack_copy is not pack


def test_voxel_grid_creation_and_porosity():
    matrix = np.zeros((20, 20, 20), dtype=bool)
    matrix[5:15, 5:15, 5:15] = True  # 10x10x10 = 1000 solid voxels out of 8000

    grid = VoxelGrid(matrix=matrix, voxel_size=0.5)
    assert grid.shape == (20, 20, 20)
    assert np.isclose(grid.volume, 20 * 0.5 * 20 * 0.5 * 20 * 0.5)
    assert np.isclose(grid.porosity(), 1.0 - (1000.0 / 8000.0))
    assert np.isclose(grid.solid_fraction(), 1000.0 / 8000.0)
    assert grid.pore_mask.shape == (20, 20, 20)


def test_rasterize_spheres_periodic():
    coords = np.array([[0.0, 0.0, 0.0]])  # Centered right at corner
    radii = np.array([4.0])
    box_size = np.array([20.0, 20.0, 20.0])
    voxel_size = 1.0

    matrix = rasterize_spheres(coords, radii, box_size, voxel_size, periodic=(True, True, True))
    assert matrix.shape == (20, 20, 20)

    # In periodic mode, sphere at (0,0,0) should wrap to all 8 corners:
    # (0,0,0), (19,0,0), (0,19,0), (19,19,0), etc.
    assert matrix[0, 0, 0] is np.True_ or matrix[0, 0, 0] == True
    assert matrix[19, 0, 0] is np.True_ or matrix[19, 0, 0] == True
    assert matrix[0, 19, 0] is np.True_ or matrix[0, 19, 0] == True
    assert matrix[19, 19, 19] is np.True_ or matrix[19, 19, 19] == True

    # Center of domain should be empty pore space
    assert matrix[10, 10, 10] == False
