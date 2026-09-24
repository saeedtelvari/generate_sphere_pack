"""Unit tests for compaction and cementation diagenesis."""

import numpy as np
import pytest
from poropack.core.representations import GrainPack, VoxelGrid
from poropack.diagenesis.compaction import apply_compaction
from poropack.diagenesis.cementation import apply_cementation, apply_dissolution


def test_compaction():
    coords = np.array([[10.0, 10.0, 20.0], [10.0, 10.0, 40.0]])
    radii = np.array([5.0, 5.0])
    box_size = np.array([50.0, 50.0, 60.0])

    pack = GrainPack(coords=coords, radii=radii, box_size=box_size)
    compacted = apply_compaction(pack, vertical_strain=0.10)

    # Z-coordinates should be scaled by (1 - 0.10) = 0.90
    assert np.isclose(compacted.coords[0, 2], 18.0)
    assert np.isclose(compacted.coords[1, 2], 36.0)
    assert np.isclose(compacted.box_size[2], 54.0)
    # Radii should remain constant
    assert np.array_equal(compacted.radii, pack.radii)


def test_cementation_and_dissolution():
    # Create empty pore space with single solid block in center
    matrix = np.zeros((30, 30, 30), dtype=bool)
    matrix[10:20, 10:20, 10:20] = True  # initial 1000 solid voxels
    grid = VoxelGrid(matrix=matrix, voxel_size=1.0)
    initial_porosity = grid.porosity()

    # Apply 5% cementation
    cemented = apply_cementation(grid, cement_fraction=0.05)
    assert cemented.porosity() < initial_porosity
    assert abs((initial_porosity - cemented.porosity()) - 0.05) < 0.01

    # Apply dissolution on cemented rock
    dissolved = apply_dissolution(cemented, porosity_increase=0.03, mode="pore_surface")
    assert dissolved.porosity() > cemented.porosity()

    # Apply grain_cores dissolution and dissolve alias
    from poropack.diagenesis.cementation import dissolve
    dissolved_cores = dissolve(cemented, porosity_increase=0.02, mode="grain_cores")
    assert dissolved_cores.porosity() > cemented.porosity()
