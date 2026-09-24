"""Unit tests for petrophysical analysis and diagnostics."""

import numpy as np
import pytest
from poropack.core.representations import VoxelGrid
from poropack.petrophysics.porosity import analyze_porosity
from poropack.petrophysics.morphology import (
    chord_length_distribution,
    specific_surface_area,
    two_point_correlation,
)
from poropack.petrophysics.transport import formation_factor, kozeny_carman
from poropack.petrophysics.pnm import extract_pore_network


def test_porosity_analysis():
    # Grid with a continuous pore channel through Z
    matrix = np.ones((20, 20, 20), dtype=bool)
    matrix[8:12, 8:12, :] = False  # Straight channel percolating in Z

    grid = VoxelGrid(matrix=matrix, voxel_size=1.0)
    res = analyze_porosity(grid, flow_axis=2)

    assert res.total > 0.0
    assert res.percolating_z is True
    assert res.percolating_x is False
    assert res.effective == res.total  # Channel is 100% percolating


def test_morphology_and_transport():
    # Synthetic porous media: sphere in box
    from poropack.core.representations import GrainPack
    pack = GrainPack(coords=np.array([[15.0, 15.0, 15.0]]), radii=np.array([8.0]), box_size=np.array([30.0, 30.0, 30.0]))
    grid = pack.to_voxel(voxel_size=1.0)

    # Specific surface area
    sv = specific_surface_area(grid)
    assert sv > 0.0

    # Two-point correlation
    r, s2 = two_point_correlation(grid, max_radius=8.0, n_bins=10)
    assert len(r) == 10
    # S2(0) should be close to porosity
    assert abs(s2[0] - grid.porosity()) < 0.05

    # Kozeny-Carman
    k = kozeny_carman(grid, unit="physical")
    assert k > 0.0

    # Formation factor
    F = formation_factor(grid, cementation_exponent=2.0)
    assert F >= 1.0

    # Chords
    chords = chord_length_distribution(grid, axis=0)
    assert len(chords) > 0


def test_extract_pore_network():
    # Create simple grid with two connected spherical pores
    from poropack.core.representations import GrainPack
    from poropack.generators.rsa import RSAGenerator
    gen = RSAGenerator(box_size=(40.0, 40.0, 40.0), psd=[6.0], random_state=42)
    pack = gen.generate(target_porosity=0.70, max_attempts=500)
    grid = pack.to_voxel(voxel_size=2.0)

    # Test extraction via PoreSpy SNOW2
    net = extract_pore_network(grid, accuracy="standard")
    assert "pore.coords" in net
    assert "throat.conns" in net
