"""Unit tests for PeriodicSpatialGrid."""

import numpy as np
import pytest
from poropack.core.spatial import PeriodicSpatialGrid


def test_spatial_grid_insertion_and_collision():
    grid = PeriodicSpatialGrid(box_size=(100.0, 100.0, 100.0), r_max=10.0, min_throat=1.0)

    # Insert a sphere at (50, 50, 50) with radius 10
    idx = grid.insert(50.0, 50.0, 50.0, 10.0)
    assert idx == 0
    assert len(grid) == 1

    # Candidate inside or overlapping (distance = 15 < 10 + 10 + 1 = 21)
    assert grid.check_collision(50.0, 65.0, 50.0, 10.0) is True

    # Candidate touching within min_throat: distance = 20.5 < 10 + 10 + 1 = 21
    assert grid.check_collision(50.0, 70.5, 50.0, 10.0) is True

    # Candidate safely separated: distance = 25 > 21
    assert grid.check_collision(50.0, 75.0, 50.0, 10.0) is False


def test_spatial_grid_periodic_toroidal_collision():
    grid = PeriodicSpatialGrid(
        box_size=(100.0, 100.0, 100.0),
        r_max=10.0,
        periodic=(True, True, True),
    )

    # Insert sphere near x = 2.0 with radius 5.0
    grid.insert(2.0, 50.0, 50.0, 5.0)

    # Candidate at x = 98.0 with radius 5.0
    # Minimum image distance across periodic boundary is |(98 - 100) - 2| = |-4| = 4.0
    # Allowed distance is 5.0 + 5.0 = 10.0 -> Must collide!
    assert grid.check_collision(98.0, 50.0, 50.0, 5.0) is True

    # Candidate at x = 85.0: minimum distance is 17.0 > 10.0 -> Free!
    assert grid.check_collision(85.0, 50.0, 50.0, 5.0) is False
