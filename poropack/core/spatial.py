"""
poropack.core.spatial
=====================
High-performance O(1) dynamic spatial hash grid for hard-sphere collision detection
with toroidal Periodic Boundary Conditions (PBC).
"""

from __future__ import annotations

from typing import List, Optional, Tuple
import numpy as np


class PeriodicSpatialGrid:
    """
    Uniform 3D Spatial Hash Grid with Toroidal Periodic Boundary Conditions.
    Provides O(1) dynamic particle insertion and O(1) collision queries.

    Parameters
    ----------
    box_size : tuple or ndarray of float
        (Lx, Ly, Lz) dimensions of domain.
    r_max : float
        Maximum radius among all candidate spheres.
    min_throat : float, default=0.0
        Minimum clearance / throat distance required between sphere surfaces.
    periodic : tuple of bool, default=(True, True, True)
        Whether to wrap coordinates periodically along (x, y, z).
    """

    def __init__(
        self,
        box_size: Tuple[float, float, float] | np.ndarray,
        r_max: float,
        min_throat: float = 0.0,
        periodic: Tuple[bool, bool, bool] = (True, True, True),
    ):
        self.box_size = np.asarray(box_size, dtype=np.float64)
        self.r_max = float(r_max)
        self.min_throat = float(min_throat)
        self.periodic = tuple(bool(p) for p in periodic)

        # Minimum physical cell size required to guarantee that all spheres within
        # distance (r_cand + r_neighbor + min_throat) <= (2*r_max + min_throat)
        # reside in the 3x3x3 neighbor cell stencil.
        min_cell_size = 2.0 * self.r_max + self.min_throat

        # Enforce at least 3 cells along periodic dimensions to guarantee no toroidal self-intersection
        dims = np.floor(self.box_size / min_cell_size).astype(np.int32)
        for d in range(3):
            if self.periodic[d]:
                dims[d] = max(3, dims[d])
            else:
                dims[d] = max(1, dims[d])

        self.dims = dims
        self.cell_size = self.box_size / self.dims
        self.n_cells = int(np.prod(self.dims))

        # Precompute 27 neighbor stencil offsets in 3D
        self.offsets = [
            (dx, dy, dz)
            for dx in (-1, 0, 1)
            for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)
        ]

        # Flat cell-linked bucket list: cell_idx -> list of sphere indices
        self.grid: List[List[int]] = [[] for _ in range(self.n_cells)]
        self.coords: List[Tuple[float, float, float]] = []
        self.radii: List[float] = []

    def __len__(self) -> int:
        return len(self.coords)

    def _get_cell_idx(self, x: float, y: float, z: float) -> Tuple[int, int, int, int]:
        """Convert continuous coordinate to 3D and 1D flat cell index."""
        cx = int(x / self.cell_size[0])
        cy = int(y / self.cell_size[1])
        cz = int(z / self.cell_size[2])

        # Floored modulo for periodic wrap-around
        if self.periodic[0]:
            cx = cx % self.dims[0]
        else:
            cx = max(0, min(self.dims[0] - 1, cx))

        if self.periodic[1]:
            cy = cy % self.dims[1]
        else:
            cy = max(0, min(self.dims[1] - 1, cy))

        if self.periodic[2]:
            cz = cz % self.dims[2]
        else:
            cz = max(0, min(self.dims[2] - 1, cz))

        flat_idx = (cx * self.dims[1] + cy) * self.dims[2] + cz
        return cx, cy, cz, flat_idx

    def check_collision(
        self,
        x: float,
        y: float,
        z: float,
        r_cand: float,
        ignore_idx: Optional[int] = None,
    ) -> bool:
        """
        Check if candidate sphere (x, y, z, r_cand) collides with any existing sphere.
        Runs in O(1) time by querying only the 27 neighboring cells.

        Returns True if collision / overlap occurs, False if space is available.
        """
        cx, cy, cz, _ = self._get_cell_idx(x, y, z)
        Lx, Ly, Lz = self.box_size
        half_Lx, half_Ly, half_Lz = 0.5 * Lx, 0.5 * Ly, 0.5 * Lz
        pbc_x, pbc_y, pbc_z = self.periodic

        for dx, dy, dz in self.offsets:
            ncx = cx + dx
            ncy = cy + dy
            ncz = cz + dz

            # Check boundary conditions per axis
            if pbc_x:
                ncx = ncx % self.dims[0]
            elif ncx < 0 or ncx >= self.dims[0]:
                continue

            if pbc_y:
                ncy = ncy % self.dims[1]
            elif ncy < 0 or ncy >= self.dims[1]:
                continue

            if pbc_z:
                ncz = ncz % self.dims[2]
            elif ncz < 0 or ncz >= self.dims[2]:
                continue

            cell_flat = (ncx * self.dims[1] + ncy) * self.dims[2] + ncz

            for s_idx in self.grid[cell_flat]:
                if ignore_idx is not None and s_idx == ignore_idx:
                    continue

                sx, sy, sz = self.coords[s_idx]
                limit_dist = r_cand + self.radii[s_idx] + self.min_throat

                # Minimum image convention with early axis-aligned coordinate pruning
                dist_x = abs(x - sx)
                if pbc_x and dist_x > half_Lx:
                    dist_x = Lx - dist_x
                if dist_x >= limit_dist:
                    continue

                dist_y = abs(y - sy)
                if pbc_y and dist_y > half_Ly:
                    dist_y = Ly - dist_y
                if dist_y >= limit_dist:
                    continue

                dist_z = abs(z - sz)
                if pbc_z and dist_z > half_Lz:
                    dist_z = Lz - dist_z
                if dist_z >= limit_dist:
                    continue

                # Full 3D Euclidean distance check
                if (dist_x * dist_x + dist_y * dist_y + dist_z * dist_z) < (limit_dist * limit_dist):
                    return True  # Collision detected

        return False  # Free space

    def insert(self, x: float, y: float, z: float, r: float) -> int:
        """
        Dynamically insert an accepted sphere into the spatial hash grid in O(1) time.

        Returns
        -------
        int
            Index of the newly inserted sphere.
        """
        # Periodic coordinate wrapping
        if self.periodic[0]:
            x = x % self.box_size[0]
        if self.periodic[1]:
            y = y % self.box_size[1]
        if self.periodic[2]:
            z = z % self.box_size[2]

        _, _, _, cell_flat = self._get_cell_idx(x, y, z)
        idx = len(self.coords)
        self.coords.append((x, y, z))
        self.radii.append(r)
        self.grid[cell_flat].append(idx)
        return idx

    def get_arrays(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return coordinates and radii as NumPy arrays."""
        if len(self.coords) == 0:
            return np.empty((0, 3), dtype=np.float64), np.empty(0, dtype=np.float64)
        return np.array(self.coords, dtype=np.float64), np.array(self.radii, dtype=np.float64)
