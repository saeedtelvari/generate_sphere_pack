"""
poropack.core.representations
=============================
Lagrangian (GrainPack) and Eulerian (VoxelGrid) porous media representations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple, Union
import numpy as np
import pandas as pd


def rasterize_spheres(
    coords: np.ndarray,
    radii: np.ndarray,
    box_size: Tuple[float, float, float] | np.ndarray,
    voxel_size: float,
    periodic: Tuple[bool, bool, bool] = (True, True, True),
) -> np.ndarray:
    """
    Fast sub-volume bounding-box rasterization of 3D spheres into a binary voxel matrix.
    Supports native toroidal periodic boundary wrap-around in x, y, and z.

    Parameters
    ----------
    coords : np.ndarray
        (N, 3) sphere centroid coordinates.
    radii : np.ndarray
        (N,) sphere radii.
    box_size : tuple or ndarray of float
        (Lx, Ly, Lz) dimensions of the bounding domain.
    voxel_size : float
        Physical isotropic spacing per voxel.
    periodic : tuple of bool
        Whether to wrap spheres periodically across (x, y, z) boundaries.

    Returns
    -------
    matrix : np.ndarray of bool
        3D binary array where True (1) is solid grain and False (0) is pore space.
    """
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    box_size = np.asarray(box_size, dtype=np.float64)

    # Grid dimensions (number of voxels along each axis)
    dims = np.maximum(1, np.round(box_size / voxel_size).astype(np.int32))
    Nx, Ny, Nz = int(dims[0]), int(dims[1]), int(dims[2])
    matrix = np.zeros((Nx, Ny, Nz), dtype=bool)

    if len(coords) == 0:
        return matrix

    pbc_x, pbc_y, pbc_z = periodic

    for (x, y, z), r in zip(coords, radii):
        # Determine bounding box in voxel indices
        r_vox = r / voxel_size
        r_vox_ceil = int(np.ceil(r_vox))

        cx = x / voxel_size
        cy = y / voxel_size
        cz = z / voxel_size

        i_min = int(np.floor(cx - r_vox_ceil))
        i_max = int(np.ceil(cx + r_vox_ceil))
        j_min = int(np.floor(cy - r_vox_ceil))
        j_max = int(np.ceil(cy + r_vox_ceil))
        k_min = int(np.floor(cz - r_vox_ceil))
        k_max = int(np.ceil(cz + r_vox_ceil))

        # Local integer coordinates
        ii = np.arange(i_min, i_max + 1)
        jj = np.arange(j_min, j_max + 1)
        kk = np.arange(k_min, k_max + 1)

        # Distances to sphere center
        dx = (ii + 0.5) - cx
        dy = (jj + 0.5) - cy
        dz = (kk + 0.5) - cz

        # 3D distance squared via outer addition
        dist_sq_xy = dx[:, None] ** 2 + dy[None, :] ** 2
        # (len(ii), len(jj), len(kk))
        dist_sq = dist_sq_xy[:, :, None] + (dz[None, None, :] ** 2)
        inside = dist_sq <= (r_vox ** 2)

        if not np.any(inside):
            continue

        # Get relative indices where points are inside sphere
        local_i, local_j, local_k = np.where(inside)
        global_i = ii[local_i]
        global_j = jj[local_j]
        global_k = kk[local_k]

        # Handle Periodic vs Non-Periodic boundary clipping
        if pbc_x:
            global_i = global_i % Nx
        else:
            mask_x = (global_i >= 0) & (global_i < Nx)
            global_i = global_i[mask_x]
            global_j = global_j[mask_x]
            global_k = global_k[mask_x]

        if pbc_y:
            global_j = global_j % Ny
        else:
            mask_y = (global_j >= 0) & (global_j < Ny)
            global_i = global_i[mask_y]
            global_j = global_j[mask_y]
            global_k = global_k[mask_y]

        if pbc_z:
            global_k = global_k % Nz
        else:
            mask_z = (global_k >= 0) & (global_k < Nz)
            global_i = global_i[mask_z]
            global_j = global_j[mask_z]
            global_k = global_k[mask_z]

        matrix[global_i, global_j, global_k] = True

    return matrix


@dataclass
class GrainPack:
    """
    Lagrangian discrete representation of 3D spherical grain assemblies.

    Attributes
    ----------
    coords : np.ndarray
        (N, 3) float64 array of grain centroid positions.
    radii : np.ndarray
        (N,) float64 array of grain radii.
    box_size : np.ndarray
        (3,) float64 dimensions [Lx, Ly, Lz].
    periodic : tuple of bool
        Periodic boundary indicators for (x, y, z).
    attributes : dict of str -> np.ndarray
        Optional metadata (e.g. contact indices, mineralogy IDs).
    """

    coords: np.ndarray
    radii: np.ndarray
    box_size: np.ndarray
    periodic: Tuple[bool, bool, bool] = (True, True, True)
    attributes: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.coords = np.asarray(self.coords, dtype=np.float64)
        if self.coords.ndim == 1 and len(self.coords) == 0:
            self.coords = np.empty((0, 3), dtype=np.float64)
        elif self.coords.ndim != 2 or self.coords.shape[1] != 3:
            raise ValueError(f"coords must be shape (N, 3), got {self.coords.shape}")

        self.radii = np.asarray(self.radii, dtype=np.float64)
        if self.radii.ndim != 1 or len(self.radii) != len(self.coords):
            raise ValueError(f"radii must be 1D array of length {len(self.coords)}")

        self.box_size = np.asarray(self.box_size, dtype=np.float64)
        if self.box_size.shape != (3,):
            raise ValueError(f"box_size must be 3-element array, got {self.box_size.shape}")

    def __len__(self) -> int:
        return len(self.coords)

    @property
    def volume(self) -> float:
        """Total bounding volume Lx * Ly * Lz."""
        return float(np.prod(self.box_size))

    @property
    def solid_volume(self) -> float:
        """Sum of individual sphere volumes: sum(4/3 * pi * r^3)."""
        return float(np.sum(4.0 / 3.0 * np.pi * (self.radii ** 3)))

    def analytical_porosity(self) -> float:
        """
        Analytical porosity assuming zero overlap: 1 - solid_volume / box_volume.
        Note: For overlapping packs, use .to_voxel().porosity().
        """
        return 1.0 - (self.solid_volume / self.volume)

    def wrap_coordinates(self) -> GrainPack:
        """Returns a copy of GrainPack with coordinates modulo box_size where periodic."""
        new_coords = self.coords.copy()
        for d in range(3):
            if self.periodic[d]:
                new_coords[:, d] = np.mod(new_coords[:, d], self.box_size[d])
        return GrainPack(
            coords=new_coords,
            radii=self.radii.copy(),
            box_size=self.box_size.copy(),
            periodic=self.periodic,
            attributes={k: v.copy() if hasattr(v, "copy") else v for k, v in self.attributes.items()},
        )

    def copy(self) -> GrainPack:
        """Create a deep copy of GrainPack."""
        return GrainPack(
            coords=self.coords.copy(),
            radii=self.radii.copy(),
            box_size=self.box_size.copy(),
            periodic=self.periodic,
            attributes={k: v.copy() if hasattr(v, "copy") else v for k, v in self.attributes.items()},
        )

    def to_voxel(self, voxel_size: float) -> VoxelGrid:
        """
        Rasterize the sphere pack into a discrete Eulerian VoxelGrid.

        Parameters
        ----------
        voxel_size : float
            Isotropic voxel spacing.

        Returns
        -------
        VoxelGrid
            Eulerian representation of the porous media.
        """
        matrix = rasterize_spheres(
            self.coords,
            self.radii,
            self.box_size,
            voxel_size=voxel_size,
            periodic=self.periodic,
        )
        return VoxelGrid(
            matrix=matrix,
            voxel_size=voxel_size,
            origin=np.zeros(3, dtype=np.float64),
            periodic=self.periodic,
        )

    def rasterize(self, voxel_size: float) -> VoxelGrid:
        """Alias for to_voxel: Rasterize the sphere pack into a discrete Eulerian VoxelGrid."""
        return self.to_voxel(voxel_size)

    def to_dataframe(self) -> pd.DataFrame:
        """Export grain coordinates and radii as a pandas DataFrame."""
        df = pd.DataFrame(
            {
                "X": self.coords[:, 0],
                "Y": self.coords[:, 1],
                "Z": self.coords[:, 2],
                "R": self.radii,
            }
        )
        for k, v in self.attributes.items():
            if isinstance(v, (list, np.ndarray)) and len(v) == len(self):
                df[k] = v
        return df

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        box_size: Optional[Tuple[float, float, float] | np.ndarray] = None,
        periodic: Tuple[bool, bool, bool] = (True, True, True),
    ) -> GrainPack:
        """Construct GrainPack from a DataFrame containing X, Y, Z, R columns."""
        x_col = "X" if "X" in df.columns else "x"
        y_col = "Y" if "Y" in df.columns else "y"
        z_col = "Z" if "Z" in df.columns else "z"
        r_col = "R" if "R" in df.columns else "r"

        coords = df[[x_col, y_col, z_col]].to_numpy(dtype=np.float64)
        radii = df[r_col].to_numpy(dtype=np.float64)

        if box_size is None:
            # Estimate box size based on max coordinates + radii
            max_c = np.max(coords + radii[:, None], axis=0)
            box_size = np.ceil(max_c)

        # Retain extra columns as attributes
        extra_cols = [c for c in df.columns if c not in (x_col, y_col, z_col, r_col)]
        attrs = {c: df[c].to_numpy() for c in extra_cols}

        return cls(coords=coords, radii=radii, box_size=box_size, periodic=periodic, attributes=attrs)

    def to_csv(self, filepath: str) -> None:
        """Save sphere coordinates and radii to CSV."""
        self.to_dataframe().to_csv(filepath, index=False)


@dataclass
class VoxelGrid:
    """
    Eulerian 3D discrete representation of porous media.

    Attributes
    ----------
    matrix : np.ndarray of bool
        3D array: True (1) = solid grain / cement, False (0) = pore space.
    voxel_size : float
        Isotropic voxel dimension h.
    origin : np.ndarray
        (3,) physical coordinates of origin [x0, y0, z0].
    periodic : tuple of bool
        Periodic boundary indicators for (x, y, z).
    """

    matrix: np.ndarray
    voxel_size: float = 1.0
    origin: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    periodic: Tuple[bool, bool, bool] = (True, True, True)

    def __post_init__(self):
        mat = np.asarray(self.matrix)
        if mat.dtype == bool:
            self.matrix = mat
        elif np.issubdtype(mat.dtype, np.integer):
            self.matrix = mat.astype(np.uint8)
        else:
            self.matrix = mat.astype(bool)

        if self.matrix.ndim != 3:
            raise ValueError(f"matrix must be 3D array, got ndim={self.matrix.ndim}")
        self.voxel_size = float(self.voxel_size)
        if self.voxel_size <= 0:
            raise ValueError("voxel_size must be strictly positive")
        self.origin = np.asarray(self.origin, dtype=np.float64)

    @property
    def voxels(self) -> np.ndarray:
        """Alias for matrix matching CONTEXT.md domain standards."""
        return self.matrix

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Grid dimensions (Nx, Ny, Nz)."""
        return self.matrix.shape

    @property
    def physical_size(self) -> np.ndarray:
        """Physical box dimensions Lx, Ly, Lz = Nx*h, Ny*h, Nz*h."""
        return np.array(self.matrix.shape, dtype=np.float64) * self.voxel_size

    @property
    def volume(self) -> float:
        """Total physical volume of the voxel domain."""
        return float(np.prod(self.physical_size))

    @property
    def pore_mask(self) -> np.ndarray:
        """Boolean mask where True = pore space (phase 0), False = solid."""
        if self.matrix.dtype == bool:
            return ~self.matrix
        return self.matrix == 0

    @property
    def solid_mask(self) -> np.ndarray:
        """Boolean mask where True = solid phase (phase > 0), False = pore."""
        if self.matrix.dtype == bool:
            return self.matrix
        return self.matrix > 0

    def porosity(self) -> float:
        """Exact numerical porosity: fraction of voxels occupied by pore space."""
        return float(np.count_nonzero(self.pore_mask)) / self.matrix.size

    def solid_fraction(self) -> float:
        """Exact numerical solid volume fraction: 1 - porosity."""
        return float(np.count_nonzero(self.solid_mask)) / self.matrix.size

    def to_porespy(self) -> np.ndarray:
        """
        Return the boolean 3D array in PoreSpy's convention.
        In PoreSpy, True usually denotes the pore phase and False denotes solid.
        """
        return self.pore_mask

    def to_openpnm(self) -> Any:
        """
        Export network representation to OpenPNM using PoreSpy SNOW2.
        """
        from poropack.petrophysics.pnm import extract_pore_network
        return extract_pore_network(self)

    def copy(self) -> VoxelGrid:
        """Deep copy of VoxelGrid."""
        return VoxelGrid(
            matrix=self.matrix.copy(),
            voxel_size=self.voxel_size,
            origin=self.origin.copy(),
            periodic=self.periodic,
        )
