"""
poropack.petrophysics.porosity
==============================
Porosity partitioning: Total, Effective (percolating), and Isolated porosity analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
from scipy import ndimage

from poropack.core.representations import VoxelGrid


@dataclass
class PorosityResult:
    """Detailed porosity breakdown."""
    total: float
    effective: float
    isolated: float
    percolating_x: bool
    percolating_y: bool
    percolating_z: bool

    @property
    def total_porosity(self) -> float:
        return self.total

    @property
    def effective_porosity(self) -> float:
        return self.effective

    @property
    def isolated_porosity(self) -> float:
        return self.isolated

    def summary(self) -> str:
        return (
            f"Total Porosity:     {self.total:.4f} ({self.total*100:.2f}%)\n"
            f"Effective Porosity: {self.effective:.4f} ({self.effective*100:.2f}%)\n"
            f"Isolated Porosity:  {self.isolated:.4f} ({self.isolated*100:.2f}%)\n"
            f"Percolating Axes:   X={self.percolating_x}, Y={self.percolating_y}, Z={self.percolating_z}"
        )


def analyze_porosity(
    grid: VoxelGrid,
    connectivity: int = 6,
    flow_axis: Optional[int] = None,
) -> PorosityResult:
    """
    Partition total porosity into percolating (effective) and dead-end (isolated) components.

    Parameters
    ----------
    grid : VoxelGrid
        Input porous media.
    connectivity : {6, 26}, default=6
        Neighbor connectivity for 3D void clustering.
    flow_axis : int in {0, 1, 2}, optional
        If specified, effective porosity is strictly defined as percolating along this axis.
        If None, defined as percolating across any of the orthogonal axes.

    Returns
    -------
    PorosityResult
    """
    pore_mask = ~grid.matrix
    total_voxels = grid.matrix.size
    total_pore_voxels = int(np.count_nonzero(pore_mask))

    if total_pore_voxels == 0:
        return PorosityResult(0.0, 0.0, 0.0, False, False, False)

    # 3D Structuring element for clustering
    if connectivity == 6:
        struct = ndimage.generate_binary_structure(3, 1)
    else:
        struct = ndimage.generate_binary_structure(3, 3)

    labeled_pores, num_features = ndimage.label(pore_mask, structure=struct)

    Nx, Ny, Nz = grid.shape

    # Check percolation across orthogonal boundaries
    # Face clusters at 0 and L
    x0_labels = set(np.unique(labeled_pores[0, :, :])) - {0}
    x1_labels = set(np.unique(labeled_pores[Nx - 1, :, :])) - {0}
    perc_x_labels = x0_labels.intersection(x1_labels)
    percolates_x = len(perc_x_labels) > 0

    y0_labels = set(np.unique(labeled_pores[:, 0, :])) - {0}
    y1_labels = set(np.unique(labeled_pores[:, Ny - 1, :])) - {0}
    perc_y_labels = y0_labels.intersection(y1_labels)
    percolates_y = len(perc_y_labels) > 0

    z0_labels = set(np.unique(labeled_pores[:, :, 0])) - {0}
    z1_labels = set(np.unique(labeled_pores[:, :, Nz - 1])) - {0}
    perc_z_labels = z0_labels.intersection(z1_labels)
    percolates_z = len(perc_z_labels) > 0

    if flow_axis == 0:
        effective_labels = perc_x_labels
    elif flow_axis == 1:
        effective_labels = perc_y_labels
    elif flow_axis == 2:
        effective_labels = perc_z_labels
    else:
        effective_labels = perc_x_labels.union(perc_y_labels).union(perc_z_labels)

    if len(effective_labels) > 0:
        # Sum voxels belonging to percolating clusters
        effective_mask = np.isin(labeled_pores, list(effective_labels))
        effective_voxels = int(np.count_nonzero(effective_mask))
    else:
        effective_voxels = 0

    isolated_voxels = total_pore_voxels - effective_voxels

    total_phi = total_pore_voxels / total_voxels
    effective_phi = effective_voxels / total_voxels
    isolated_phi = isolated_voxels / total_voxels

    return PorosityResult(
        total=total_phi,
        effective=effective_phi,
        isolated=isolated_phi,
        percolating_x=percolates_x,
        percolating_y=percolates_y,
        percolating_z=percolates_z,
    )
