"""
poropack.io.visualization
=========================
2D and 3D visualization utilities:
- plot_xray: Radiographic projection
- plot_slice: 2D cross-sectional slice
- plot_orthogonal_slices: Tri-planar orthogonal view (XY, XZ, YZ)
"""

from __future__ import annotations

from typing import Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np

from poropack.core.representations import VoxelGrid


def plot_xray(
    grid: VoxelGrid,
    axis: int = 2,
    depth: Optional[int] = None,
    cmap: str = "bone",
    ax: Optional[plt.Axes] = None,
    show: bool = False,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Generate an X-ray radiographic projection of the porous structure.

    Parameters
    ----------
    grid : VoxelGrid
        Porous media volume.
    axis : {0, 1, 2}, default=2
        Axis of projection.
    depth : int, optional
        Maximum voxel depth to penetrate. Defaults to full domain.
    cmap : str, default='bone'
        Matplotlib colormap.
    ax : plt.Axes, optional
        Existing axes to draw on.
    show : bool, default=False
        Whether to call plt.show().

    Returns
    -------
    fig, ax
    """
    matrix = grid.matrix
    if depth is not None:
        slices = [slice(None)] * 3
        slices[axis] = slice(0, depth)
        matrix = matrix[tuple(slices)]

    # Project solid density along chosen axis
    proj = np.mean(matrix.astype(np.float32), axis=axis)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    im = ax.imshow(proj.T, cmap=cmap, origin="lower")
    ax.set_title(f"X-Ray Radiograph (Axis {axis})")
    ax.axis("off")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Solid Density")

    if show:
        plt.show()
    return fig, ax


def plot_slice(
    grid: VoxelGrid,
    slice_idx: Optional[int] = None,
    axis: int = 2,
    cmap: str = "gray_r",
    ax: Optional[plt.Axes] = None,
    show: bool = False,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot a 2D cross-sectional slice through the porous volume.
    """
    Nx, Ny, Nz = grid.shape
    axis_lens = [Nx, Ny, Nz]

    if slice_idx is None:
        slice_idx = axis_lens[axis] // 2

    slices = [slice(None)] * 3
    slices[axis] = slice_idx
    slice_2d = grid.matrix[tuple(slices)]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    im = ax.imshow(slice_2d.T, cmap=cmap, origin="lower")
    ax.set_title(f"Slice at {['X', 'Y', 'Z'][axis]} = {slice_idx}")
    ax.axis("off")

    if show:
        plt.show()
    return fig, ax


def plot_orthogonal_slices(
    grid: VoxelGrid,
    point: Optional[Tuple[int, int, int]] = None,
    cmap: str = "gray_r",
    figsize: Tuple[int, int] = (12, 4),
    show: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """
    Plot orthogonal tri-planar cross-sections intersecting at point (ix, iy, iz).
    """
    Nx, Ny, Nz = grid.shape
    if point is None:
        ix, iy, iz = Nx // 2, Ny // 2, Nz // 2
    else:
        ix, iy, iz = point

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # XY slice (Z = iz)
    axes[0].imshow(grid.matrix[:, :, iz].T, cmap=cmap, origin="lower")
    axes[0].axvline(ix, color="red", linestyle="--", alpha=0.5)
    axes[0].axhline(iy, color="green", linestyle="--", alpha=0.5)
    axes[0].set_title(f"XY Plane (Z={iz})")
    axes[0].axis("off")

    # XZ slice (Y = iy)
    axes[1].imshow(grid.matrix[:, iy, :].T, cmap=cmap, origin="lower")
    axes[1].axvline(ix, color="red", linestyle="--", alpha=0.5)
    axes[1].axhline(iz, color="blue", linestyle="--", alpha=0.5)
    axes[1].set_title(f"XZ Plane (Y={iy})")
    axes[1].axis("off")

    # YZ slice (X = ix)
    axes[2].imshow(grid.matrix[ix, :, :].T, cmap=cmap, origin="lower")
    axes[2].axvline(iy, color="green", linestyle="--", alpha=0.5)
    axes[2].axhline(iz, color="blue", linestyle="--", alpha=0.5)
    axes[2].set_title(f"YZ Plane (X={ix})")
    axes[2].axis("off")

    plt.tight_layout()
    if show:
        plt.show()
    return fig, axes
