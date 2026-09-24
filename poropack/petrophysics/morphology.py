"""
poropack.petrophysics.morphology
================================
Morphological characterization of porous media:
- Specific Surface Area (Sv) via Marching Cubes isosurfaces
- Two-Point Spatial Autocorrelation Function S2(r) via 3D FFT
- Directional Chord Length Distributions
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple
import numpy as np
from scipy import ndimage

from poropack.core.representations import VoxelGrid


def specific_surface_area(
    grid: VoxelGrid,
    method: str = "marching_cubes",
) -> float:
    """
    Calculate the specific surface area (interfacial area per unit bulk volume: Sv = A / V).

    Parameters
    ----------
    grid : VoxelGrid
        Porous media grid.
    method : {'marching_cubes', 'voxel_faces'}, default='marching_cubes'
        Surface measurement algorithm.

    Returns
    -------
    float
        Specific surface area Sv in units of 1 / [voxel_size] (e.g. 1/um or 1/m).
    """
    try:
        from skimage import measure
        has_skimage = True
    except ImportError:
        has_skimage = False

    h = grid.voxel_size
    V_total = grid.volume

    if method == "marching_cubes" and has_skimage:
        # Marching cubes on solid boolean field
        matrix_float = grid.matrix.astype(np.float32)
        # Pad with zeros at borders if not strictly periodic to close boundary triangles
        verts, faces, normals, values = measure.marching_cubes(
            matrix_float, level=0.5, spacing=(h, h, h)
        )
        surface_area = float(measure.mesh_surface_area(verts, faces))
        return surface_area / V_total

    else:
        # Fallback: Count unshared voxel face boundaries between pore and solid
        # Face area is h^2
        diff_x = np.diff(grid.matrix.astype(np.int8), axis=0) != 0
        diff_y = np.diff(grid.matrix.astype(np.int8), axis=1) != 0
        diff_z = np.diff(grid.matrix.astype(np.int8), axis=2) != 0
        total_faces = (
            np.count_nonzero(diff_x)
            + np.count_nonzero(diff_y)
            + np.count_nonzero(diff_z)
        )
        surface_area = total_faces * (h ** 2)
        return float(surface_area / V_total)


def two_point_correlation(
    grid: VoxelGrid,
    max_radius: Optional[float] = None,
    n_bins: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the 3D two-point spatial autocorrelation function S2(r) of the pore space via 3D FFT.
    S2(0) = porosity, and S2(r -> inf) = porosity^2.

    Parameters
    ----------
    grid : VoxelGrid
        Input porous volume.
    max_radius : float, optional
        Maximum physical correlation lag distance r to evaluate.
        Defaults to min(Lx, Ly, Lz) / 3.
    n_bins : int, default=50
        Number of radial distance bins.

    Returns
    -------
    r_centers : np.ndarray
        Radial separation distances.
    S2 : np.ndarray
        Autocorrelation probability values.
    """
    pore_field = (~grid.matrix).astype(np.float64)
    Nx, Ny, Nz = grid.shape
    h = grid.voxel_size

    if max_radius is None:
        max_radius = min(Nx, Ny, Nz) * h / 3.0

    # 1. 3D Autocorrelation via Wiener-Khinchin theorem
    fft_f = np.fft.fftn(pore_field)
    power_spec = np.abs(fft_f) ** 2
    autocorr = np.real(np.fft.ifftn(power_spec)) / pore_field.size

    # Shift zero-lag to center
    autocorr_shifted = np.fft.fftshift(autocorr)

    # 2. Compute radial coordinates from center
    cx, cy, cz = Nx // 2, Ny // 2, Nz // 2
    x = (np.arange(Nx) - cx) * h
    y = (np.arange(Ny) - cy) * h
    z = (np.arange(Nz) - cz) * h

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    # 3. Azimuthal radial binning
    r_flat = R.ravel()
    corr_flat = autocorr_shifted.ravel()

    valid_mask = r_flat <= max_radius
    r_valid = r_flat[valid_mask]
    corr_valid = corr_flat[valid_mask]

    bins = np.linspace(0.0, max_radius, n_bins + 1)
    bin_idx = np.digitize(r_valid, bins) - 1

    r_centers = 0.5 * (bins[:-1] + bins[1:])
    s2_vals = np.zeros(n_bins, dtype=np.float64)

    for b in range(n_bins):
        mask_b = bin_idx == b
        if np.any(mask_b):
            s2_vals[b] = np.mean(corr_valid[mask_b])
        else:
            s2_vals[b] = np.nan

    return r_centers, s2_vals


def chord_length_distribution(
    grid: VoxelGrid,
    axis: int = 0,
    phase: str = "pore",
) -> np.ndarray:
    """
    Extract 1D chord length distribution along a given axis.

    Parameters
    ----------
    grid : VoxelGrid
        Input grid.
    axis : {0, 1, 2}
        Direction of 1D line probes.
    phase : {'pore', 'solid'}, default='pore'
        Phase whose chord lengths are measured.

    Returns
    -------
    np.ndarray
        Array of physical chord lengths.
    """
    mask = (~grid.matrix) if phase == "pore" else grid.matrix
    h = grid.voxel_size

    # Move target axis to the last dimension
    transposed = np.moveaxis(mask, axis, -1)
    chords: list[float] = []

    # Iterate over 1D rays
    flat_rays = transposed.reshape(-1, transposed.shape[-1])
    for ray in flat_rays:
        # Run-length encoding of True segments
        padded = np.pad(ray, (1, 1), mode="constant", constant_values=False)
        diff = np.diff(padded.astype(np.int8))
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        lengths = (ends - starts) * h
        chords.extend(lengths)

    return np.array(chords, dtype=np.float64)
