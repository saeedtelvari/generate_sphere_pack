"""
poropack.diagenesis.cementation
===============================
Morphological cementation (quartz syntaxial overgrowth, pore-filling calcite)
and secondary dissolution / vug creation using Euclidean Distance Transforms (EDT).
"""

from __future__ import annotations

from typing import Optional, Tuple
import numpy as np
from scipy import ndimage

from poropack.core.representations import VoxelGrid


def apply_cementation(
    grid: VoxelGrid,
    cement_fraction: Optional[float] = None,
    distance_threshold: Optional[float] = None,
) -> VoxelGrid:
    """
    Simulate diagenetic cementation via Euclidean Distance Transform (EDT) dilation.
    Fills narrow pore throats preferentially before larger pore bodies.

    Parameters
    ----------
    grid : VoxelGrid
        Input voxel grid representing the porous medium before cementation.
    cement_fraction : float, optional
        Fraction of bulk volume to deposit as secondary mineral cement (e.g. 0.08 for 8% cement).
    distance_threshold : float, optional
        Explicit physical distance (in physical units) from grain surfaces to fill with cement.

    Returns
    -------
    VoxelGrid
        Cemented porous media with reduced porosity and occluded pore throats.
    """
    if cement_fraction is None and distance_threshold is None:
        raise ValueError("Must specify either cement_fraction or distance_threshold")

    # True = solid, False = pore. We compute distance from solid boundaries into pore space:
    pore_mask = ~grid.matrix

    if not np.any(pore_mask):
        return grid.copy()

    # Compute Euclidean distance transform in physical units
    # sampling=voxel_size scales distances accurately
    dt = ndimage.distance_transform_edt(pore_mask, sampling=grid.voxel_size)

    new_matrix = grid.matrix.copy()

    if distance_threshold is not None:
        cement_mask = pore_mask & (dt <= distance_threshold)
        new_matrix[cement_mask] = True
    else:
        # Determine threshold to match target cement volume fraction
        target_cement_voxels = int(np.round(float(cement_fraction) * grid.matrix.size))
        pore_dt_values = dt[pore_mask]
        sorted_dt = np.sort(pore_dt_values)

        if target_cement_voxels >= len(sorted_dt):
            # Cements entire pore space
            new_matrix[:] = True
        elif target_cement_voxels > 0:
            crit_dt = sorted_dt[target_cement_voxels - 1]
            cement_mask = pore_mask & (dt <= crit_dt)
            new_matrix[cement_mask] = True

    return VoxelGrid(
        matrix=new_matrix,
        voxel_size=grid.voxel_size,
        origin=grid.origin.copy(),
        periodic=grid.periodic,
    )


def apply_dissolution(
    grid: VoxelGrid,
    porosity_increase: Optional[float] = None,
    mode: str = "pore_surface",
    random_state: Optional[int | np.random.Generator] = None,
    dissolution_fraction: Optional[float] = None,
) -> VoxelGrid:
    """
    Simulate secondary dissolution / leaching to create secondary micro-porosity and vugs.

    Parameters
    ----------
    grid : VoxelGrid
        Porous rock grid.
    porosity_increase : float, optional
        Additional porosity to generate (e.g. 0.05 for +5% secondary porosity).
    mode : {'pore_surface', 'stochastic_vugs', 'grain_boundary'}, default='pore_surface'
        - 'pore_surface' / 'grain_boundary': Expands pore surfaces via solid erosion (acid leaching).
        - 'stochastic_vugs': Creates disconnected secondary vugs/cavities.
    random_state : int or Generator, optional
        Seed or random generator.
    dissolution_fraction : float, optional
        Alias for porosity_increase.

    Returns
    -------
    VoxelGrid
    """
    if porosity_increase is None:
        if dissolution_fraction is not None:
            porosity_increase = dissolution_fraction
        else:
            raise ValueError("Must specify porosity_increase (or dissolution_fraction)")

    if porosity_increase <= 0.0:
        return grid.copy()

    target_extra_voxels = int(np.round(porosity_increase * grid.matrix.size))
    new_matrix = grid.matrix.copy()
    solid_mask = grid.solid_mask

    if not np.any(solid_mask):
        return grid.copy()

    if mode in ("pore_surface", "grain_boundary", "matrix"):
        # Distance from pore space into solid grains (surface leaching)
        dt_solid = ndimage.distance_transform_edt(solid_mask, sampling=grid.voxel_size)
        solid_dt_values = dt_solid[solid_mask]
        sorted_dt = np.sort(solid_dt_values)

        if target_extra_voxels >= len(sorted_dt):
            new_matrix[:] = False
        else:
            crit_dt = sorted_dt[target_extra_voxels - 1]
            dissolved = solid_mask & (dt_solid <= crit_dt)
            new_matrix[dissolved] = False

    elif mode in ("grain_cores", "cores"):
        # Intra-granular leaching of grain cores (where dt_solid is largest)
        dt_solid = ndimage.distance_transform_edt(solid_mask, sampling=grid.voxel_size)
        solid_dt_values = dt_solid[solid_mask]
        sorted_dt = np.sort(solid_dt_values)

        if target_extra_voxels >= len(sorted_dt):
            new_matrix[:] = False
        else:
            idx = max(0, len(sorted_dt) - target_extra_voxels)
            crit_dt = sorted_dt[idx]
            dissolved = solid_mask & (dt_solid >= crit_dt)
            new_matrix[dissolved] = False

    elif mode in ("stochastic_vugs", "vugs"):
        rng = np.random.default_rng(random_state)
        # Random spherical seeds inside solid phase
        solid_indices = np.argwhere(solid_mask)
        if len(solid_indices) > 0:
            n_vugs = max(1, target_extra_voxels // 50)
            choice_idx = rng.choice(len(solid_indices), size=n_vugs, replace=False)
            vug_centers = solid_indices[choice_idx]

            # Mark seeds and dilate
            seed_mask = np.zeros_like(grid.matrix, dtype=bool)
            seed_mask[vug_centers[:, 0], vug_centers[:, 1], vug_centers[:, 2]] = True

            vug_dt = ndimage.distance_transform_edt(~seed_mask, sampling=grid.voxel_size)
            sorted_vug_dt = np.sort(vug_dt[solid_mask])
            k = min(len(sorted_vug_dt) - 1, target_extra_voxels)
            crit_vug_dt = sorted_vug_dt[k]

            dissolved = solid_mask & (vug_dt <= crit_vug_dt)
            new_matrix[dissolved] = False
    else:
        raise ValueError(f"Unknown dissolution mode: {mode}")

    return VoxelGrid(
        matrix=new_matrix,
        voxel_size=grid.voxel_size,
        origin=grid.origin.copy(),
        periodic=grid.periodic,
    )


def dissolve(
    grid: VoxelGrid,
    porosity_increase: Optional[float] = None,
    mode: str = "matrix",
    random_state: Optional[int | np.random.Generator] = None,
    dissolution_fraction: Optional[float] = None,
) -> VoxelGrid:
    """Alias for apply_dissolution matching Ticket 07 specification."""
    return apply_dissolution(
        grid=grid,
        porosity_increase=porosity_increase,
        mode=mode,
        random_state=random_state,
    )
