"""
poropack.diagenesis.compaction
==============================
Post-depositional mechanical compaction modeling.
Applies vertical strain tensors to Lagrangian grain assemblies, inducing realistic grain interpenetration.
"""

from __future__ import annotations

from typing import Optional, Union
import numpy as np
from scipy import ndimage

from poropack.core.representations import GrainPack, VoxelGrid


def apply_compaction(
    target: Union[GrainPack, VoxelGrid],
    vertical_strain: Optional[float] = None,
    poisson_ratio: float = 0.0,
    strain_z: Optional[float] = None,
) -> Union[GrainPack, VoxelGrid]:
    """
    Simulate post-depositional uniaxial or semi-confined mechanical compaction.
    Supports both Lagrangian GrainPack and Eulerian VoxelGrid representations.

    Parameters
    ----------
    target : GrainPack or VoxelGrid
        Original uncompacted grain assembly or binary voxel grid.
    vertical_strain : float, optional
        Uniaxial compressive strain along z-axis (e.g. 0.05 to 0.20 for 5% to 20% vertical compression).
    poisson_ratio : float, default=0.0
        Poisson's ratio nu for lateral expansion (nu = 0 for standard uniaxial 1D strain).
    strain_z : float, optional
        Alias for vertical_strain.

    Returns
    -------
    GrainPack or VoxelGrid
        Compacted assembly or voxel grid with squashed vertical dimension.
    """
    if vertical_strain is None:
        if strain_z is not None:
            vertical_strain = strain_z
        else:
            raise ValueError("Must specify vertical_strain (or strain_z)")

    if not (0.0 <= vertical_strain < 1.0):
        raise ValueError(f"vertical_strain must be in [0, 1), got {vertical_strain}")

    ezz = float(vertical_strain)
    nu = float(poisson_ratio)

    exx = nu * ezz
    eyy = nu * ezz

    if isinstance(target, GrainPack):
        new_coords = target.coords.copy()
        new_coords[:, 0] *= (1.0 + exx)
        new_coords[:, 1] *= (1.0 + eyy)
        new_coords[:, 2] *= (1.0 - ezz)

        new_box = target.box_size.copy()
        new_box[0] *= (1.0 + exx)
        new_box[1] *= (1.0 + eyy)
        new_box[2] *= (1.0 - ezz)

        attrs = {k: v.copy() if hasattr(v, "copy") else v for k, v in target.attributes.items()}
        attrs["vertical_strain"] = ezz
        attrs["poisson_ratio"] = nu

        return GrainPack(
            coords=new_coords,
            radii=target.radii.copy(),
            box_size=new_box,
            periodic=target.periodic,
            attributes=attrs,
        )

    elif isinstance(target, VoxelGrid):
        zoom_factors = (1.0 + exx, 1.0 + eyy, 1.0 - ezz)
        compacted_matrix = ndimage.zoom(
            target.matrix.astype(np.float32),
            zoom=zoom_factors,
            order=0,
            mode="nearest",
        ) >= 0.5

        return VoxelGrid(
            matrix=compacted_matrix,
            voxel_size=target.voxel_size,
            origin=target.origin.copy(),
            periodic=target.periodic,
        )
    else:
        raise TypeError(f"Expected GrainPack or VoxelGrid, got {type(target)}")
