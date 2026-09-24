"""
poropack.petrophysics.pnm
=========================
Pore Network Model (PNM) extraction bridge to PoreSpy and OpenPNM.
"""

from __future__ import annotations

from typing import Any, Dict, Optional
import numpy as np

from poropack.core.representations import VoxelGrid


def extract_pore_network(
    grid: VoxelGrid,
    voxel_size: Optional[float] = None,
    accuracy: str = "high",
) -> Dict[str, Any]:
    """
    Extract a topological pore network from VoxelGrid using PoreSpy's SNOW2 algorithm.

    Parameters
    ----------
    grid : VoxelGrid
        Porous media grid.
    voxel_size : float, optional
        Voxel resolution. Defaults to grid.voxel_size.
    accuracy : {'high', 'standard'}, default='high'
        SNOW algorithm accuracy mode.

    Returns
    -------
    dict
        OpenPNM-compatible dictionary containing pore and throat geometries:
        'pore.coords', 'pore.diameter', 'throat.conns', 'throat.diameter', etc.
    """
    try:
        import porespy as ps
    except ImportError as e:
        raise ImportError(
            "Pore network extraction requires 'porespy'. Install it via 'pip install porespy'."
        ) from e

    h = grid.voxel_size if voxel_size is None else float(voxel_size)

    # PoreSpy 3.x convention: phases array where pore is 1 and solid is 0
    phases = (~grid.matrix).astype(int)

    # Execute SNOW2 network extraction
    snow_out = ps.networks.snow2(
        phases=phases,
        voxel_size=h,
        accuracy=accuracy,
    )

    return snow_out.network
