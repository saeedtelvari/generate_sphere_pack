"""
poropack.petrophysics.transport
===============================
Transport property estimators:
- Kozeny-Carman analytical permeability
- Archie formation factor and tortuosity
"""

from __future__ import annotations

from typing import Optional
import numpy as np

from poropack.core.representations import VoxelGrid
from poropack.petrophysics.morphology import specific_surface_area
from poropack.petrophysics.porosity import analyze_porosity


def kozeny_carman(
    grid: VoxelGrid,
    kozeny_constant: float = 5.0,
    unit: str = "physical",
) -> float:
    """
    Estimate absolute permeability using the Kozeny-Carman equation:
    k = phi_eff^3 / (c * Sv^2)

    Parameters
    ----------
    grid : VoxelGrid
        Porous media volume.
    kozeny_constant : float, default=5.0
        Kozeny geometric tortuosity constant (typically ~5.0 for sphere packs).
    unit : {'physical', 'darcy', 'millidarcy'}, default='physical'
        Output units. If 'physical', returns length^2 (e.g. um^2).
        If 'darcy', assumes voxel_size is in microns (1 Darcy = 0.986923 um^2).

    Returns
    -------
    float
        Estimated absolute permeability.
    """
    poro_res = analyze_porosity(grid)
    phi_eff = poro_res.effective if poro_res.effective > 0 else poro_res.total

    if phi_eff <= 0.0:
        return 0.0

    sv = specific_surface_area(grid)
    if sv <= 0.0:
        return 0.0

    # k in [voxel_size]^2
    k_raw = (phi_eff ** 3) / (kozeny_constant * (sv ** 2))

    if unit.lower() == "darcy":
        # 1 um^2 = 1.01325 Darcy
        return k_raw / 0.986923
    elif unit.lower() in ("millidarcy", "md"):
        return (k_raw / 0.986923) * 1000.0
    else:
        return k_raw


def formation_factor(
    grid: VoxelGrid,
    cementation_exponent: float = 2.0,
    tortuosity_factor: float = 1.0,
) -> float:
    """
    Estimate electrical Formation Resistivity Factor (F) via Archie's first law:
    F = a * phi^(-m)

    Parameters
    ----------
    grid : VoxelGrid
        Porous rock grid.
    cementation_exponent : float, default=2.0
        Archie's cementation exponent m (typically 1.8 - 2.2 for sandstones).
    tortuosity_factor : float, default=1.0
        Structural coefficient a.

    Returns
    -------
    float
        Electrical formation factor F.
    """
    poro = analyze_porosity(grid).effective
    if poro <= 0.0:
        return float("inf")
    return float(tortuosity_factor * (poro ** (-cementation_exponent)))
