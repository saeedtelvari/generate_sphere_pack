"""
poropack.petrophysics
=====================
Integrated digital rock physics and petrophysics suite:
- analyze_porosity: Total, effective, and isolated void space breakdown
- specific_surface_area: Interfacial area per volume via Marching Cubes
- two_point_correlation: Spatial autocorrelation S2(r) via 3D FFT
- chord_length_distribution: 1D directional chord lengths
- kozeny_carman: Analytical permeability estimation
- formation_factor: Archie electrical formation factor
- extract_pore_network: SNOW2 topological pore-throat extraction
"""

from poropack.petrophysics.porosity import analyze_porosity, PorosityResult
from poropack.petrophysics.morphology import (
    specific_surface_area,
    two_point_correlation,
    chord_length_distribution,
)
from poropack.petrophysics.transport import kozeny_carman, formation_factor
from poropack.petrophysics.pnm import extract_pore_network

__all__ = [
    "analyze_porosity",
    "PorosityResult",
    "specific_surface_area",
    "two_point_correlation",
    "chord_length_distribution",
    "kozeny_carman",
    "formation_factor",
    "extract_pore_network",
]
