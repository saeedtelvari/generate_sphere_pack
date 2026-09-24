"""
poropack.generators
===================
Porous media generation engines:
- RSAGenerator: Accelerated Spatial-Hashed Random Sequential Adsorption
- SedimentationGenerator: 3D Ballistic Gravitational Sedimentation (Drop-and-Roll)
- DensePackingGenerator: Force-Biased Algorithm (FBA) with FIRE optimization
- GaussianRandomField: 3D Level-Cut Spectral Fourier Random Fields
- QSGSGenerator: Quartet Structure Generation Set
"""

from poropack.generators.rsa import RSAGenerator
from poropack.generators.sedimentation import SedimentationGenerator
from poropack.generators.dense_packing import DensePackingGenerator
from poropack.generators.continuum import GaussianRandomField, QSGSGenerator

__all__ = [
    "RSAGenerator",
    "SedimentationGenerator",
    "DensePackingGenerator",
    "GaussianRandomField",
    "QSGSGenerator",
]
