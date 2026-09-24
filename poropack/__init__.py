"""
poropack
========
A high-performance Python library for generating, diagenetically modifying,
and petrophysically characterizing 3D multi-purpose porous media for Digital Rock Physics.
"""

__version__ = "2.0.0"

# Core representations & distributions
from poropack.core.representations import GrainPack, VoxelGrid, rasterize_spheres
from poropack.core.spatial import PeriodicSpatialGrid
from poropack.core.psd import (
    BasePSD,
    LogNormalPSD,
    WeibullPSD,
    TruncatedGaussianPSD,
    DiscretePSD,
    EmpiricalSievePSD,
)

# Generators
from poropack.generators.rsa import RSAGenerator
from poropack.generators.sedimentation import SedimentationGenerator
from poropack.generators.dense_packing import DensePackingGenerator
from poropack.generators.continuum import GaussianRandomField, QSGSGenerator

# Diagenesis
from poropack.diagenesis.compaction import apply_compaction
from poropack.diagenesis.cementation import apply_cementation, apply_dissolution, dissolve

# Petrophysics
from poropack.petrophysics.porosity import analyze_porosity, PorosityResult
from poropack.petrophysics.morphology import (
    specific_surface_area,
    two_point_correlation,
    chord_length_distribution,
)
from poropack.petrophysics.transport import kozeny_carman, formation_factor
from poropack.petrophysics.pnm import extract_pore_network

# I/O & Visualization
from poropack.io.exporters import (
    export_csv,
    export_parquet,
    export_tiff,
    export_raw,
    export_vtk_vti,
    export_vtk_vtp,
    export_mrst_mat,
)
from poropack.io.visualization import (
    plot_xray,
    plot_slice,
    plot_orthogonal_slices,
)

__all__ = [
    "__version__",
    # Core
    "GrainPack",
    "VoxelGrid",
    "rasterize_spheres",
    "PeriodicSpatialGrid",
    "BasePSD",
    "LogNormalPSD",
    "WeibullPSD",
    "TruncatedGaussianPSD",
    "DiscretePSD",
    "EmpiricalSievePSD",
    # Generators
    "RSAGenerator",
    "SedimentationGenerator",
    "DensePackingGenerator",
    "GaussianRandomField",
    "QSGSGenerator",
    # Diagenesis
    "apply_compaction",
    "apply_cementation",
    "apply_dissolution",
    "dissolve",
    # Petrophysics
    "analyze_porosity",
    "PorosityResult",
    "specific_surface_area",
    "two_point_correlation",
    "chord_length_distribution",
    "kozeny_carman",
    "formation_factor",
    "extract_pore_network",
    # I/O
    "export_csv",
    "export_parquet",
    "export_tiff",
    "export_raw",
    "export_vtk_vti",
    "export_vtk_vtp",
    "export_mrst_mat",
    "plot_xray",
    "plot_slice",
    "plot_orthogonal_slices",
]
