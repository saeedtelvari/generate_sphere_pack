"""
poropack.io
===========
Export and visualization interfaces:
- export_tiff: Multi-page TIFF stacks
- export_raw: Raw binary format with MHD header
- export_vtk_vti: Structured image VTK XML for ParaView
- export_vtk_vtp: Particle and contact network PolyData VTK XML
- export_mrst_mat: MATLAB MRST grid and rock struct
- plot_xray, plot_slice, plot_orthogonal_slices: Matplotlib visualizers
"""

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
