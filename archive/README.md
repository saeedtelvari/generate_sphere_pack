# Legacy Code & Experiments Archive

This directory preserves early exploratory scripts and notebooks from the original repository (2020-2022) for historical reference and provenance.

### Contents:
- **`examples_legacy.ipynb`**: Original demonstration notebook showing the brute-force rejection algorithm (`method_I` and `method_II`).
- **`functions_legacy.py`**: Initial procedural implementation before OOP refactoring.
- **`history_legacy.ipynb`** & **`main_functional_legacy.ipynb`**: Exploratory prototyping notebooks.
- **`ordered_sphere_pack_legacy.ipynb`**: Preliminary experiments with regular sphere packing geometries.
- **`images/`**: Historical rendering artifacts.

### Modern Replacement:
All capabilities in this archive have been superseded by the production-ready **`poropack`** package:
- $O(1)$ spatial hash grid collision detection (`poropack.core.spatial.PeriodicSpatialGrid`) replacing $O(N^2)$ brute-force loops.
- 4 physics-grounded generation engines: RSA, Ballistic Sedimentation ("Drop-and-Roll"), Dense Random Packing (FBA + FIRE), and Continuum Spectral GRF / QSGS.
- Continuous geological PSDs (LogNormal with Folk-Ward sorting, Weibull, Sieve curves with Hatch-Choate volume-to-number conversion).
- Geological diagenesis pipeline (compaction, cementation, dissolution).
- Comprehensive petrophysics suite (porosity partitioning, $S_v$, $S_2(r)$, chords, Kozeny-Carman, SNOW2 PNM).
- Multiformat exporters (3D TIFF, ParaView VTK, MATLAB MRST, Parquet, CSV).
- For modern interactive tutorials, see [`notebooks/`](../notebooks/) and [`examples_v2.ipynb`](../examples_v2.ipynb).
- Backward compatibility is maintained at the root via [`genrandsp.py`](../genrandsp.py).
