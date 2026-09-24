"""
poropack.generators.rsa
=======================
High-performance Random Sequential Adsorption (RSA) hard-sphere generator
accelerated by the PeriodicSpatialGrid (O(N) total packing time).
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence, Tuple, Union
import numpy as np

from poropack.core.psd import BasePSD, DiscretePSD
from poropack.core.representations import GrainPack
from poropack.core.spatial import PeriodicSpatialGrid


class RSAGenerator:
    """
    Accelerated Random Sequential Adsorption (RSA) hard-sphere generator.

    Parameters
    ----------
    box_size : tuple or ndarray of float
        (Lx, Ly, Lz) dimensions of the packing domain.
    psd : BasePSD or sequence of float
        Particle size distribution or list of discrete radii.
    min_throat : float, default=0.0
        Minimum clearance / throat distance between sphere surfaces.
    periodic : tuple of bool, default=(True, True, True)
        Whether to wrap coordinates periodically along (x, y, z).
    random_state : int or Generator, optional
        Seed or NumPy random Generator for reproducible generation.
    """

    def __init__(
        self,
        box_size: Tuple[float, float, float] | np.ndarray,
        psd: Union[BasePSD, Sequence[float]],
        min_throat: float = 0.0,
        periodic: Tuple[bool, bool, bool] = (True, True, True),
        random_state: Optional[int | np.random.Generator] = None,
        throat_cutoff: Optional[float] = None,
    ):
        self.box_size = np.asarray(box_size, dtype=np.float64)
        if isinstance(psd, BasePSD):
            self.psd = psd
        else:
            self.psd = DiscretePSD(psd)

        if throat_cutoff is not None:
            min_throat = throat_cutoff
        self.min_throat = float(min_throat)
        self.periodic = tuple(bool(p) for p in periodic)

        if isinstance(random_state, np.random.Generator):
            self.rng = random_state
        else:
            self.rng = np.random.default_rng(random_state)

        r_min, r_max = self.psd.get_bounds()
        self.r_max = r_max
        self.r_min = r_min
        self.volume = float(np.prod(self.box_size))

    def _sample_position(self, r: float) -> Tuple[float, float, float]:
        """Generate a random coordinate respecting boundary conditions."""
        pos = np.zeros(3, dtype=np.float64)
        for d in range(3):
            if self.periodic[d]:
                pos[d] = self.rng.uniform(0.0, self.box_size[d])
            else:
                pos[d] = self.rng.uniform(r, self.box_size[d] - r)
        return float(pos[0]), float(pos[1]), float(pos[2])

    def generate(
        self,
        target_porosity: float,
        max_attempts: int = 200000,
        progress_callback: Optional[Callable[[int, int, float], None]] = None,
    ) -> GrainPack:
        """
        Generate spheres until the target porosity criterion is met or jamming limit is reached.

        Parameters
        ----------
        target_porosity : float
            Desired porosity (e.g. 0.65). Note that monodisperse hard spheres jam at ~0.618.
        max_attempts : int, default=200000
            Maximum consecutive candidate trial insertions before halting.
        progress_callback : callable, optional
            Function called with (attempt_count, sphere_count, current_porosity).

        Returns
        -------
        GrainPack
            Generated hard-sphere assembly.
        """
        grid = PeriodicSpatialGrid(
            box_size=self.box_size,
            r_max=self.r_max,
            min_throat=self.min_throat,
            periodic=self.periodic,
        )

        current_solid_vol = 0.0
        current_poro = 1.0
        attempts = 0

        while current_poro > target_porosity and attempts < max_attempts:
            attempts += 1
            r = float(self.psd.sample(1, rng=self.rng)[0])
            x, y, z = self._sample_position(r)

            if not grid.check_collision(x, y, z, r):
                grid.insert(x, y, z, r)
                current_solid_vol += 4.0 / 3.0 * np.pi * (r ** 3)
                current_poro = 1.0 - (current_solid_vol / self.volume)

            if progress_callback and attempts % 1000 == 0:
                progress_callback(attempts, len(grid), current_poro)

        coords, radii = grid.get_arrays()
        return GrainPack(
            coords=coords,
            radii=radii,
            box_size=self.box_size,
            periodic=self.periodic,
            attributes={"attempts": attempts},
        )

    def generate_by_count(
        self,
        n_spheres: int,
        max_attempts: int = 500000,
        progress_callback: Optional[Callable[[int, int, float], None]] = None,
    ) -> GrainPack:
        """
        Generate a specific number of non-overlapping spheres.

        Parameters
        ----------
        n_spheres : int
            Total number of spheres to place.
        max_attempts : int, default=500000
            Maximum total candidate trials.
        progress_callback : callable, optional
            Progress monitoring callback.

        Returns
        -------
        GrainPack
        """
        grid = PeriodicSpatialGrid(
            box_size=self.box_size,
            r_max=self.r_max,
            min_throat=self.min_throat,
            periodic=self.periodic,
        )

        attempts = 0
        current_solid_vol = 0.0

        while len(grid) < n_spheres and attempts < max_attempts:
            attempts += 1
            r = float(self.psd.sample(1, rng=self.rng)[0])
            x, y, z = self._sample_position(r)

            if not grid.check_collision(x, y, z, r):
                grid.insert(x, y, z, r)
                current_solid_vol += 4.0 / 3.0 * np.pi * (r ** 3)

            if progress_callback and attempts % 1000 == 0:
                poro = 1.0 - (current_solid_vol / self.volume)
                progress_callback(attempts, len(grid), poro)

        coords, radii = grid.get_arrays()
        return GrainPack(
            coords=coords,
            radii=radii,
            box_size=self.box_size,
            periodic=self.periodic,
            attributes={"attempts": attempts},
        )
