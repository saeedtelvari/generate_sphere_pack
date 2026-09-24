"""
genrandsp.py
============
Backward-compatible wrapper for poropack.
Maintains the original generate_spherepack class and method_I / method_II APIs,
powered by the high-performance PeriodicSpatialGrid and RSAGenerator under the hood.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple, Union
import numpy as np
import pandas as pd

from poropack.core.psd import DiscretePSD
from poropack.core.representations import GrainPack
from poropack.core.spatial import PeriodicSpatialGrid
from poropack.generators.rsa import RSAGenerator


class generate_spherepack:
    """
    Legacy wrapper for 3D randomly generated sphere packs.
    Preserves original method signatures while leveraging poropack's O(1) spatial indexing.
    """

    def __init__(
        self,
        dimens: Optional[Sequence[float]] = None,
        Rs: Optional[Sequence[float]] = None,
        min_throat: float = 0.0,
        overlapping: bool = False,
        scale_factor: float = 1.0,
        loop_limits: int = 20000,
        box_dim: Optional[Sequence[float]] = None,
        radii: Optional[Sequence[float]] = None,
        throat_limit: Optional[float] = None,
    ):
        if dimens is None and box_dim is not None:
            dimens = box_dim
        if Rs is None and radii is not None:
            Rs = radii
        if throat_limit is not None:
            min_throat = throat_limit

        if dimens is None:
            dimens = [100.0, 100.0, 100.0]
        if Rs is None:
            Rs = [10.0, 5.0]

        self.dimens = np.asarray(dimens, dtype=np.float64)
        self.dx = float(self.dimens[0])
        self.dy = float(self.dimens[1])
        self.dz = float(self.dimens[2])
        self.Rs = list(Rs)
        self.rmax = float(np.max(self.Rs))
        self.rmin = float(np.min(self.Rs))
        self.min_throat = float(min_throat)
        self.overlapping = bool(overlapping)
        self.scale_factor = float(scale_factor)
        self.loop_limits = int(loop_limits)

        self.loop_counter = 0
        self.sphere_counter = 0
        self.sphere_coords: List[Tuple[float, float, float, float]] = []
        self.volume = self.dx * self.dy * self.dz
        self.spheres = {key: 0 for key in self.Rs}
        self.Rs_limit = 2000

    def is_overlapping(self, sphere1: Sequence[float], sphere2: Sequence[float], min_throat: float) -> bool:
        """Legacy pairwise collision check."""
        allowed_dist = sphere1[3] + sphere2[3] + min_throat
        distance = np.sqrt(
            (sphere1[0] - sphere2[0]) ** 2
            + (sphere1[1] - sphere2[1]) ** 2
            + (sphere1[2] - sphere2[2]) ** 2
        )
        return distance < allowed_dist

    def is_overlapping_all(self, x: float, y: float, z: float, r: float) -> List[bool]:
        """Legacy comparison against all previously generated spheres."""
        if not self.sphere_coords:
            return []
        temp = np.tile([x, y, z, r], len(self.sphere_coords)).reshape([len(self.sphere_coords), 4])
        minths = np.tile(self.min_throat, len(self.sphere_coords)).reshape([len(self.sphere_coords), 1])
        return list(map(self.is_overlapping, self.sphere_coords, temp, minths))

    def calc_porosity(self, sphere_coords: Sequence[Sequence[float]], volume: float) -> float:
        """Calculate theoretical porosity from non-overlapping sphere volume sum."""
        if len(sphere_coords) == 0:
            return 1.0
        coords_arr = np.asarray(sphere_coords)
        all_radius = coords_arr[:, 3]
        return 1.0 - np.sum(4.0 / 3.0 * np.pi * (all_radius ** 3)) / volume

    def generate_rand_sphere(self) -> Tuple[float, float, float, float]:
        """Generate a random coordinate and discrete radius."""
        r = float(np.random.choice(self.Rs))
        x = float(np.round(np.random.uniform(self.rmax, self.dx - self.rmax), decimals=6))
        y = float(np.round(np.random.uniform(self.rmax, self.dy - self.rmax), decimals=6))
        z = float(np.round(np.random.uniform(self.rmax, self.dz - self.rmax), decimals=6))
        return x, y, z, r

    def update_loop_counter(self, reset: bool = False):
        self.loop_counter = 0 if reset else self.loop_counter + 1

    def update_sphere_counter(self, reset: bool = False):
        self.sphere_counter = 0 if reset else self.sphere_counter + 1

    def print_progress(self, timestep: int = 500):
        if self.loop_counter % timestep == 0 and self.loop_counter > 0:
            poro = self.calc_porosity(self.sphere_coords, self.volume)
            print(f"Loop #{self.loop_counter:0>5}, n_spheres={self.sphere_counter:0>5}, porosity={poro:.4f}")

    def reset_all(self):
        self.sphere_coords = []
        self.loop_counter = 0
        self.sphere_counter = 0
        self.spheres = {key: 0 for key in self.Rs}

    def method_I(self, porosity: Optional[float] = None, target_porosity: Optional[float] = None) -> np.ndarray:
        """
        Generate spheres until target porosity is reached.
        Uses accelerated PeriodicSpatialGrid for O(1) collision detection.
        """
        if porosity is None:
            porosity = target_porosity if target_porosity is not None else 0.5
        self.reset_all()

        if self.overlapping:
            # Overlapping random generation
            while self.calc_porosity(self.sphere_coords, self.volume) >= porosity and self.loop_counter < self.loop_limits:
                self.update_loop_counter()
                x, y, z, r = self.generate_rand_sphere()
                self.sphere_coords.append((x, y, z, r))
                self.sphere_counter += 1
                self.spheres[r] = self.spheres.get(r, 0) + 1
                self.print_progress(500)
        else:
            # High-performance spatial-hashed RSA
            psd = DiscretePSD(self.Rs)
            gen = RSAGenerator(
                box_size=self.dimens,
                psd=psd,
                min_throat=self.min_throat,
                periodic=(False, False, False),  # match legacy wall-bounded mode
            )
            pack = gen.generate(target_porosity=porosity, max_attempts=self.loop_limits)
            self.loop_counter = pack.attributes.get("attempts", len(pack))
            self.sphere_counter = len(pack)
            for (x, y, z), r in zip(pack.coords, pack.radii):
                self.sphere_coords.append((float(x), float(y), float(z), float(r)))
                self.spheres[r] = self.spheres.get(r, 0) + 1

        print(f"\nFinal state: n_spheres={self.sphere_counter}, porosity={self.calc_porosity(self.sphere_coords, self.volume):.4f}")
        print(f"Number of spheres of each size : {self.spheres}")
        return np.array(self.sphere_coords)

    def method_II(self, num_spheres: Optional[int] = None, n_spheres: Optional[int] = None) -> np.ndarray:
        """
        Generate a specific number of spheres.
        Uses accelerated PeriodicSpatialGrid for O(1) collision detection.
        """
        if num_spheres is None:
            num_spheres = n_spheres if n_spheres is not None else 100
        self.reset_all()

        if self.overlapping:
            while self.sphere_counter < num_spheres and self.loop_counter < self.loop_limits:
                self.update_loop_counter()
                x, y, z, r = self.generate_rand_sphere()
                self.sphere_coords.append((x, y, z, r))
                self.sphere_counter += 1
                self.spheres[r] = self.spheres.get(r, 0) + 1
                self.print_progress(500)
        else:
            psd = DiscretePSD(self.Rs)
            gen = RSAGenerator(
                box_size=self.dimens,
                psd=psd,
                min_throat=self.min_throat,
                periodic=(False, False, False),
            )
            pack = gen.generate_by_count(n_spheres=num_spheres, max_attempts=self.loop_limits)
            self.loop_counter = pack.attributes.get("attempts", len(pack))
            self.sphere_counter = len(pack)
            for (x, y, z), r in zip(pack.coords, pack.radii):
                self.sphere_coords.append((float(x), float(y), float(z), float(r)))
                self.spheres[r] = self.spheres.get(r, 0) + 1

        print(f"\nFinal state: n_spheres={self.sphere_counter}, porosity={self.calc_porosity(self.sphere_coords, self.volume):.4f}")
        print(f"Number of spheres of each size : {self.spheres}")
        return np.array(self.sphere_coords)
