"""
poropack.generators.dense_packing
=================================
Dense Random Packing (DRP) of hard spheres breaking the RSA jamming limit.
Uses the Force-Biased Algorithm (FBA) with Hookean potential and FIRE optimization.
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence, Tuple, Union
import numpy as np
import scipy.spatial as sp

from poropack.core.psd import BasePSD, DiscretePSD
from poropack.core.representations import GrainPack


class DensePackingGenerator:
    """
    Collective rearrangement engine generating Dense Random Packings (DRP)
    of hard spheres (solid fraction >= 0.640, porosity <= 0.360) breaking the RSA limit.

    Parameters
    ----------
    box_size : tuple or ndarray of float
        (Lx, Ly, Lz) dimensions of the packing domain.
    psd : BasePSD or sequence of float, optional
        Particle size distribution. If None, generates monodisperse packing.
    spring_constant : float, default=1.0
        Hookean repulsive spring stiffness.
    random_state : int or Generator, optional
        Seed or random generator.
    """

    def __init__(
        self,
        box_size: Tuple[float, float, float] | np.ndarray = (100.0, 100.0, 100.0),
        psd: Optional[Union[BasePSD, Sequence[float]]] = None,
        spring_constant: float = 1.0,
        random_state: Optional[int | np.random.Generator] = None,
    ):
        self.box_size = np.asarray(box_size, dtype=np.float64)
        self.volume = float(np.prod(self.box_size))
        self.k = float(spring_constant)

        if psd is None:
            self.psd = None
        elif isinstance(psd, BasePSD):
            self.psd = psd
        else:
            self.psd = DiscretePSD(psd)

        if isinstance(random_state, np.random.Generator):
            self.rng = random_state
        else:
            self.rng = np.random.default_rng(random_state)

    def _compute_contacts_and_forces(
        self,
        coords: np.ndarray,
        radii: np.ndarray,
        r_max: float,
    ) -> Tuple[np.ndarray, float, float, int]:
        """
        Vectorized evaluation of periodic minimum-image Hookean repulsive forces.
        """
        tree = sp.cKDTree(coords, boxsize=self.box_size)
        pairs = tree.query_pairs(2.0 * r_max, output_type="ndarray")

        if len(pairs) == 0:
            return np.zeros_like(coords), 0.0, 0.0, 0

        i, j = pairs[:, 0], pairs[:, 1]

        # Minimum Image Convention displacement
        dr = coords[j] - coords[i]
        dr -= self.box_size * np.round(dr / self.box_size)
        dist_sq = np.sum(dr ** 2, axis=1)
        dist = np.sqrt(dist_sq)

        # Overlaps
        sigma = radii[i] + radii[j]
        overlap = sigma - dist
        active = (overlap > 1e-12) & (dist > 1e-12)

        if not np.any(active):
            return np.zeros_like(coords), 0.0, 0.0, 0

        i_act = i[active]
        j_act = j[active]
        dr_act = dr[active]
        dist_act = dist[active, None]
        ov_act = overlap[active, None]

        # Hookean contact force: F_ij = k * delta * (dr / dist)
        f_ij = self.k * ov_act * (dr_act / dist_act)

        # Vectorized scatter-add
        forces = np.zeros_like(coords)
        np.add.at(forces, i_act, -f_ij)
        np.add.at(forces, j_act, f_ij)

        energy = float(0.5 * self.k * np.sum(ov_act ** 2))
        max_overlap = float(np.max(ov_act))
        contact_count = int(np.count_nonzero(active))

        return forces, energy, max_overlap, contact_count

    def _relax_fire(
        self,
        coords: np.ndarray,
        radii: np.ndarray,
        max_iters: int = 600,
        tol_overlap_rel: float = 5e-5,
    ) -> Tuple[np.ndarray, float, int]:
        """Fast Inertial Relaxation Engine (FIRE) energy minimization."""
        r_mean = float(np.mean(radii))
        r_max = float(np.max(radii))
        v = np.zeros_like(coords)

        # FIRE hyper-parameters
        dt = 0.02
        dt_max = 0.10
        dt_min = 1e-4
        alpha = 0.10
        n_pos = 0
        N_min = 5
        f_inc = 1.10
        f_dec = 0.50
        f_alpha = 0.99

        for it in range(max_iters):
            forces, energy, max_ov, n_cont = self._compute_contacts_and_forces(coords, radii, r_max)

            if n_cont == 0 or (max_ov / r_mean) < tol_overlap_rel:
                return coords, max_ov, it

            P = float(np.sum(forces * v))
            f_norm = float(np.linalg.norm(forces))
            v_norm = float(np.linalg.norm(v))

            if f_norm > 1e-12 and v_norm > 1e-12:
                v = (1.0 - alpha) * v + alpha * (forces / f_norm) * v_norm

            if P > 0:
                n_pos += 1
                if n_pos > N_min:
                    dt = min(dt * f_inc, dt_max)
                    alpha *= f_alpha
            else:
                n_pos = 0
                v[:] = 0.0
                dt = max(dt * f_dec, dt_min)
                alpha = 0.10

            v += forces * dt
            coords = np.mod(coords + v * dt, self.box_size)

        return coords, max_ov, max_iters

    def generate(
        self,
        n_spheres: int = 1000,
        target_solid_fraction: float = 0.642,
        progress_callback: Optional[Callable[[int, float, float], None]] = None,
    ) -> GrainPack:
        """
        Executes staged inflation from dilute seed to dense jammed packing.

        Parameters
        ----------
        n_spheres : int, default=1000
            Number of spheres in the dense assembly.
        target_solid_fraction : float, default=0.642
            Target solid volume fraction (e.g. 0.642 for monodisperse RCP).
            Corresponding porosity is 1 - target_solid_fraction.
        progress_callback : callable, optional
            Receives (stage_number, current_solid_fraction, max_relative_overlap).

        Returns
        -------
        GrainPack
            Dense random close packing with strict zero overlap.
        """
        if self.psd is None:
            base_radii = np.ones(n_spheres, dtype=np.float64)
        else:
            base_radii = self.psd.sample(n_spheres, rng=self.rng)

        base_sphere_vol = np.sum(4.0 / 3.0 * np.pi * (base_radii ** 3))

        # 1. Initialize random positions at dilute density phi_0 = 0.30
        coords = self.rng.uniform(0.0, self.box_size, size=(n_spheres, 3))
        phi_curr = 0.30
        scale = ((phi_curr * self.volume) / base_sphere_vol) ** (1.0 / 3.0)
        radii = base_radii * scale

        coords, max_ov, iters = self._relax_fire(coords, radii, max_iters=500)

        # 2. Staged Inflation Loop
        stage = 0
        while phi_curr < target_solid_fraction:
            stage += 1
            if phi_curr < 0.50:
                dphi = 0.040
            elif phi_curr < 0.58:
                dphi = 0.020
            elif phi_curr < 0.62:
                dphi = 0.010
            elif phi_curr < 0.638:
                dphi = 0.005
            else:
                dphi = 0.002

            phi_next = min(phi_curr + dphi, target_solid_fraction)
            scale = ((phi_next * self.volume) / base_sphere_vol) ** (1.0 / 3.0)
            radii = base_radii * scale
            phi_curr = phi_next

            coords, max_ov, iters = self._relax_fire(coords, radii, max_iters=800, tol_overlap_rel=5e-5)

            if progress_callback:
                progress_callback(stage, phi_curr, max_ov / float(np.mean(radii)))

        # 3. Guaranteed Hard-Sphere Deflation (Zero Overlap Enforcement)
        tree = sp.cKDTree(coords, boxsize=self.box_size)
        pairs = tree.query_pairs(2.0 * float(np.max(radii)), output_type="ndarray")

        if len(pairs) > 0:
            i, j = pairs[:, 0], pairs[:, 1]
            dr = coords[j] - coords[i]
            dr -= self.box_size * np.round(dr / self.box_size)
            dist = np.linalg.norm(dr, axis=1)
            sigma = radii[i] + radii[j]
            overlap = sigma - dist
            actual_max_overlap = float(np.max(overlap)) if np.any(overlap > 0) else 0.0
        else:
            actual_max_overlap = 0.0

        if actual_max_overlap > 0.0:
            radii_hard = np.maximum(1e-6, radii - (actual_max_overlap / 2.0))
        else:
            radii_hard = radii.copy()

        return GrainPack(
            coords=coords,
            radii=radii_hard,
            box_size=self.box_size,
            periodic=(True, True, True),
            attributes={"dense_packing_stage": stage},
        )
