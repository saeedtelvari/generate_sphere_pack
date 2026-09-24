"""
poropack.generators.sedimentation
=================================
3D Ballistic Gravitational Sedimentation ("Drop-and-Roll") of hard spheres.
Solves exact analytical trajectories of grains rolling into stable 3-point mechanical equilibrium
with Periodic Boundary Conditions in (x, y) and an impermeable floor at z = 0.
"""

from __future__ import annotations

from typing import Callable, List, Optional, Tuple, Union
import numpy as np

from poropack.core.psd import BasePSD, DiscretePSD
from poropack.core.representations import GrainPack


class SedimentationGenerator:
    """
    Simulates gravitational ballistic deposition (Drop-and-Roll) of hard spheres
    into stable 3-point mechanical equilibrium with Periodic Boundary Conditions in (x, y).

    Parameters
    ----------
    box_size : tuple of float
        (Lx, Ly, Lz) dimensions of the deposition domain.
    psd : BasePSD or sequence of float
        Particle size distribution or list of discrete radii.
    random_state : int or Generator, optional
        Seed or random generator.
    """

    def __init__(
        self,
        box_size: Tuple[float, float, float] | np.ndarray,
        psd: Union[BasePSD, Sequence[float]],
        random_state: Optional[int | np.random.Generator] = None,
        periodic_xy: bool = True,
    ):
        self.box_size = np.asarray(box_size, dtype=np.float64)
        self.Lx = float(self.box_size[0])
        self.Ly = float(self.box_size[1])
        self.Lz = float(self.box_size[2])
        self.periodic_xy = bool(periodic_xy)

        if isinstance(psd, BasePSD):
            self.psd = psd
        else:
            self.psd = DiscretePSD(psd)

        if isinstance(random_state, np.random.Generator):
            self.rng = random_state
        else:
            self.rng = np.random.default_rng(random_state)

        r_min, r_max = self.psd.get_bounds()
        self.r_max = float(r_max)
        self.r_min = float(r_min)

    def _pbc_disp(self, dx: float, dy: float) -> Tuple[float, float]:
        """Minimum image displacement in periodic (x, y)."""
        dx -= self.Lx * np.round(dx / self.Lx)
        dy -= self.Ly * np.round(dy / self.Ly)
        return float(dx), float(dy)

    def _commit_grain(
        self,
        coords: List[np.ndarray],
        radii: List[float],
        contacts: List[List[int]],
        grid: dict,
        dx_cell: float,
        dy_cell: float,
        Nx_cell: int,
        Ny_cell: int,
        xp: float,
        yp: float,
        zp: float,
        r_p: float,
        c_list: List[int],
    ) -> None:
        final_z = max(float(r_p), float(zp))
        if final_z == float(r_p) and -1 not in c_list:
            c_list.append(-1)
        new_idx = len(coords)
        coords.append(np.array([xp % self.Lx, yp % self.Ly, final_z], dtype=np.float64))
        radii.append(r_p)
        contacts.append(c_list)
        cx = int(np.floor((xp % self.Lx) / dx_cell)) % Nx_cell
        cy = int(np.floor((yp % self.Ly) / dy_cell)) % Ny_cell
        key = (cx, cy)
        if key not in grid:
            grid[key] = []
        grid[key].append(new_idx)

    def generate(
        self,
        target_porosity: Optional[float] = None,
        n_spheres: Optional[int] = None,
        max_hops_per_grain: int = 50,
        progress_callback: Optional[Callable[[int, float], None]] = None,
    ) -> GrainPack:
        """
        Deposit grains sequentially under gravity into stable multi-contact equilibrium.

        Parameters
        ----------
        target_porosity : float, optional
            Target void fraction to stop deposition (e.g. 0.40).
        n_spheres : int, optional
            Fixed number of grains to deposit.
        max_hops_per_grain : int, default=50
            Maximum rolling transitions before committing grain.
        progress_callback : callable, optional
            Callback receiving (n_deposited, current_porosity).

        Returns
        -------
        GrainPack
            Sedimented grain assembly with contact attributes.
        """
        if target_porosity is None and n_spheres is None:
            raise ValueError("Must specify either target_porosity or n_spheres")

        # 2D cell-linked list for fast column queries in (x, y)
        d_cell = 2.0 * self.r_max
        Nx_cell = max(3, int(np.floor(self.Lx / d_cell)))
        Ny_cell = max(3, int(np.floor(self.Ly / d_cell)))
        dx_cell = self.Lx / Nx_cell
        dy_cell = self.Ly / Ny_cell

        # Grid cells: (cx, cy) -> list of grain indices
        grid: dict[Tuple[int, int], List[int]] = {}

        def query_neighbors(qx: float, qy: float) -> List[int]:
            cx = int(np.floor(qx / dx_cell)) % Nx_cell
            cy = int(np.floor(qy / dy_cell)) % Ny_cell
            candidates: List[int] = []
            for ox in (-1, 0, 1):
                for oy in (-1, 0, 1):
                    key = ((cx + ox) % Nx_cell, (cy + oy) % Ny_cell)
                    if key in grid:
                        candidates.extend(grid[key])
            return candidates

        coords: List[np.ndarray] = []
        radii: List[float] = []
        contacts: List[List[int]] = []

        total_domain_vol = self.Lx * self.Ly * self.Lz
        current_solid_vol = 0.0

        def should_continue() -> bool:
            if n_spheres is not None:
                return len(coords) < n_spheres
            if target_porosity is not None:
                current_poro = 1.0 - (current_solid_vol / total_domain_vol)
                return current_poro > target_porosity
            return False

        while should_continue():
            r_p = float(self.psd.sample(1, rng=self.rng)[0])
            xp = float(self.rng.uniform(0.0, self.Lx))
            yp = float(self.rng.uniform(0.0, self.Ly))
            zp = self.Lz

            stage = 1  # 1: Vertical drop, 2: 1-point roll, 3: 2-point valley roll
            contact_A = -1
            contact_B = -1
            hops = 0
            committed = False

            while hops < max_hops_per_grain and not committed:
                hops += 1

                # ==========================================================
                # STAGE 1: Vertical Ballistic Drop
                # ==========================================================
                if stage == 1:
                    best_z = r_p  # Substrate contact level
                    best_idx = -1

                    for idx in query_neighbors(xp, yp):
                        ci = coords[idx]
                        ri = radii[idx]
                        dx, dy = self._pbc_disp(xp - ci[0], yp - ci[1])
                        dxy2 = dx * dx + dy * dy
                        R_sum = r_p + ri
                        if dxy2 < (R_sum * R_sum):
                            z_cand = ci[2] + np.sqrt(R_sum * R_sum - dxy2)
                            if z_cand > best_z:
                                best_z = z_cand
                                best_idx = idx

                    if best_idx == -1:
                        # Rests directly on flat substrate floor
                        self._commit_grain(
                            coords, radii, contacts, grid,
                            dx_cell, dy_cell, Nx_cell, Ny_cell,
                            xp, yp, r_p, r_p, [-1]
                        )
                        committed = True
                        break
                    else:
                        contact_A = best_idx
                        cA = coords[contact_A]
                        dx, dy = self._pbc_disp(xp - cA[0], yp - cA[1])
                        if (dx * dx + dy * dy) < 1e-12:
                            dx = 1e-6 * r_p  # Apex perturbation
                        xp = cA[0] + dx
                        yp = cA[1] + dy
                        zp = best_z
                        stage = 2
                        continue

                # ==========================================================
                # STAGE 2: 1-Point Steepest Descent Roll on Sphere A
                # ==========================================================
                elif stage == 2:
                    cA = coords[contact_A]
                    rA = radii[contact_A]
                    RPA = r_p + rA

                    dx, dy = self._pbc_disp(xp - cA[0], yp - cA[1])
                    dz = zp - cA[2]
                    phi = np.arctan2(dy, dx)
                    theta0 = float(np.arccos(np.clip(dz / RPA, -1.0, 1.0)))

                    event_theta = np.pi / 2.0  # Equator detachment
                    event_type = "equator"
                    event_target = -1

                    # Check substrate floor collision along meridian
                    if (cA[2] - r_p) < RPA and (cA[2] - r_p) > -RPA:
                        cos_th_f = (r_p - cA[2]) / RPA
                        if -1.0 <= cos_th_f <= 1.0:
                            th_f = float(np.arccos(cos_th_f))
                            if theta0 < th_f < event_theta:
                                event_theta = th_f
                                event_type = "floor"

                    # Intersect neighbor spheres along meridian path
                    for idx in query_neighbors(xp, yp):
                        if idx == contact_A:
                            continue
                        ci = coords[idx]
                        ri = radii[idx]
                        RPB = r_p + ri
                        dx_AB, dy_AB = self._pbc_disp(cA[0] - ci[0], cA[1] - ci[1])
                        dz_AB = cA[2] - ci[2]
                        d_vec = np.array([dx_AB, dy_AB, dz_AB])
                        d2 = float(np.sum(d_vec ** 2))

                        if d2 > (RPA + RPB) ** 2 or d2 < (RPA - RPB) ** 2:
                            continue

                        K = (RPB ** 2 - RPA ** 2 - d2) / (2.0 * RPA)
                        A_c = d_vec[0] * np.cos(phi) + d_vec[1] * np.sin(phi)
                        B_c = d_vec[2]
                        R_amp = np.sqrt(A_c * A_c + B_c * B_c)

                        if R_amp < 1e-12 or abs(K) > R_amp:
                            continue

                        ratio = float(np.clip(K / R_amp, -1.0, 1.0))
                        psi = float(np.arctan2(A_c, B_c))
                        acos_val = float(np.arccos(ratio))

                        for th in ((psi + acos_val) % (2.0 * np.pi), (psi - acos_val) % (2.0 * np.pi)):
                            if theta0 + 1e-6 < th < event_theta:
                                event_theta = th
                                event_type = "sphere"
                                event_target = idx

                    # Update position along meridian arc
                    xp = cA[0] + RPA * np.sin(event_theta) * np.cos(phi)
                    yp = cA[1] + RPA * np.sin(event_theta) * np.sin(phi)
                    zp = cA[2] + RPA * np.cos(event_theta)

                    if event_type == "floor" or zp <= r_p:
                        self._commit_grain(
                            coords, radii, contacts, grid,
                            dx_cell, dy_cell, Nx_cell, Ny_cell,
                            xp, yp, r_p, r_p, [contact_A, -1]
                        )
                        committed = True
                        break
                    elif event_type == "equator":
                        stage = 1  # Detach and free-fall vertically
                        continue
                    elif event_type == "sphere":
                        contact_B = event_target
                        stage = 3  # Crevasse / valley roll
                        continue

                # ==========================================================
                # STAGE 3: 2-Point Valley Roll on Grains A and B
                # ==========================================================
                elif stage == 3:
                    cA = coords[contact_A]
                    rA = radii[contact_A]
                    cB = coords[contact_B]
                    rB = radii[contact_B]
                    RPA = r_p + rA
                    RPB = r_p + rB

                    dx_AB, dy_AB = self._pbc_disp(cB[0] - cA[0], cB[1] - cA[1])
                    dz_AB = cB[2] - cA[2]
                    D_vec = np.array([dx_AB, dy_AB, dz_AB])
                    D = float(np.linalg.norm(D_vec))

                    if D < 1e-12 or D > (RPA + RPB):
                        stage = 1
                        continue

                    n_AB = D_vec / D
                    dA = (D * D + RPA * RPA - RPB * RPB) / (2.0 * D)
                    rho2 = RPA * RPA - dA * dA
                    if rho2 <= 0.0:
                        stage = 1
                        continue

                    rho = np.sqrt(rho2)
                    c_circ = cA + dA * n_AB

                    # Gravity projection onto circular locus
                    k_down = np.array([0.0, 0.0, -1.0])
                    v_down = k_down - np.dot(k_down, n_AB) * n_AB
                    norm_v = float(np.linalg.norm(v_down))
                    if norm_v < 1e-8:
                        stage = 1
                        continue

                    u_down = v_down / norm_v
                    u_perp = np.cross(n_AB, u_down)
                    u_perp /= float(np.linalg.norm(u_perp))

                    # Angle on circle alpha0
                    dx_p, dy_p = self._pbc_disp(xp - c_circ[0], yp - c_circ[1])
                    dz_p = zp - c_circ[2]
                    r_rel = np.array([dx_p, dy_p, dz_p])
                    alpha0 = float(np.arctan2(np.dot(r_rel, u_perp), np.dot(r_rel, u_down)))

                    event_alpha = 0.0  # Valley minimum
                    event_type = "bottom"
                    event_target = -1

                    # Intersect candidate 3rd grains C
                    for idx in query_neighbors(xp, yp):
                        if idx in (contact_A, contact_B):
                            continue
                        ci = coords[idx]
                        ri = radii[idx]
                        RPC = r_p + ri

                        dx_C, dy_C = self._pbc_disp(c_circ[0] - ci[0], c_circ[1] - ci[1])
                        dz_C = c_circ[2] - ci[2]
                        dC = np.array([dx_C, dy_C, dz_C])
                        dC2 = float(np.sum(dC ** 2))

                        if dC2 > (rho + RPC) ** 2:
                            continue

                        KC = (RPC ** 2 - rho ** 2 - dC2) / (2.0 * rho)
                        AC = float(np.dot(dC, u_down))
                        BC = float(np.dot(dC, u_perp))
                        RC_amp = np.sqrt(AC * AC + BC * BC)

                        if RC_amp < 1e-12 or abs(KC) > RC_amp:
                            continue

                        ratio = float(np.clip(KC / RC_amp, -1.0, 1.0))
                        psiC = float(np.arctan2(BC, AC))
                        acos_val = float(np.arccos(ratio))

                        for a in (psiC + acos_val, psiC - acos_val):
                            a = (a + np.pi) % (2.0 * np.pi) - np.pi
                            if alpha0 > 0.0 and (1e-6 < a < alpha0):
                                if a > event_alpha:
                                    event_alpha = a
                                    event_type = "sphere"
                                    event_target = idx
                            elif alpha0 < 0.0 and (alpha0 < a < -1e-6):
                                if a < event_alpha:
                                    event_alpha = a
                                    event_type = "sphere"
                                    event_target = idx

                    # Update position on circle
                    pos = c_circ + rho * (np.cos(event_alpha) * u_down + np.sin(event_alpha) * u_perp)
                    xp, yp, zp = float(pos[0]), float(pos[1]), float(pos[2])

                    # Floor boundary check
                    if zp <= r_p:
                        self._commit_grain(
                            coords, radii, contacts, grid,
                            dx_cell, dy_cell, Nx_cell, Ny_cell,
                            xp, yp, r_p, r_p, [contact_A, contact_B, -1]
                        )
                        committed = True
                        break

                    if event_type == "sphere":
                        contact_C = event_target
                        # Mechanical force equilibrium check
                        vA = self._unit_normal(xp, yp, zp, coords[contact_A])
                        vB = self._unit_normal(xp, yp, zp, coords[contact_B])
                        vC = self._unit_normal(xp, yp, zp, coords[contact_C])

                        M = np.column_stack([vA, vB, vC])
                        if abs(float(np.linalg.det(M))) > 1e-6:
                            forces = np.linalg.solve(M, np.array([0.0, 0.0, 1.0]))
                            if np.all(forces > 1e-4):
                                # Stable 3-point tripod
                                self._commit_grain(
                                    coords, radii, contacts, grid,
                                    dx_cell, dy_cell, Nx_cell, Ny_cell,
                                    xp, yp, zp, r_p, [contact_A, contact_B, contact_C]
                                )
                                committed = True
                                break
                            else:
                                min_idx = int(np.argmin(forces))
                                if min_idx == 0:
                                    contact_A = contact_C
                                elif min_idx == 1:
                                    contact_B = contact_C
                                stage = 3
                                continue
                        else:
                            contact_A = contact_C
                            stage = 3
                            continue
                    else:
                        # Reached valley bottom without third contact
                        self._commit_grain(
                            coords, radii, contacts, grid,
                            dx_cell, dy_cell, Nx_cell, Ny_cell,
                            xp, yp, zp, r_p, [contact_A, contact_B]
                        )
                        committed = True
                        break

            if not committed:
                # Commit at current position as fallback
                self._commit_grain(
                    coords, radii, contacts, grid,
                    dx_cell, dy_cell, Nx_cell, Ny_cell,
                    xp, yp, zp, r_p, [contact_A] if contact_A >= 0 else [-1]
                )

            current_solid_vol += 4.0 / 3.0 * np.pi * (r_p ** 3)
            if progress_callback and len(coords) % 500 == 0:
                poro = 1.0 - (current_solid_vol / total_domain_vol)
                progress_callback(len(coords), poro)

        out_coords = np.array(coords, dtype=np.float64) if len(coords) > 0 else np.empty((0, 3))
        out_radii = np.array(radii, dtype=np.float64) if len(radii) > 0 else np.empty(0)

        return GrainPack(
            coords=out_coords,
            radii=out_radii,
            box_size=self.box_size,
            periodic=(self.periodic_xy, self.periodic_xy, False),  # Periodic in (x, y), closed floor at z=0
            attributes={"contacts": contacts},
        )

    def _unit_normal(self, xp: float, yp: float, zp: float, target_c: np.ndarray) -> np.ndarray:
        dx, dy = self._pbc_disp(xp - target_c[0], yp - target_c[1])
        dz = zp - target_c[2]
        vec = np.array([dx, dy, dz], dtype=np.float64)
        norm = float(np.linalg.norm(vec))
        return vec / max(1e-12, norm)
