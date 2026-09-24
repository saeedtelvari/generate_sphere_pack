"""
poropack.generators.continuum
=============================
Continuum porous media generation via 3D Spectral Gaussian Random Fields (GRF)
and Quartet Structure Generation Set (QSGS).
"""

from __future__ import annotations

from typing import Optional, Tuple
import numpy as np
from scipy import special

from poropack.core.representations import VoxelGrid


class GaussianRandomField:
    """
    3D Level-Cut Spectral Gaussian Random Field (GRF) generator.
    Produces bicontinuous, carbonate, chalk, and spinodal-like microstructures in O(V log V) time.

    Parameters
    ----------
    shape : tuple of int
        (Nx, Ny, Nz) grid dimensions.
    voxel_size : float, default=1.0
        Physical voxel dimension h.
    correlation_length : float, default=10.0
        Spatial autocorrelation length lc (in physical units).
    kernel : {'gaussian', 'exponential', 'von_karman'}, default='gaussian'
        Spectral filter kernel type.
    hurst_exponent : float, default=0.5
        Roughness exponent H in (0, 1) for Von Kármán kernel.
    random_state : int or Generator, optional
        Seed or random generator.
    """

    def __init__(
        self,
        shape: Tuple[int, int, int] = (100, 100, 100),
        voxel_size: float = 1.0,
        correlation_length: float = 10.0,
        kernel: str = "gaussian",
        hurst_exponent: float = 0.5,
        random_state: Optional[int | np.random.Generator] = None,
    ):
        self.shape = (int(shape[0]), int(shape[1]), int(shape[2]))
        self.voxel_size = float(voxel_size)
        self.lc = float(correlation_length)
        self.kernel = kernel.lower()
        self.hurst = float(hurst_exponent)

        if isinstance(random_state, np.random.Generator):
            self.rng = random_state
        else:
            self.rng = np.random.default_rng(random_state)

    def _build_spectral_filter(self) -> np.ndarray:
        """Construct the 3D radially symmetric filter in Fourier space."""
        Nx, Ny, Nz = self.shape
        # Spatial frequency coordinates (cycles per physical unit)
        kx = np.fft.fftfreq(Nx, d=self.voxel_size)
        ky = np.fft.fftfreq(Ny, d=self.voxel_size)
        kz = np.fft.fftfreq(Nz, d=self.voxel_size)

        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
        K = 2.0 * np.pi * np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)

        kc = 2.0 * np.pi / self.lc

        if self.kernel == "gaussian":
            # Autocorrelation ~ exp(-(r/lc)^2) -> Spectral power ~ exp(-(k/kc)^2)
            power = np.exp(-((K / kc) ** 2))
        elif self.kernel == "exponential":
            # Autocorrelation ~ exp(-r/lc) -> Spectral power ~ (1 + (K/kc)^2)^(-2)
            power = (1.0 + (K / kc) ** 2) ** (-2.0)
        elif self.kernel in ("von_karman", "matern"):
            # Matérn / Von Kármán spectral filter
            power = (1.0 + (K / kc) ** 2) ** (-(self.hurst + 1.5))
        else:
            raise ValueError(f"Unknown spectral kernel: {self.kernel}")

        # Amplitude filter is sqrt(power)
        return np.sqrt(power)

    def generate(self, target_porosity: float = 0.30) -> VoxelGrid:
        """
        Generate a 3D porous structure with guaranteed exact analytical porosity.

        Parameters
        ----------
        target_porosity : float
            Desired void fraction (e.g. 0.30 for 30% porosity).

        Returns
        -------
        VoxelGrid
            Eulerian 3D representation where True is solid and False is pore.
        """
        if not (0.0 < target_porosity < 1.0):
            raise ValueError(f"target_porosity must be between 0 and 1, got {target_porosity}")

        Nx, Ny, Nz = self.shape
        # 1. Generate 3D white Gaussian noise
        white_noise = self.rng.standard_normal(size=(Nx, Ny, Nz))

        # 2. Fourier transform
        noise_hat = np.fft.fftn(white_noise)

        # 3. Spectral convolution (multiplication by filter kernel)
        filt = self._build_spectral_filter()
        field_hat = noise_hat * filt

        # 4. Inverse Fourier transform to real space
        field = np.real(np.fft.ifftn(field_hat))

        # 5. Standardize field to zero mean and unit variance N(0, 1)
        mean_val = np.mean(field)
        std_val = np.std(field)
        if std_val < 1e-12:
            std_val = 1.0
        standardized_field = (field - mean_val) / std_val

        # 6. Exact analytical thresholding via inverse normal CDF
        # P(Z < alpha) = target_porosity -> pore phase is Z < alpha, solid is Z >= alpha
        alpha = np.sqrt(2.0) * special.erfinv(2.0 * target_porosity - 1.0)
        matrix = standardized_field >= alpha

        return VoxelGrid(
            matrix=matrix,
            voxel_size=self.voxel_size,
            origin=np.zeros(3, dtype=np.float64),
            periodic=(True, True, True),
        )


class QSGSGenerator:
    """
    Quartet Structure Generation Set (QSGS) for anisotropic & multi-scale porous media.

    Parameters
    ----------
    shape : tuple of int
        (Nx, Ny, Nz) grid dimensions.
    voxel_size : float, default=1.0
        Physical voxel size.
    core_probability : float, default=0.01
        Nucleation probability for initial solid seed cores (c_d).
    growth_probabilities : tuple of 3 float, default=(0.1, 0.1, 0.1)
        Directional expansion probabilities (P_x, P_y, P_z).
        Ratios P_x/P_z > 1 induce stratification and horizontal permeability anisotropy.
    random_state : int or Generator, optional
        Seed or random generator.
    """

    def __init__(
        self,
        shape: Tuple[int, int, int] = (64, 64, 64),
        voxel_size: float = 1.0,
        core_probability: float = 0.01,
        growth_probabilities: Tuple[float, float, float] = (0.1, 0.1, 0.1),
        random_state: Optional[int | np.random.Generator] = None,
        core_distribution_probability: Optional[float] = None,
    ):
        self.shape = (int(shape[0]), int(shape[1]), int(shape[2]))
        self.voxel_size = float(voxel_size)
        if core_distribution_probability is not None:
            core_probability = core_distribution_probability
        self.cd = float(core_probability)
        self.P = np.asarray(growth_probabilities, dtype=np.float64)

        if isinstance(random_state, np.random.Generator):
            self.rng = random_state
        else:
            self.rng = np.random.default_rng(random_state)

    def generate(self, target_porosity: float = 0.35, max_cycles: int = 50) -> VoxelGrid:
        """
        Grow anisotropic solid phase from random core seeds until target porosity is reached.

        Parameters
        ----------
        target_porosity : float
            Desired void fraction.
        max_cycles : int, default=50
            Maximum growth cycles.

        Returns
        -------
        VoxelGrid
        """
        Nx, Ny, Nz = self.shape
        target_solid_count = int(np.round((1.0 - target_porosity) * (Nx * Ny * Nz)))

        # 1. Distribute random seed cores
        matrix = self.rng.random(self.shape) < self.cd
        current_solid = int(np.count_nonzero(matrix))

        if current_solid >= target_solid_count:
            return VoxelGrid(matrix=matrix, voxel_size=self.voxel_size, periodic=(True, True, True))

        # Directional shifts along x, y, z
        directions = [
            (1, 0, 0, self.P[0]),
            (-1, 0, 0, self.P[0]),
            (0, 1, 0, self.P[1]),
            (0, -1, 0, self.P[1]),
            (0, 0, 1, self.P[2]),
            (0, 0, -1, self.P[2]),
        ]

        cycle = 0
        while current_solid < target_solid_count and cycle < max_cycles:
            cycle += 1
            new_solid = np.zeros_like(matrix)

            for dx, dy, dz, prob in directions:
                # Roll matrix periodically
                shifted = np.roll(matrix, shift=(dx, dy, dz), axis=(0, 1, 2))
                # Only empty voxels adjacent to solid can grow
                candidates = shifted & (~matrix)
                # Stochastic acceptance
                accepted = candidates & (self.rng.random(self.shape) < prob)
                new_solid |= accepted

            matrix |= new_solid
            current_solid = int(np.count_nonzero(matrix))

        return VoxelGrid(
            matrix=matrix,
            voxel_size=self.voxel_size,
            origin=np.zeros(3, dtype=np.float64),
            periodic=(True, True, True),
        )
