"""
poropack.core.psd
=================
Continuous geological and synthetic Particle Size Distributions (PSD).
Supports volume-weighted vs. number-weighted analytical conversions (Hatch-Choate).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional, Sequence, Tuple
import numpy as np
from scipy import special, interpolate


class BasePSD(ABC):
    """Abstract base class for grain size distributions."""

    @abstractmethod
    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Sample grain radii from the distribution."""
        pass

    @abstractmethod
    def get_bounds(self) -> Tuple[float, float]:
        """Return (r_min, r_max) bounds of the distribution."""
        pass


class LogNormalPSD(BasePSD):
    """
    Log-Normal Particle Size Distribution.

    Parameters
    ----------
    d50 : float
        Median grain diameter.
    sigma_phi : float, default=0.5
        Folk-Ward sedimentological sorting coefficient in phi units (phi = -log2(d)).
        - Well sorted: sigma_phi < 0.50
        - Moderately sorted: 0.71 <= sigma_phi <= 1.00
        - Poorly sorted: sigma_phi > 1.00
    volume_weighted : bool, default=True
        Whether the input d50 represents a volume-weighted sieve distribution (common in geology).
        If True, the distribution is converted to a number-weighted sampling density
        via the analytical Hatch-Choate equation: ln(d50_N) = ln(d50_V) - 3 * (ln sigma_g)^2.
    bounds : tuple of float, optional
        (r_min, r_max) truncation bounds for grain radii.
    """

    def __init__(
        self,
        d50: float,
        sigma_phi: float = 0.5,
        volume_weighted: bool = True,
        bounds: Optional[Tuple[float, float]] = None,
    ):
        self.d50_input = float(d50)
        self.sigma_phi = float(sigma_phi)
        self.volume_weighted = bool(volume_weighted)

        # Geometric standard deviation: sigma_g = 2^(sigma_phi)
        self.sigma_g = 2.0 ** self.sigma_phi
        self.ln_sigma = np.log(self.sigma_g)

        # Convert volume-weighted median to number-weighted median
        if self.volume_weighted:
            self.mu_ln_d = np.log(self.d50_input) - 3.0 * (self.ln_sigma ** 2)
        else:
            self.mu_ln_d = np.log(self.d50_input)

        # Convert diameter to radius parameter: ln(r) = ln(d/2) = ln(d) - ln(2)
        self.mu_ln_r = self.mu_ln_d - np.log(2.0)
        self.sigma_ln_r = self.ln_sigma

        # Default bounds: 3.5 geometric standard deviations around median
        if bounds is None:
            r_median = np.exp(self.mu_ln_r)
            r_min = max(1e-4, r_median / (self.sigma_g ** 3.5))
            r_max = r_median * (self.sigma_g ** 3.5)
            self.bounds = (float(r_min), float(r_max))
        else:
            self.bounds = (float(bounds[0]), float(bounds[1]))

    def get_bounds(self) -> Tuple[float, float]:
        return self.bounds

    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        """Sample grain radii via truncated log-normal distribution."""
        if rng is None:
            rng = np.random.default_rng()

        r_min, r_max = self.bounds
        # Compute CDF bounds
        cdf_min = 0.5 * (1.0 + special.erf((np.log(r_min) - self.mu_ln_r) / (self.sigma_ln_r * np.sqrt(2.0))))
        cdf_max = 0.5 * (1.0 + special.erf((np.log(r_max) - self.mu_ln_r) / (self.sigma_ln_r * np.sqrt(2.0))))

        # Inverse transform sampling within [cdf_min, cdf_max]
        u = rng.uniform(cdf_min, cdf_max, size=size)
        ln_r = self.mu_ln_r + self.sigma_ln_r * np.sqrt(2.0) * special.erfinv(2.0 * u - 1.0)
        return np.exp(ln_r)


class WeibullPSD(BasePSD):
    """
    Weibull / Rosin-Rammler Particle Size Distribution.

    Parameters
    ----------
    scale : float
        Weibull characteristic radius scale parameter (lambda).
    shape : float, default=2.0
        Weibull shape / uniformity parameter (k).
    bounds : tuple of float, optional
        (r_min, r_max) truncation bounds.
    """

    def __init__(
        self,
        scale: float,
        shape: float = 2.0,
        bounds: Optional[Tuple[float, float]] = None,
    ):
        self.scale = float(scale)
        self.shape = float(shape)
        if bounds is None:
            r_min = self.scale * ((-np.log(1.0 - 0.001)) ** (1.0 / self.shape))
            r_max = self.scale * ((-np.log(1.0 - 0.999)) ** (1.0 / self.shape))
            self.bounds = (float(r_min), float(r_max))
        else:
            self.bounds = (float(bounds[0]), float(bounds[1]))

    def get_bounds(self) -> Tuple[float, float]:
        return self.bounds

    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        if rng is None:
            rng = np.random.default_rng()

        r_min, r_max = self.bounds
        cdf_min = 1.0 - np.exp(-((r_min / self.scale) ** self.shape))
        cdf_max = 1.0 - np.exp(-((r_max / self.scale) ** self.shape))

        u = rng.uniform(cdf_min, cdf_max, size=size)
        return self.scale * ((-np.log(1.0 - u)) ** (1.0 / self.shape))


class TruncatedGaussianPSD(BasePSD):
    """
    Truncated Gaussian Particle Size Distribution.

    Parameters
    ----------
    mean_radius : float
        Mean sphere radius mu.
    std_radius : float
        Standard deviation sigma.
    bounds : tuple of float, optional
        (r_min, r_max) truncation bounds.
    """

    def __init__(
        self,
        mean_radius: float,
        std_radius: float,
        bounds: Optional[Tuple[float, float]] = None,
    ):
        self.mean = float(mean_radius)
        self.std = float(std_radius)
        if bounds is None:
            self.bounds = (max(1e-4, self.mean - 3.0 * self.std), self.mean + 3.0 * self.std)
        else:
            self.bounds = (float(bounds[0]), float(bounds[1]))

    def get_bounds(self) -> Tuple[float, float]:
        return self.bounds

    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        if rng is None:
            rng = np.random.default_rng()

        r_min, r_max = self.bounds
        cdf_min = 0.5 * (1.0 + special.erf((r_min - self.mean) / (self.std * np.sqrt(2.0))))
        cdf_max = 0.5 * (1.0 + special.erf((r_max - self.mean) / (self.std * np.sqrt(2.0))))

        u = rng.uniform(cdf_min, cdf_max, size=size)
        return self.mean + self.std * np.sqrt(2.0) * special.erfinv(2.0 * u - 1.0)


class DiscretePSD(BasePSD):
    """
    Discrete radius distribution (e.g. for backward compatibility with [12, 8, 4]).

    Parameters
    ----------
    radii : sequence of float
        Allowed discrete sphere radii.
    probabilities : sequence of float, optional
        Relative or normalized probabilities for each radius.
    """

    def __init__(
        self,
        radii: Sequence[float],
        probabilities: Optional[Sequence[float]] = None,
    ):
        self.radii = np.asarray(radii, dtype=np.float64)
        if len(self.radii) == 0:
            raise ValueError("radii must contain at least one value")

        if probabilities is None:
            self.probs = np.ones(len(self.radii), dtype=np.float64) / len(self.radii)
        else:
            p = np.asarray(probabilities, dtype=np.float64)
            self.probs = p / np.sum(p)

    def get_bounds(self) -> Tuple[float, float]:
        return float(np.min(self.radii)), float(np.max(self.radii))

    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        if rng is None:
            rng = np.random.default_rng()
        return rng.choice(self.radii, size=size, p=self.probs)


class EmpiricalSievePSD(BasePSD):
    """
    Empirical cumulative sieve distribution interpolated from laboratory data.

    Parameters
    ----------
    diameters : sequence of float
        Sieve mesh opening sizes (e.g. in microns or mm).
    percent_passing : sequence of float
        Cumulative percentage of mass/volume passing through each sieve (0 to 100).
    volume_weighted : bool, default=False
        If True, converts volume-weighted cumulative sieve curve to number frequency
        via f_N(r) proportional to f_V(r) / r^3 before constructing the inverse CDF.
    """

    def __init__(
        self,
        diameters: Sequence[float],
        percent_passing: Sequence[float],
        volume_weighted: bool = False,
    ):
        d = np.asarray(diameters, dtype=np.float64)
        p = np.asarray(percent_passing, dtype=np.float64) / 100.0

        sort_idx = np.argsort(d)
        d_sorted = d[sort_idx]
        p_sorted = p[sort_idx]

        # Convert diameter to radius
        r_sorted = d_sorted / 2.0
        self.bounds = (float(r_sorted[0]), float(r_sorted[-1]))

        if volume_weighted and len(r_sorted) >= 3:
            # Fine discretization across radius range
            r_fine = np.linspace(r_sorted[0], r_sorted[-1], 500)
            f_v_interp = interpolate.PchipInterpolator(r_sorted, p_sorted)
            F_V = f_v_interp(r_fine)
            f_V_density = np.maximum(0.0, np.gradient(F_V, r_fine))
            # Number density: f_N(r) proportional to f_V(r) / r^3
            f_N_density = f_V_density / (r_fine ** 3)
            F_N = np.cumsum(f_N_density) * (r_fine[1] - r_fine[0])
            if F_N[-1] > 0:
                F_N /= F_N[-1]
            self.p_min = float(F_N[0])
            self.p_max = float(F_N[-1])
            self.inv_cdf = interpolate.PchipInterpolator(F_N, r_fine)
        else:
            self.p_min = float(p_sorted[0])
            self.p_max = float(p_sorted[-1])
            self.inv_cdf = interpolate.PchipInterpolator(p_sorted, r_sorted)

    def get_bounds(self) -> Tuple[float, float]:
        return self.bounds

    def sample(self, size: int = 1, rng: Optional[np.random.Generator] = None) -> np.ndarray:
        if rng is None:
            rng = np.random.default_rng()
        u = rng.uniform(self.p_min, self.p_max, size=size)
        samples = self.inv_cdf(u)
        return np.clip(samples, self.bounds[0], self.bounds[1])
