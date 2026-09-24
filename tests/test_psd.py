"""Unit tests for continuous geological PSDs."""

import numpy as np
import pytest
from poropack.core.psd import (
    DiscretePSD,
    EmpiricalSievePSD,
    LogNormalPSD,
    TruncatedGaussianPSD,
    WeibullPSD,
)


def test_log_normal_psd():
    psd_vol = LogNormalPSD(d50=20.0, sigma_phi=0.4, volume_weighted=True)
    psd_num = LogNormalPSD(d50=20.0, sigma_phi=0.4, volume_weighted=False)

    samples_vol = psd_vol.sample(1000, rng=np.random.default_rng(42))
    samples_num = psd_num.sample(1000, rng=np.random.default_rng(42))

    assert len(samples_vol) == 1000
    assert len(samples_num) == 1000

    # Due to Hatch-Choate frequency shift, volume-weighted input sampled by number
    # shifts to smaller radii compared to purely number-weighted
    assert np.median(samples_vol) < np.median(samples_num)

    r_min, r_max = psd_vol.get_bounds()
    assert np.all(samples_vol >= r_min)
    assert np.all(samples_vol <= r_max)


def test_weibull_and_gaussian_psd():
    w = WeibullPSD(scale=10.0, shape=2.5)
    samples_w = w.sample(500)
    assert len(samples_w) == 500
    assert np.all(samples_w > 0)

    g = TruncatedGaussianPSD(mean_radius=15.0, std_radius=2.0)
    samples_g = g.sample(500)
    assert len(samples_g) == 500
    assert np.all(samples_g >= g.bounds[0])
    assert np.all(samples_g <= g.bounds[1])


def test_discrete_and_empirical_psd():
    d = DiscretePSD(radii=[12.0, 8.0, 4.0])
    samples_d = d.sample(100)
    for s in samples_d:
        assert s in (12.0, 8.0, 4.0)

    # Sieve sizes: 10, 20, 50 um with 10%, 50%, 90% passing
    sieve = EmpiricalSievePSD(diameters=[10.0, 20.0, 50.0], percent_passing=[10.0, 50.0, 90.0])
    samples_s = sieve.sample(200)
    assert len(samples_s) == 200
    assert np.all(samples_s >= 5.0)  # radius = 10 / 2
    assert np.all(samples_s <= 25.0) # radius = 50 / 2
