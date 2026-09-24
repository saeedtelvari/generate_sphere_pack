"""Unit tests for all poropack generation engines."""

import numpy as np
import pytest
from poropack.core.psd import DiscretePSD, LogNormalPSD
from poropack.generators.rsa import RSAGenerator
from poropack.generators.sedimentation import SedimentationGenerator
from poropack.generators.dense_packing import DensePackingGenerator
from poropack.generators.continuum import GaussianRandomField, QSGSGenerator


def test_rsa_generator():
    psd = DiscretePSD([6.0, 4.0])
    gen = RSAGenerator(box_size=(60.0, 60.0, 60.0), psd=psd, random_state=42)

    pack = gen.generate(target_porosity=0.75, max_attempts=5000)
    assert len(pack) > 10
    assert pack.analytical_porosity() <= 0.75

    # Check zero overlap among spheres
    coords = pack.coords
    radii = pack.radii
    L = pack.box_size
    for i in range(len(pack)):
        for j in range(i + 1, len(pack)):
            dr = coords[j] - coords[i]
            # Minimum image
            dr -= L * np.round(dr / L)
            dist = np.linalg.norm(dr)
            assert dist >= (radii[i] + radii[j] - 1e-8)


def test_sedimentation_generator():
    psd = DiscretePSD([4.0])
    sed = SedimentationGenerator(box_size=(40.0, 40.0, 60.0), psd=psd, random_state=42)

    # Deposit 40 spheres
    pack = sed.generate(n_spheres=40)
    assert len(pack) == 40
    assert "contacts" in pack.attributes

    # Check that spheres don't fall below substrate z = r
    assert np.all(pack.coords[:, 2] >= pack.radii - 1e-8)


def test_dense_packing_generator():
    # Small test assembly to test FIRE convergence and breaking the RSA limit
    gen = DensePackingGenerator(
        box_size=(30.0, 30.0, 30.0),
        random_state=42,
    )
    # Inflate to solid fraction 0.640 (porosity 0.360)
    pack = gen.generate(n_spheres=120, target_solid_fraction=0.640)
    assert len(pack) == 120

    # Strict hard sphere check (no overlap)
    coords = pack.coords
    radii = pack.radii
    L = pack.box_size
    for i in range(len(pack)):
        for j in range(i + 1, len(pack)):
            dr = coords[j] - coords[i]
            dr -= L * np.round(dr / L)
            dist = np.linalg.norm(dr)
            assert dist >= (radii[i] + radii[j] - 1e-6)

    # Solid fraction far exceeds RSA jamming limit (~0.35) and achieves dense packing
    phi_solid = 1.0 - pack.analytical_porosity()
    assert phi_solid >= 0.58


def test_gaussian_random_field():
    grf = GaussianRandomField(
        shape=(32, 32, 32),
        voxel_size=1.0,
        correlation_length=5.0,
        kernel="gaussian",
        random_state=42,
    )
    grid = grf.generate(target_porosity=0.30)
    assert grid.shape == (32, 32, 32)
    # Exact analytical thresholding should match target porosity within small discrete sampling tolerance
    assert abs(grid.porosity() - 0.30) < 0.02


def test_qsgs_generator():
    qsgs = QSGSGenerator(
        shape=(30, 30, 30),
        voxel_size=1.0,
        core_probability=0.02,
        growth_probabilities=(0.15, 0.15, 0.05),  # anisotropic (P_x, P_y > P_z)
        random_state=42,
    )
    grid = qsgs.generate(target_porosity=0.40)
    assert grid.shape == (30, 30, 30)
    assert abs(grid.porosity() - 0.40) < 0.05
