"""Unit tests for backward compatibility of genrandsp.py."""

import numpy as np
import pytest
from genrandsp import generate_spherepack


def test_legacy_wrapper_method_I_and_II():
    gsp = generate_spherepack(
        dimens=[50, 50, 50],
        Rs=[8, 5, 3],
        min_throat=0.5,
        overlapping=False,
        loop_limits=2000,
    )

    # Test method_I (target porosity)
    coords1 = gsp.method_I(porosity=0.85)
    assert isinstance(coords1, np.ndarray)
    assert coords1.ndim == 2
    assert coords1.shape[1] == 4
    assert len(coords1) > 0

    # Test method_II (fixed sphere count)
    coords2 = gsp.method_II(num_spheres=15)
    assert isinstance(coords2, np.ndarray)
    assert len(coords2) == 15
    assert coords2.shape[1] == 4
