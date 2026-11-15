import numpy as np
import pytest
from higgs_review.kinematics import diphoton_mass, acceptance_mask


def test_back_to_back_central():
    assert np.isclose(diphoton_mass(62.5,0,0,62.5,0,np.pi),125.0,atol=1e-9)


def test_exchange_and_global_rotation_invariant():
    a = diphoton_mass(50,0.3,0.2,38,-0.8,2.9)
    b = diphoton_mass(38,-0.8,2.9+0.8,50,0.3,0.2+0.8)
    assert np.isclose(a,b)


def test_zero_collinear_mass():
    assert diphoton_mass(35,0,0,35,0,0) == 0


def test_nonnegative_validation():
    with pytest.raises(ValueError):
        diphoton_mass(-1,0,0,50,0,1)


def test_fiducial_crack_veto():
    assert bool(acceptance_mask(60,50,0.2,0.3,125))
    assert not bool(acceptance_mask(60,50,1.44,0.3,125))


def test_relative_pt_cut():
    assert not bool(acceptance_mask(30,30,0.2,0.3,125))
