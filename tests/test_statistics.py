import numpy as np
import pytest
from higgs_review.statistics import background_bin_counts, signal_bin_fractions, fit_local_signal, sideband_mask
from higgs_review.simulation import toy_pseudodata, exponential_samples

EDGES = np.linspace(100,160,61)


def test_constant_rate_exact_integral():
    b = background_bin_counts(EDGES, np.log(13.0), 0)
    assert np.allclose(b, 13.0)


def test_exponential_integral():
    b = background_bin_counts(np.array([100.,160.]), np.log(2), -.02)
    exact = 2*np.exp(-.02*(100-125))*np.expm1(-.02*60)/(-.02)
    assert np.isclose(b[0],exact)


def test_signal_template_normalized():
    sig = signal_bin_fractions(EDGES,125,1.7)
    assert np.isclose(sig.sum(),1.0)
    assert np.argmax(sig) in (24,25)


def test_bad_input_rejected():
    with pytest.raises(ValueError):
        background_bin_counts([100,102,101],0,0)
    with pytest.raises(ValueError):
        fit_local_signal([-1]*60,EDGES)


def test_null_vs_alt_reasonable_synthetic():
    n, _ = np.histogram(toy_pseudodata(background=8000,signal=600,seed=10),bins=EDGES)
    r=fit_local_signal(n,EDGES)
    assert r.converged and r.q0 > 0 and r.signal_yield > 0 and r.z_local_asymptotic > 1


def test_deterministic_pseudodata():
    assert np.array_equal(toy_pseudodata(seed=12),toy_pseudodata(seed=12))
    assert not np.array_equal(toy_pseudodata(seed=12),toy_pseudodata(seed=13))


def test_sideband_blinding():
    m=sideband_mask(EDGES)
    assert not m[25]
    assert m[10]
    assert np.sum(~m)==10


def test_exponential_sampling_range():
    a=exponential_samples(1000,seed=5)
    assert np.min(a)>=100 and np.max(a)<=160
