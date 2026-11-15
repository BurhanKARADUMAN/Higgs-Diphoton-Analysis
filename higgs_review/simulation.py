"""Synthetic pseudodata for CI and public demos. NOT ATLAS collision data."""
from __future__ import annotations
import numpy as np


def exponential_samples(n, *, low=100, high=160, slope=-0.029, seed=41):
    rng = np.random.default_rng(seed)
    z = rng.uniform(size=int(n))
    if abs(slope) < 1e-12:
        return low + (high - low) * z
    return low + np.log1p(z * np.expm1(slope * (high - low))) / slope


def toy_pseudodata(*, background=12000, signal=400, mass=125.1, sigma=1.7, seed=42):
    rng = np.random.default_rng(seed)
    # Sample from a falling exponential using an independent RNG stream.
    bkg = exponential_samples(background, seed=seed + 17)
    sig = rng.normal(mass, sigma, int(signal))
    sig = sig[(sig >= 100) & (sig <= 160)]
    return rng.permutation(np.concatenate((bkg, sig)))
