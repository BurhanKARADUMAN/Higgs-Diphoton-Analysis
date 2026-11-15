"""Exploratory *local* binned Poisson likelihood; not an ATLAS measurement.

The fixed Gaussian hypothesis, exponential background, no systematic nuisance,
and no trials correction are deliberately explicit. This is NOT a discovery tool.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, ndtr


@dataclass(frozen=True)
class FitResult:
    z_local_asymptotic: float
    q0: float
    signal_yield: float
    bg_yield: float
    log_background_rate_at_125: float
    background_slope_per_GeV: float
    converged: bool
    observed_counts: int
    mass_hypothesis_GeV: float
    resolution_sigma_GeV: float


def background_bin_counts(edges, log_rate, slope):
    """Expected exponential background counts; integrate each bin exactly."""
    e = np.asarray(edges, float)
    if e.ndim != 1 or len(e) < 2 or np.any(np.diff(e) <= 0):
        raise ValueError('Strictly increasing bin edges required')
    width = np.diff(e)
    mid = 0.5 * (e[1:] + e[:-1])
    x = slope * width / 2.0
    correction = np.ones_like(x)
    nz = np.abs(x) > 1e-8
    correction[nz] = np.sinh(x[nz]) / x[nz]
    return np.exp(log_rate + slope * (mid - 125.0)) * width * correction


def signal_bin_fractions(edges, mass=125.09, sigma=1.7):
    """Normal template fractions; normalized *within* histogram range."""
    e = np.asarray(edges, float)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('sigma must be positive')
    v = np.diff(ndtr((e - mass) / sigma))
    if np.sum(v) <= 0:
        raise ValueError('Signal mass outside histogram range')
    return v / v.sum()


def poisson_nll(counts, mu):
    n = np.asarray(counts, float)
    u = np.asarray(mu, float)
    if np.any(u <= 0) or np.any(~np.isfinite(u)):
        return np.inf
    return float(np.sum(u - n * np.log(u) + gammaln(n + 1.0)))


def fit_local_signal(counts, edges, *, mass=125.09, sigma=1.7):
    """Fit a fixed-mass hypothesis using a simplistic 1D local likelihood.

    Call only after predefining mass, model, binning and all selections. The
    asymptotic Z=sqrt(q0) is exploratory and *not* a calibrated local/global p.
    """
    n = np.asarray(counts, float)
    edges = np.asarray(edges, float)
    if n.ndim != 1 or len(n) != len(edges) - 1 or np.any(n < 0) or np.any(~np.isfinite(n)):
        raise ValueError('counts must be nonnegative, finite, and match edges')
    if n.sum() < 10:
        raise ValueError('Too few counts for this illustrative fit')
    fractions = signal_bin_fractions(edges, mass, sigma)
    width_avg = float(np.mean(np.diff(edges)))
    baseline = max(float(np.mean(n)) / width_avg, 1e-6)
    bounds_null = [(-20.0, 20.0), (-0.2, 0.2)]
    bounds_alt = bounds_null + [(0.0, max(1000.0, float(n.sum()) * 2.0))]
    def nullfn(p):
        return poisson_nll(n, background_bin_counts(edges, p[0], p[1]))
    def altfn(p):
        return poisson_nll(n, background_bin_counts(edges, p[0], p[1]) + p[2] * fractions)
    init = [np.log(baseline), -0.025]
    h0 = minimize(nullfn, init, method='L-BFGS-B', bounds=bounds_null)
    h1_runs = [minimize(altfn, init + [s], method='L-BFGS-B', bounds=bounds_alt)
               for s in (0.0, max(1.0, 0.03 * n.sum()), max(1.0, 0.1 * n.sum()))]
    h1 = min(h1_runs, key=lambda r: r.fun)
    if not h0.success or not h1.success:
        raise RuntimeError(f'Fit failed: null={h0.message}; alt={h1.message}')
    q0 = max(0.0, 2 * (h0.fun - h1.fun))
    b = background_bin_counts(edges, h1.x[0], h1.x[1])
    return FitResult(z_local_asymptotic=float(np.sqrt(q0)), q0=float(q0),
                     signal_yield=float(h1.x[2]), bg_yield=float(b.sum()),
                     log_background_rate_at_125=float(h1.x[0]),
                     background_slope_per_GeV=float(h1.x[1]), converged=True,
                     observed_counts=int(n.sum()), mass_hypothesis_GeV=float(mass),
                     resolution_sigma_GeV=float(sigma))


def sideband_mask(edges, blind=(120.0, 130.0)):
    e = np.asarray(edges, float)
    centers = (e[:-1] + e[1:]) / 2.0
    return (centers < blind[0]) | (centers >= blind[1])
