"""Massless diphoton kinematics (explicit GeV units)."""
from __future__ import annotations
import numpy as np


def diphoton_mass(pt1, eta1, phi1, pt2, eta2, phi2):
    """Return m_gg in GeV for massless photons with transverse momenta in GeV.

    m^2 = 2 pT1 pT2 [cosh(eta1-eta2)-cos(phi1-phi2)].
    Negative numerical roundoff is clipped to zero. Never infer input units.
    """
    pt1, pt2 = np.asarray(pt1, float), np.asarray(pt2, float)
    if np.any(~np.isfinite(pt1)) or np.any(~np.isfinite(pt2)) or np.any(pt1 < 0) or np.any(pt2 < 0):
        raise ValueError('pT must be finite, nonnegative, and expressed in GeV')
    deta = np.asarray(eta1, float) - np.asarray(eta2, float)
    dphi = np.asarray(phi1, float) - np.asarray(phi2, float)
    m2 = 2.0 * pt1 * pt2 * (np.cosh(deta) - np.cos(dphi))
    return np.sqrt(np.maximum(m2, 0.0))


def acceptance_mask(pt1, pt2, eta1, eta2, mass, *, eta_max=2.37, crack=(1.37, 1.52)):
    """An illustrative subset of ATLAS-style photon kinematic selections.

    Inputs are GeV; isolated / tight / trigger criteria are applied in the ROOT
    reader, not this function. pT ordering is the caller's responsibility.
    """
    p1, p2, e1, e2, m = np.broadcast_arrays(
        np.asarray(pt1, float), np.asarray(pt2, float),
        np.asarray(eta1, float), np.asarray(eta2, float), np.asarray(mass, float)
    )
    a1, a2 = np.abs(e1), np.abs(e2)
    good_eta = (a1 < eta_max) & (a2 < eta_max) & ~((a1 > crack[0]) & (a1 < crack[1])) & ~((a2 > crack[0]) & (a2 < crack[1]))
    return (np.isfinite(m) & (m >= 100) & (m <= 160) & (p1 > 25) & (p2 > 25) &
            (p1 / np.maximum(m, 1e-9) >= 0.35) & (p2 / np.maximum(m, 1e-9) >= 0.25) & good_eta)
