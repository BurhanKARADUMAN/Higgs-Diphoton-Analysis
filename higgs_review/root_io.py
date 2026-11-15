"""Opt-in ATLAS 2020 ROOT reader. No ROOT data are bundled or auto-downloaded.

Requires pip install -e '.[root]'. Photon pT and isolation units are explicit.
Some datasets / ROOT releases differ: verify TTree and branch metadata first.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from .kinematics import diphoton_mass, acceptance_mask

REQUIRED = ('photon_pt', 'photon_eta', 'photon_phi', 'photon_isTightID', 'photon_etcone20', 'trigP')


def iter_diphoton_masses(files, *, tree='mini', units='MeV', step_size='50 MB'):
    """Yield selected m_gg arrays from actual ATLAS-style open data ROOT files.

    Trigger, TightID, etcone20 isolation <= 4 GeV, fiducial and relative-pT
    cuts match a documented subset of the original ROOT macro. Real data only.
    """
    if units not in ('MeV', 'GeV'):
        raise ValueError('Explicit units must be MeV or GeV')
    try:
        import uproot
        import awkward as ak
    except ImportError as exc:
        raise RuntimeError("Install ROOT input dependencies with pip install -e '.[root]'") from exc
    scale = 1000.0 if units == 'MeV' else 1.0
    paths = [Path(p) for p in files]
    if not paths or any(not p.is_file() for p in paths):
        raise FileNotFoundError('Provide local ROOT files; no data is bundled')
    for path in paths:
        with uproot.open(path) as f:
            if tree not in f:
                raise KeyError(f'{path}: TTree {tree!r} not found')
            t = f[tree]
            missing = sorted(set(REQUIRED) - set(t.keys()))
            if missing:
                raise KeyError(f'{path}: required branches missing: {missing}')
            for arrays in t.iterate(REQUIRED, step_size=step_size, library='ak'):
                has2 = ak.num(arrays['photon_pt']) >= 2
                for key in ('photon_eta', 'photon_phi', 'photon_isTightID', 'photon_etcone20'):
                    has2 = has2 & (ak.num(arrays[key]) >= 2)
                if not bool(ak.any(has2)):
                    continue
                a = arrays[has2]
                # Explicitly sort the first two photons by pT for a stable convention.
                # The official analysis may use different photon ordering conventions.
                pt = ak.to_numpy(a['photon_pt'][:, :2]) / scale
                eta = ak.to_numpy(a['photon_eta'][:, :2])
                phi = ak.to_numpy(a['photon_phi'][:, :2])
                tight = ak.to_numpy(a['photon_isTightID'][:, :2]).astype(bool)
                iso = ak.to_numpy(a['photon_etcone20'][:, :2]) / scale
                trigger = ak.to_numpy(a['trigP']).astype(bool)
                lead = np.argmax(pt, axis=1)
                sub = 1 - lead
                idx = np.arange(len(pt))
                pt1, pt2 = pt[idx, lead], pt[idx, sub]
                eta1, eta2 = eta[idx, lead], eta[idx, sub]
                phi1, phi2 = phi[idx, lead], phi[idx, sub]
                mass = diphoton_mass(pt1, eta1, phi1, pt2, eta2, phi2)
                good = acceptance_mask(pt1, pt2, eta1, eta2, mass)
                good &= trigger & tight[:, 0] & tight[:, 1]
                good &= (iso[:, 0] < 4.0) & (iso[:, 1] < 4.0)
                yield np.asarray(mass[good], dtype=float)
