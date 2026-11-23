"""The published fingerprints (FPs/AD_CFP.mat, FPs/AD_NFP.mat) and applying a fingerprint to new subjects.

Python counterpart of getPMs/applyFP.m.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from . import data


def node_labels() -> pd.DataFrame:
    """A424 node table (index 0..423): abbreviation, name, AA-7/24/50 network memberships."""
    t = pd.read_csv(data.path("atlas/A424_Labels_AA-AAc_main_maps.csv"), encoding="utf-8-sig")
    return pd.DataFrame({"node": t["AA-424-Abbreviation"], "name": t["AA-424-Name"],
                         "AA7": t["AA-7-Abbreviation"], "AA24": t["AA-24-Abbreviation"],
                         "AA50": t["AA-50-Abbreviation"]})


def load_published() -> dict:
    """{"cfp": {"pos"|"neg": {"AA24"|"AA50"|"AA424": m x m}},
        "nfp": {"pos"|"neg": DataFrame(424 x [node labels..., nS, neNRS, niNRS])}}"""
    import scipy.io as sio
    c = sio.loadmat(data.path("FPs/AD_CFP.mat"), squeeze_me=True, struct_as_record=False)["CFPs"]
    cfp = {s: {k.replace("_NRSv", ""): np.asarray(getattr(getattr(c, s), k), float)
               for k in getattr(c, s)._fieldnames} for s in ("pos", "neg")}
    out = {"cfp": cfp}
    try:
        import matio                                  # mat-io reads the MATLAB table objects in AD_NFP.mat
        n = matio.load_from_mat(str(data.path("FPs/AD_NFP.mat")))["NFPs"]
        lab = node_labels()
        out["nfp"] = {s: pd.concat([lab, n[s][0, 0].reset_index(drop=True)], axis=1) for s in ("pos", "neg")}
    except ImportError:
        out["nfp"] = None
    return out


@dataclass
class Applied:
    score: np.ndarray              # per subject, arbitrary units
    r: float = np.nan              # Pearson r of score with the outcome (scale-free)
    p: float = np.nan
    pred: np.ndarray | None = None  # cross-validated predictions on the outcome's scale


def apply_fingerprint(features, pos, neg, behav=None, k: int = 10, seed=None) -> Applied:
    """Apply a fingerprint to new subjects.

    features: subjects x features, or subjects x m x m (laid out like the fingerprint).
    pos/neg:  fingerprint weights with the shape of ONE subject's features.
    behav:    optional outcome; gives r/p and k-fold predictions calibrated in the new data (the score -> outcome
              line is fitted on training folds of the new data only).

    Inputs must match how the fingerprint was built. For the published FPs/AD_CFP.mat / AD_NFP.mat: A424
    parcellation with HCP-style preprocessing (ICA-FIX + global signal regression); the same feature definition and
    AA level (e.g. AA-24 NRS for "AA24"); CHANGE scores (week 1 minus baseline), not single-session connectivity;
    and the outcome defined as in the original analysis. Otherwise the score is computed but not meaningful.
    """
    X = np.asarray(features, float)
    n = X.shape[0]
    X = np.nan_to_num(X.reshape(n, -1))
    P = np.nan_to_num(np.asarray(pos, float).reshape(-1))
    N = np.nan_to_num(np.asarray(neg, float).reshape(-1))
    score = X @ P - X @ N
    res = Applied(score)
    if behav is None:
        return res
    y = np.asarray(behav, float).ravel()
    ok = ~np.isnan(y)
    res.r, res.p = (float(v) for v in stats.pearsonr(score[ok], y[ok]))
    idx = np.flatnonzero(ok)
    k = min(k, len(idx))
    folds = np.random.default_rng(seed).permutation(np.arange(len(idx)) % k)
    pred = np.full(n, np.nan)
    for f in range(k):
        te, tr = idx[folds == f], idx[folds != f]
        b = np.polyfit(score[tr], y[tr], 1)
        pred[te] = b[0] * score[te] + b[1]
    res.pred = pred
    return res
