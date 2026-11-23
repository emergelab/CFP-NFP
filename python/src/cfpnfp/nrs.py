"""Connectome and network-restricted strength (NRS) -- Python port of getPMs/funaNRS.m.

compat="matlab2020" (default) reproduces funaNRS. compat="v2":
  * internal strength averages over the k-1 other nodes of a network (the zero self-connection is not counted);
  * optional z-scoring uses the off-diagonal edges only.
"""
from __future__ import annotations

import numpy as np

COMPAT = ("matlab2020", "v2", "fixed")


def connectome(ts: np.ndarray, is_corr: bool = False, zscore: bool = False, compat: str = "matlab2020") -> np.ndarray:
    """ts: time x nodes (or a nodes x nodes correlation matrix if is_corr). Returns Fisher-z matrix, diagonal 0."""
    if is_corr:
        M = np.array(ts, dtype=float, copy=True)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            M = np.corrcoef(np.asarray(ts, dtype=float), rowvar=False)
    np.fill_diagonal(M, 0.0)
    M = np.clip(M, -0.99, 0.99)
    M[~np.isfinite(M)] = 0.0
    M = np.arctanh(M)
    if zscore:
        if compat in ("v2", "fixed"):
            off = M[~np.eye(len(M), dtype=bool)]
            M = (M - off.mean()) / off.std(ddof=1)
            np.fill_diagonal(M, 0.0)
        else:
            M = (M - M.mean()) / M.std(ddof=1)
    return M


def nrs(M: np.ndarray, labels: np.ndarray, compat: str = "matlab2020"):
    """Network-restricted strength at one atlas level.

    M: n x n Fisher-z connectome (diagonal 0). labels: length-n module ids 1..m.
    Returns (niNRS[n], neNRS[n], NRS[m x m]); NRS[t, s] = mean edge weight between modules t and s,
    NRS[s, s] = internal strength of module s.
    """
    if compat not in COMPAT:
        raise ValueError(compat)
    labels = np.asarray(labels).astype(int)
    n, m = len(labels), int(labels.max())
    P = np.zeros((n, m))
    P[np.arange(n), labels - 1] = 1.0
    k = P.sum(0)                                  # module sizes
    S = M @ P                                     # node -> module sums (n x m)
    own = S[np.arange(n), labels - 1]
    kn = k[labels - 1]
    tot = S.sum(1)
    if m == n:                                    # every node its own module -> nodal strength (GBC)
        ni = tot / (n - 1 if compat in ("v2", "fixed") else n)
        return ni, np.zeros(n), M.copy()
    with np.errstate(invalid="ignore", divide="ignore"):
        if compat in ("v2", "fixed"):
            ni = own / (kn - 1)
            ni[kn == 1] = np.nan
        else:
            ni = own / kn
        ne = (tot - own) / (n - kn)
        B = P.T @ S                               # module x module sums
        D = np.outer(k, k)
        if compat in ("v2", "fixed"):
            np.fill_diagonal(D, k * (k - 1))
        N = B / D
    return ni, ne, N


def triu_colmajor(N: np.ndarray) -> np.ndarray:
    """Upper triangle (with diagonal) in MATLAB column-major order, i.e. NRS(logical(triu(ones(m))))."""
    m = len(N)
    idx = np.triu(np.ones((m, m), dtype=bool))
    return N.T[idx.T]                             # column-major walk of idx


def from_triu_colmajor(v: np.ndarray, m: int) -> np.ndarray:
    idx = np.triu(np.ones((m, m), dtype=bool))
    out = np.zeros((m, m))
    out.T[idx.T] = v
    return out


def subject_features(ts, levels_labels: dict, is_corr=False, zscore=False, compat="matlab2020"):
    """All NRS outputs for one subject. levels_labels: {name: labels}. Returns {name: (ni, ne, N)}."""
    M = connectome(ts, is_corr=is_corr, zscore=zscore, compat=compat)
    return {name: nrs(M, lab, compat=compat) for name, lab in levels_labels.items()}
