"""Connectome predictive models -- Python port of getPMs/funPM.m + funpPM.m.

Features are subjects x features (2-D). Models: "cpm" (sum), "wcpm" (r^2-weighted sum), "rcpm" (elastic net
with MATLAB lasso's Alpha/Lambda convention, standardised predictors).

compat="matlab2020" (default) reproduces funPM / funpPM:
  * wcpm: training subjects' summary scores are r^2-weighted sums; held-out subjects' are plain sums;
  * covariates (if given) are regressed out of the features across all subjects before cross-validation;
  * permutation p = (#null r > mean observed r + 1) / (iterations + 1).
Options (a list, or "v2" for all of them):
  "test_weights"  held-out subjects use the same r^2 weights (= funPM2 / funpPM2)
  "fold_covars"   covariates regressed out within each training fold and applied to its test fold
  "ties"          permutation p counts null r values equal to the observed mean (>=)
"""
from __future__ import annotations

from dataclasses import dataclass, field

FIXES = frozenset({"test_weights", "fold_covars", "ties"})   # self_edges lives in nrs.py (features)


def fixes_for(compat) -> frozenset:
    """compat: "matlab2020" (no fixes), "fixed" (all), or an iterable of FIXES names (ablation)."""
    if compat == "matlab2020":
        return frozenset()
    if compat in ("v2", "fixed"):
        return FIXES
    f = frozenset(compat) & FIXES
    return f

import numpy as np
from scipy import stats


def kfold_labels(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    """Random, (near) equal-sized folds labelled 1..k (what crossvalind('Kfold', n, k) produces)."""
    return rng.permutation(np.arange(n) % k) + 1


def _corr_p(X: np.ndarray, y: np.ndarray, spearman: bool = False):
    """Pearson (or Spearman) r and two-sided p of every column of X with y (MATLAB corr semantics)."""
    if spearman:
        X = stats.rankdata(X, axis=0)
        y = stats.rankdata(y)
    n = len(y)
    if np.isnan(X).any():                         # 'rows','complete' per column: rare, do it the slow way
        r = np.empty(X.shape[1]); p = np.empty(X.shape[1])
        for j in range(X.shape[1]):
            ok = ~np.isnan(X[:, j])
            r[j], p[j] = _corr_p(X[ok, j:j + 1], y[ok])
        return r, p
    Xc = X - X.mean(0)
    yc = y - y.mean()
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (Xc * yc[:, None]).sum(0) / np.sqrt((Xc ** 2).sum(0) * (yc ** 2).sum())
        r = np.clip(r, -1, 1)
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
    p[np.abs(r) == 1] = 0.0
    return r, p


def _residualise(train_X, train_C, test_X=None, test_C=None):
    A = np.column_stack([np.ones(len(train_C)), train_C])
    beta, *_ = np.linalg.lstsq(A, train_X, rcond=None)
    out = train_X - A @ beta
    if test_X is None:
        return out
    At = np.column_stack([np.ones(len(test_C)), test_C])
    return out, test_X - At @ beta


def _polyfit1(x, y):
    """MATLAB polyfit(x, y, 1) incl. its least-squares behaviour when x is constant (no features selected)."""
    A = np.column_stack([x, np.ones_like(x)])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return b


def _enet(X, y, alpha, lam):
    """MATLAB lasso(X, y, 'Alpha', alpha, 'Lambda', lam) with default standardisation -> (coef, intercept)."""
    from sklearn.linear_model import ElasticNet
    mu, sd = X.mean(0), X.std(0, ddof=0)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd
    m = ElasticNet(alpha=lam, l1_ratio=alpha, fit_intercept=True, max_iter=100000, tol=1e-6).fit(Z, y)
    coef = m.coef_ / sd
    return coef, m.intercept_ - mu @ coef


@dataclass
class PMResult:
    r: float
    p: float
    pos: np.ndarray            # feature weights, 0..1
    neg: np.ndarray
    coef: np.ndarray
    pred: np.ndarray | None = None
    mse_ratio: float | None = None        # mean squared prediction error / var(y): 1 = no better than the mean
    extra: dict = field(default_factory=dict)


def pm(features, behav, model="wcpm", k=10, thresh=0.05, spearman=False, s=1.0, v_alpha=1e-6, lam=1.0,
       folds=None, covars=None, compat="matlab2020", rng=None) -> PMResult:
    """One k-fold cross-validated predictive model (funPM.m). features: subjects x features."""
    rng = rng or np.random.default_rng()
    fx = fixes_for(compat)
    X = np.asarray(features, float)
    y = np.asarray(behav, float).ravel()
    keep = ~np.isnan(y)
    X, y = X[keep], y[keep]
    C = None if covars is None else np.asarray(covars, float).reshape(len(keep), -1)[keep]
    n, f = X.shape
    k = min(k, n)
    folds = kfold_labels(n, k, rng) if folds is None else np.asarray(folds)[keep]
    if C is not None and "fold_covars" not in fx:
        X = _residualise(X, C)
    posm = np.zeros((f, k)); negm = np.zeros((f, k)); r2w = np.zeros((f, k)); coefs = np.zeros((f, k))
    pred = np.full(n, np.nan)
    for i in range(k):
        te = folds == i + 1
        tr = ~te
        Xtr, Xte = X[tr], X[te]
        if C is not None and "fold_covars" in fx:
            Xtr, Xte = _residualise(Xtr, C[tr], Xte, C[te])
        r, p = _corr_p(Xtr, y[tr], spearman)
        pos = (r > 0) & (p < thresh)
        neg = (r < 0) & (p < thresh)
        mask = p < thresh
        posm[:, i], negm[:, i] = pos, neg
        r = np.nan_to_num(r)
        r2w[:, i] = r ** 2 * (pos | neg)
        if model in ("cpm", "wcpm"):
            w_tr = r2w[:, i] if model == "wcpm" else np.ones(f)
            w_te = w_tr if (model == "cpm" or "test_weights" in fx) else np.ones(f)
            s_tr = Xtr[:, pos] @ w_tr[pos] - Xtr[:, neg] @ w_tr[neg]
            s_te = Xte[:, pos] @ w_te[pos] - Xte[:, neg] @ w_te[neg]
            b = _polyfit1(s_tr, y[tr])
            pred[te] = b[0] * s_te + b[1]
        elif model == "rcpm":
            if mask.any():
                c, c0 = _enet(Xtr[:, mask], y[tr], v_alpha, lam)
                pred[te] = Xte[:, mask] @ c + c0
                coefs[mask, i] = c
            else:
                pred[te] = y[tr].mean()
        else:
            raise ValueError(model)
    mean_w = r2w.mean(1)
    pos_w = mean_w * (posm.sum(1) >= s * k)
    neg_w = mean_w * (negm.sum(1) >= s * k)
    coef = ((pos_w != 0) | (neg_w != 0)) * coefs.mean(1)
    if pos_w.max() > 0:
        pos_w = pos_w / pos_w.max()
    if neg_w.max() > 0:
        neg_w = neg_w / neg_w.max()
    ok = ~np.isnan(pred)
    if spearman:
        R, P = stats.spearmanr(pred[ok], y[ok])
    else:
        R, P = stats.pearsonr(pred[ok], y[ok]) if np.ptp(pred[ok]) > 0 else (np.nan, np.nan)
    return PMResult(float(R), float(P), pos_w, neg_w, coef, pred, float(np.mean((pred - y) ** 2) / np.var(y)))


def _iteration(job):
    """One repetition: true outcome, then a permuted outcome (as funpPM)."""
    X, y, args, seed, folds, perm = job
    rng = np.random.default_rng(seed)
    ft, fn = (None, None) if folds is None else folds
    a = pm(X, y, folds=ft, rng=rng, **args)
    perm = rng.permutation(len(y)) if perm is None else np.asarray(perm)
    b = pm(X, y[perm], folds=fn, rng=rng, **args)
    return a.r, b.r, a.mse_ratio, a.pos, a.neg, a.coef


def ppm(features, behav, iterations=1000, model="wcpm", k=10, thresh=0.05, spearman=False, s=1.0, d=0.0,
        v_alpha=1e-6, lam=1.0, covars=None, compat="matlab2020", seed=None, fold_sets=None,
        perm_sets=None, n_jobs=1) -> PMResult:
    """Repeated CV with a permutation null (funpPM.m).

    fold_sets / perm_sets (optional): per-iteration fold labels (true and null) and outcome permutations, so a
    run can be replayed exactly (e.g. against the MATLAB/Octave reference).
    n_jobs: iterations are independent and run in parallel processes when > 1 (-1 = all cores). Each iteration
    draws from its own generator spawned from `seed`, so results are identical for any n_jobs.
    """
    if iterations == 0:
        return pm(features, behav, model, k, thresh, spearman, s, v_alpha, lam,
                  folds=None if fold_sets is None else fold_sets[0][0], covars=covars, compat=compat,
                  rng=np.random.default_rng(seed))
    X = np.asarray(features, float)
    y = np.asarray(behav, float).ravel()
    n, f = X.shape
    seeds = np.random.SeedSequence(seed).spawn(iterations)
    args = dict(model=model, k=k, thresh=thresh, spearman=spearman, s=s, v_alpha=v_alpha, lam=lam,
                covars=covars, compat=compat)
    jobs = [(X, y, args, seeds[it], None if fold_sets is None else fold_sets[it],
             None if perm_sets is None else perm_sets[it]) for it in range(iterations)]
    if n_jobs == 1:
        out = [_iteration(j) for j in jobs]
    else:
        import os
        from concurrent.futures import ProcessPoolExecutor
        workers = os.cpu_count() if n_jobs in (-1, None) else n_jobs
        with ProcessPoolExecutor(workers) as ex:
            out = list(ex.map(_iteration, jobs, chunksize=max(1, iterations // (4 * workers))))
    Rt = np.array([o[0] for o in out]); Rn = np.array([o[1] for o in out]); mse = np.array([o[2] for o in out])
    P = np.column_stack([o[3] for o in out]); N = np.column_stack([o[4] for o in out])
    Cf = np.column_stack([o[5] for o in out])
    R = np.nanmean(Rt)
    Rn = np.nan_to_num(Rn)
    hits = (Rn >= R) if "ties" in fixes_for(compat) else (Rn > R)
    p = (hits.sum() + 1) / (iterations + 1)

    def agg(W):
        w = ((W > 0).sum(1) >= d * iterations) * W.sum(1)
        return w / w.max() if w.max() != 0 else w

    coef = Cf.sum(1)
    coef = coef / coef.max() if coef.max() != 0 else coef
    return PMResult(float(R), float(p), agg(P), agg(N), coef, None, float(np.nanmedian(mse)),
                    extra={"r_true": Rt, "r_null": Rn, "mse_ratio": mse})
