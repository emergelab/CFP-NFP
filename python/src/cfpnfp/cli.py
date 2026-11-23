"""Command line: cfpnfp {extract, nrs, pm, apply}.

  cfpnfp extract FMRIPREP_DIR OUT_DIR             A424 timeseries per subject (fMRIPrep volumetric, MNI2009cAsym)
  cfpnfp nrs TS_DIR OUT.npz [--levels AA7,AA24,AA50] [--compat matlab2020|v2] [--zscore]
                                                  NRS edge features (upper triangle incl. diagonal), niNRS, neNRS,
                                                  nodal strength; one row per subject (file order)
  cfpnfp pm FEATURES.npz KEY BEHAV.csv [--model wcpm] [--iterations 1000] [--compat matlab2020|v2] [--seed N] [--jobs N]
                                                  predictive model + permutation p (funpPM)
  cfpnfp apply FEATURES.npz KEY --fingerprint AA24|AA50|AA424 [--behav BEHAV.csv]
                                                  apply a published connectome fingerprint (FPs/AD_CFP.mat)
Timeseries files: .npy / .npz (key "ts") / .csv / .tsv, time x 424. BEHAV: one value per line, same order.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

from . import atlas, nrs, pm


def _load_ts(path):
    if path.endswith(".npz"):
        return np.load(path)["ts"]
    if path.endswith(".npy"):
        return np.load(path)
    return np.loadtxt(path, delimiter="\t" if path.endswith(".tsv") else ",")


def _behav(path):
    return np.array([float(x) if x.strip() not in ("", "nan", "NaN", "n/a") else np.nan
                     for x in open(path).read().split()])


def cmd_nrs(a):
    files = sorted(sum((glob.glob(os.path.join(a.ts_dir, p)) for p in ("*.npz", "*.npy", "*.csv", "*.tsv")), []))
    levels = a.levels.split(",")
    feats = {f"NRS_{L}": [] for L in levels} | {f"niNRS_{L}": [] for L in levels} | \
            {f"neNRS_{L}": [] for L in levels} | {"nS": []}
    for f in files:
        M = nrs.connectome(_load_ts(f), zscore=a.zscore, compat=a.compat)
        for L in levels:
            ni, ne, N = nrs.nrs(M, atlas.level(L), compat=a.compat)
            feats[f"NRS_{L}"].append(N[np.triu_indices(len(N))])
            feats[f"niNRS_{L}"].append(ni)
            feats[f"neNRS_{L}"].append(ne)
        feats["nS"].append(nrs.nrs(M, atlas.level("AA424"), compat=a.compat)[0])
    np.savez_compressed(a.out, subjects=np.array([os.path.basename(f) for f in files]),
                        **{k: np.array(v) for k, v in feats.items()})
    print(json.dumps({"subjects": len(files), "features": {k: len(v[0]) if v else 0 for k, v in feats.items()}}))


def cmd_pm(a):
    X = np.load(a.features)[a.key]
    r = pm.ppm(X, _behav(a.behav), iterations=a.iterations, model=a.model, compat=a.compat, seed=a.seed,
               n_jobs=a.jobs)
    print(json.dumps({"r": r.r, "p": r.p, "mse_ratio": r.mse_ratio, "n_pos": int((r.pos > 0).sum()),
                      "n_neg": int((r.neg > 0).sum())}))


def cmd_apply(a):
    from .fingerprints import apply_fingerprint, load_published
    fp = load_published()["cfp"]
    X = np.load(a.features)[a.key]
    m = fp["pos"][a.fingerprint].shape[0]
    if X.shape[1] != m * m:                           # upper-triangle vectors (cfpnfp nrs) -> full m x m layout
        full = np.zeros((X.shape[0], m, m))
        iu = np.triu_indices(m)
        full[:, iu[0], iu[1]] = X
        X = full
    res = apply_fingerprint(X, fp["pos"][a.fingerprint], fp["neg"][a.fingerprint],
                            _behav(a.behav) if a.behav else None, seed=a.seed)
    print(json.dumps({"score": res.score.tolist(), "r": res.r, "p": res.p,
                      "pred": None if res.pred is None else res.pred.tolist()}))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cfpnfp", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract"); e.add_argument("fmriprep_dir"); e.add_argument("out_dir")
    e.add_argument("--workers", type=int, default=8)
    n = sub.add_parser("nrs"); n.add_argument("ts_dir"); n.add_argument("out")
    n.add_argument("--levels", default="AA7,AA24,AA50"); n.add_argument("--compat", default="matlab2020")
    n.add_argument("--zscore", action="store_true")
    p = sub.add_parser("pm"); p.add_argument("features"); p.add_argument("key"); p.add_argument("behav")
    p.add_argument("--model", default="wcpm"); p.add_argument("--iterations", type=int, default=1000)
    p.add_argument("--compat", default="matlab2020"); p.add_argument("--seed", type=int)
    p.add_argument("--jobs", type=int, default=1, help="parallel processes for the iterations (-1 = all cores)")
    q = sub.add_parser("apply"); q.add_argument("features"); q.add_argument("key")
    q.add_argument("--fingerprint", required=True, choices=["AA24", "AA50", "AA424"])
    q.add_argument("--behav"); q.add_argument("--seed", type=int)
    a = ap.parse_args(argv)
    if a.cmd == "extract":
        from . import extract
        return extract.main([a.fmriprep_dir, a.out_dir, "--workers", str(a.workers)])
    return {"nrs": cmd_nrs, "pm": cmd_pm, "apply": cmd_apply}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
