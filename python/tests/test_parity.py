"""compat="matlab2020" must reproduce the repository's MATLAB functions (getPMs/funaNRS, funPM, funpPM) and
compat=["test_weights"] must reproduce getPMs/funPM2, on Example/TS.mat with identical folds/permutations.
The MATLAB side runs under GNU Octave (+ statistics package); the test is skipped if Octave is not installed."""
import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest
import scipy.io as sio

from cfpnfp import atlas, nrs as N, pm as PM, data

HERE = os.path.dirname(os.path.abspath(__file__))
LEVELS = [5, 20, 43, 137]
pytestmark = pytest.mark.skipif(shutil.which("octave") is None or data.repo_root() is None,
                                reason="needs GNU Octave and a CFP-NFP checkout")


@pytest.fixture(scope="module")
def ref(tmp_path_factory):
    work = tmp_path_factory.mktemp("parity")
    f = h5py.File(data.path("Example/TS.mat"), "r"); T = f["TS"]
    cases = ["".join(map(chr, f[r][()].ravel())) for r in T["cases"][()].ravel()]
    ts = [T[c][()].T for c in cases]
    behav = T["behav"][()].ravel()
    rng = np.random.default_rng(20201)
    n, k = len(cases), 10
    fold = lambda: rng.permutation(np.arange(n) % k) + 1          # noqa: E731
    fs, fi = np.stack([fold() for _ in range(6)]), np.stack([fold() for _ in range(10)])
    perms = np.stack([rng.permutation(n) for _ in range(5)])
    tsc = np.empty(n, dtype=object); tsc[:] = ts
    sio.savemat(work / "inputs.mat", {"ts": tsc, "behav": behav, "folds_single": fs.astype(float),
                                      "folds_iter": fi.astype(float), "perms": (perms + 1).astype(float)})
    repo = str(data.repo_root())
    cmd = f"repo='{repo}'; work='{work}'; here='{os.path.join(HERE, 'octave')}'; run('{os.path.join(HERE, 'octave', 'run_reference.m')}')"
    subprocess.run(["octave", "--no-gui", "--quiet", "--eval", cmd], check=True, capture_output=True, timeout=1800)
    r = sio.loadmat(work / "reference.mat", squeeze_me=True, struct_as_record=False)["ref"]
    return dict(ts=ts, behav=behav, fs=fs, fi=fi, perms=perms, ref=r)


def close(a, b, tol=1e-9):
    a, b = np.asarray(a, float), np.asarray(b, float)
    assert a.shape == b.shape
    d = np.abs(np.where(np.isnan(a) & np.isnan(b), 0, a - b))
    assert np.nanmax(d) <= tol, np.nanmax(d)


def features(ts, zscore=False):
    C = atlas.affiliations()
    out = []
    for x in ts:
        M = N.connectome(x, zscore=zscore)
        out.append({L: N.nrs(M, C[:, L - 1]) for L in LEVELS})
    return out


def nrs_matrix(feats, L, m):
    return np.array([N.from_triu_colmajor(N.triu_colmajor(f[L][2]), m).T.ravel() for f in feats])


def test_nrs(ref):
    for z in (0, 1):
        F = features(ref["ts"], bool(z))
        for li, L in enumerate(LEVELS):
            m = int(atlas.affiliations()[:, L - 1].max()); nt = min(m * (m + 1) // 2, 1275)
            for c, f in enumerate(F):
                close(f[L][0], getattr(ref["ref"], f"ni_z{z}")[:, li, c])
                close(f[L][1], getattr(ref["ref"], f"ne_z{z}")[:, li, c])
                close(N.triu_colmajor(f[L][2])[:nt], getattr(ref["ref"], f"tnrs_z{z}")[:nt, li, c])


@pytest.mark.parametrize("variant,compat", [("pm", "matlab2020"), ("pm2", ["test_weights"])])
def test_pm(ref, variant, compat):
    F = features(ref["ts"]); call = 0
    for L in LEVELS[:3]:
        m = int(atlas.affiliations()[:, L - 1].max()); X = nrs_matrix(F, L, m)
        for model in ("wcpm", "cpm"):
            r = PM.pm(X, ref["behav"], model=model, folds=ref["fs"][call], compat=compat); call += 1
            R = getattr(ref["ref"], f"{variant}_{model}_L{L}")
            close(r.r, R.r); close(r.pos.reshape(m, m, order="F"), R.pos); close(r.neg.reshape(m, m, order="F"), R.neg)


def test_ppm(ref):
    F = features(ref["ts"]); X = nrs_matrix(F, 20, 24); fi = ref["fi"]
    r = PM.ppm(X, ref["behav"], iterations=5, model="wcpm", d=0, compat="matlab2020",
               fold_sets=[(fi[2 * i], fi[2 * i + 1]) for i in range(5)], perm_sets=ref["perms"])
    R = ref["ref"].ppm_wcpm_L20
    close(r.r, R.r); close(r.p, R.p, 1e-12)
    close(r.pos.reshape(24, 24, order="F"), R.pos); close(r.neg.reshape(24, 24, order="F"), R.neg)


def test_apply_fingerprint_shapes():
    from cfpnfp import apply_fingerprint, load_published
    fp = load_published()["cfp"]
    X = np.random.default_rng(0).standard_normal((12, 24, 24))
    a = apply_fingerprint(X, fp["pos"]["AA24"], fp["neg"]["AA24"])
    b = apply_fingerprint(X.reshape(12, -1), fp["pos"]["AA24"].ravel(), fp["neg"]["AA24"].ravel())
    close(a.score, b.score)
