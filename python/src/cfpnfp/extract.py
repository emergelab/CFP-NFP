"""Regional A424 timeseries from fMRIPrep volumetric output (MNI152NLin2009cAsym), for funaNRS-style input.

Denoising approximates the paper's FIX + mean grayordinate timeseries regression (MGTR), which cannot be redone
on these derivatives: 24 motion parameters (6 + derivatives + squares), white matter, global signal (+ derivatives),
linear detrend and a 0.01 Hz high-pass; volumes are z-scored per region. The recipe is recorded with every output.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from . import data

ATLAS_NAME = "atlas/tpl-MNI152NLin2009cAsym_atlas-AA424Parcels2mm_res-02_desc-2mm_dseg.nii.gz"
RECIPE = {"motion": "24P (6 + temporal derivatives + squares)", "tissue": ["WhiteMatter"],
          "global_signal": "GlobalSignal + derivative (MGTR analogue)", "detrend": True, "high_pass_hz": 0.01,
          "standardize": "zscore_sample", "atlas": os.path.basename(ATLAS_NAME)}


def confounds(tsv: str) -> tuple[np.ndarray, float]:
    df = pd.read_csv(tsv, sep="\t", na_values="n/a")
    mot = df[["X", "Y", "Z", "RotX", "RotY", "RotZ"]]
    der = mot.diff().fillna(0).add_suffix("_d")
    m24 = pd.concat([mot, der, (mot ** 2).add_suffix("_2"), (der ** 2).add_suffix("_d2")], axis=1)
    gs = df[["GlobalSignal"]]
    other = pd.concat([df[["WhiteMatter"]], gs, gs.diff().fillna(0).add_suffix("_d")], axis=1)
    X = pd.concat([m24, other], axis=1).fillna(0).to_numpy()
    fd = float(df["FramewiseDisplacement"].mean(skipna=True))
    return X, fd


def region_timeseries(bold: str, conf_tsv: str, atlas: str | None = None) -> tuple[np.ndarray, dict]:
    """Mean BOLD per A424 region, then confound regression + detrend + high-pass on the regional series
    (all linear, so identical to cleaning voxels first), then z-score. Regions outside the field of view -> NaN."""
    import nibabel as nb
    from nilearn import image, signal
    atlas = atlas or str(data.path(ATLAS_NAME))
    img = nb.load(bold)
    tr = float(img.header.get_zooms()[3])
    lab = np.asarray(image.resample_to_img(atlas, img, interpolation="nearest", force_resample=True,
                                           copy_header=True).dataobj).astype(int).ravel()
    data = np.asarray(img.dataobj, dtype=np.float32).reshape(-1, img.shape[3])
    counts = np.bincount(lab, minlength=425)[1:425]
    sums = np.zeros((424, data.shape[1]))
    np.add.at(sums, lab[lab > 0] - 1, data[lab > 0])
    with np.errstate(invalid="ignore", divide="ignore"):
        raw = (sums / counts[:, None]).T                      # time x 424
    X, fd = confounds(conf_tsv)
    ok = counts > 0
    ts = np.full(raw.shape, np.nan, dtype=np.float32)
    ts[:, ok] = signal.clean(raw[:, ok], confounds=X, detrend=True, high_pass=0.01, t_r=tr,
                             standardize="zscore_sample")
    meta = dict(RECIPE, tr=tr, volumes=int(ts.shape[0]), mean_fd=fd, empty_regions=int((~ok).sum()),
                small_regions=int(((counts > 0) & (counts < 5)).sum()))
    return ts, meta


def main(argv=None):
    import argparse, glob
    from concurrent.futures import ProcessPoolExecutor
    ap = argparse.ArgumentParser(description="A424 timeseries from fMRIPrep rest runs")
    ap.add_argument("fmriprep_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    jobs = []
    for bold in sorted(glob.glob(os.path.join(a.fmriprep_dir, "sub-*", "func",
                                              "*task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz"))):
        sub = os.path.basename(bold).split("_")[0]
        out = os.path.join(a.out_dir, f"{sub}_task-rest_atlas-A424_timeseries.npz")
        if not os.path.exists(out):
            jobs.append((bold, bold.replace("_space-MNI152NLin2009cAsym_preproc.nii.gz", "_confounds.tsv"), out))
    ok = 0
    with ProcessPoolExecutor(a.workers) as ex:
        for out, res in zip([j[2] for j in jobs], ex.map(_one, jobs)):
            if isinstance(res, str):
                print("failed", os.path.basename(out), res)
            else:
                ok += 1
    print(json.dumps({"done": ok, "attempted": len(jobs)}))


def _one(job):
    bold, conf, out = job
    try:
        ts, meta = region_timeseries(bold, conf)
        np.savez_compressed(out, ts=ts, meta=json.dumps(meta))
        return True
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"


if __name__ == "__main__":
    main()
