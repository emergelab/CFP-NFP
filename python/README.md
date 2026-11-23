# cfpnfp

Python implementation of the CFP-NFP scripts in this repository: network-restricted strength (NRS; `funaNRS`),
connectome predictive models with permutation testing (`funPM`, `funpPM`, `funPM2`, `funpPM2`) and applying a
fingerprint to new subjects (`applyFP`).

    pip install "git+https://github.com/emergelab/CFP-NFP#subdirectory=python"           # core
    pip install "cfpnfp[fmri,nfp] @ git+https://github.com/emergelab/CFP-NFP#subdirectory=python"

## Use

```python
import numpy as np
from cfpnfp import atlas, nrs, pm, apply_fingerprint, load_published

M = nrs.connectome(ts)                         # ts: time x 424 (A424 regions)
ni, ne, NRS = nrs.nrs(M, atlas.level("AA24"))  # nodal internal / external strength, 24 x 24 network matrix
features = NRS[np.triu_indices(24)]            # one row per subject -> X (subjects x features)

res = pm.ppm(X, outcome, iterations=1000, model="wcpm", compat="matlab2020")   # as funpPM
res.r, res.p, res.pos, res.neg

fp = load_published()["cfp"]                   # FPs/AD_CFP.mat
a = apply_fingerprint(X_full, fp["pos"]["AA24"], fp["neg"]["AA24"], behav=outcome)   # X_full: subjects x 24 x 24
a.score, a.r, a.p, a.pred
```

Command line: `cfpnfp extract | nrs | pm | apply` (`cfpnfp -h`).

**Applying a fingerprint: inputs must match how it was built.** For the published `FPs/AD_CFP.mat` and
`FPs/AD_NFP.mat` that means: the A424 parcellation with HCP-style preprocessing (ICA-FIX and global signal
regression); the same feature definition and AA level (e.g. AA-24 NRS for `AA24_NRSv`); change scores (week 1
minus baseline), not single-session connectivity; and the outcome defined as in the original analysis (percent
improvement, with the treatment contrast when comparing arms). On differently defined features the score is still
computed but is not meaningful.

## Options

`compat` selects the behaviour:

| compat | matches |
|---|---|
| `"matlab2020"` (default) | `funaNRS`, `funPM`, `funpPM` |
| `["test_weights"]` | `funPM2`, `funpPM2` (held-out subjects scored with the training r^2 weights) |
| `["fold_covars"]` | covariates regressed out within each training fold and applied to its test fold |
| `["ties"]` | permutation p counts null r values equal to the observed mean r |
| `"v2"` | all three, plus `nrs(..., compat="fixed")`: internal strength averages over the other nodes of the network (k-1) |

Results are close across options; they matter most when you need predicted values on the outcome's scale or
use covariates.

## Performance

Implementation notes:
- `nrs`: network-restricted strengths come from one matrix product per level (node x module sums) instead of
  nested loops over module pairs; same numbers as `funaNRS` (see Tests).
- `pm`: in each fold, all feature-outcome correlations and their p-values are computed at once (vectorised
  Pearson r with the t-distribution p), not one `corr` call per feature.
- `ppm`: iterations are independent; `n_jobs` (CLI `--jobs`) runs them in parallel processes. Each iteration
  uses its own random generator spawned from `seed`, so results are identical for any `n_jobs`.

Measured on one core of a 2.3 GHz Xeon (synthetic data, `model="wcpm"`, 10 folds, true + permuted run per
iteration):

| task | time |
|---|---|
| connectome + NRS at three levels, one subject (420 x 424 timeseries) | ~5 ms |
| `ppm`, 202 subjects x 1,275 NRS features (AA-50) | ~50 ms per iteration (1,000 iterations < 1 min) |
| `ppm`, 202 subjects x 89,676 full-connectome features | ~4.6 s per iteration (1,000 iterations ~ 77 min; divide by cores with `n_jobs`) |

These are absolute timings for this implementation only; MATLAB was not benchmarked.

## Tests

`pytest python/tests` runs the MATLAB functions in `getPMs/` under GNU Octave (with the `statistics`
package) on `Example/TS.mat`, with recorded folds and permutations, and requires the Python results to match to
1e-9. Small shims in `tests/octave/shims` stand in for MATLAB toolbox functions (`corr` with p-values,
`crossvalind`, `randperm`, `nansum`, `nanmean`).

Container: `apptainer build cfpnfp.sif python/apptainer.def` (run from the `python/` folder).
