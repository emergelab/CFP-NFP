# CFP-NFP
Connectome Fingerprint scripts based on the following paper:

**A Unique Brain Connectome Fingerprint Predates and Predicts Response to Antidepressants**  
Nemati S, Akiki TJ, Roscoe J, Ju Y, Averill CL, Fouda S, Dutta A, McKie S, Krystal JH, Deakin JFW, Averill LA, Abdallah CG.  
iScience. 2020;23(1):100800. doi: [10.1016/j.isci.2019.100800](https://doi.org/10.1016/j.isci.2019.100800)

## Usage

Download this repository and run the Example/Test_getPMs, which will get NRS (funaNRS) then get the PMs (funPM.m & funpPM.m) and plot the connectome fingerprint (CFP).

To get the Akiki-Abdallah (AA) atlas and plot nodal results on the brain, use maps from: https://github.com/emergelab/hierarchical-brain-networks

## Included Fingerprints

FPs/AD_CFPs.mat (antidepressant CFP) and FPs/AD_NFP (nodal fingerprint) from Nemati, Akiki, et al. 2020, iScience . doi: [10.1016/j.isci.2019.100800](https://doi.org/10.1016/j.isci.2019.100800)

## Additional MATLAB functions

- `getPMs/funPM2.m`, `getPMs/funpPM2.m`: same interface and outputs as `funPM` / `funpPM`. They differ only for
  `model = 'wcpm'`: `funPM` scores the held-out subjects of each fold with the plain sum of the selected
  features, while `funPM2` applies the same r^2 feature weights to the held-out subjects as to the training
  subjects. Both give very similar correlations between predicted and observed scores (the `r` these scripts
  report). Use `funPM2` when you need the predicted values themselves on the scale of the outcome, for example
  prediction error, R^2, or individual predictions. `model = 'cpm'` and `'rcpm'` behave identically in both.
- `getPMs/applyFP.m`: apply a fingerprint from `FPs/` (or one you built) to new subjects. Returns each
  subject's fingerprint score, its correlation with an outcome, and optional cross-validated predictions
  calibrated in the new data (the score itself has arbitrary units).

**Applying a fingerprint: inputs must match how it was built.** For the published `FPs/AD_CFP.mat` and
`FPs/AD_NFP.mat` that means: the A424 parcellation with HCP-style preprocessing (ICA-FIX and global signal
regression); the same feature definition and AA level (e.g. AA-24 NRS for `AA24_NRSv`); change scores (week 1
minus baseline), not single-session connectivity; and the outcome defined as in the original analysis (percent
improvement, with the treatment contrast when comparing arms). On differently defined features the score is still
computed but is not meaningful.

## Python

`python/` contains `cfpnfp`, a Python implementation of the same steps (network-restricted strength, predictive
models with permutation testing, and applying fingerprints), plus A424 timeseries extraction from fMRIPrep
output. No MATLAB licence is needed.

    pip install "git+https://github.com/emergelab/CFP-NFP#subdirectory=python"

`compat="matlab2020"` reproduces `funaNRS`, `funPM` and `funpPM` exactly, and `compat=["test_weights"]`
reproduces `funPM2`; `python/tests` checks this against the MATLAB functions in this repository under GNU
Octave. See `python/README.md`.
