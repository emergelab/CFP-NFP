"""cfpnfp - Python implementation of the CFP-NFP connectome / nodal fingerprint method.

Based on Nemati S, Akiki TJ, et al. A Unique Brain Connectome Fingerprint Predates and Predicts Response to Antidepressants. iScience 2020;23(1):100800. doi:10.1016/j.isci.2019.100800

compat="matlab2020" reproduces getPMs/funaNRS.m, funPM.m and funpPM.m exactly (tests/test_parity.py).
compat="fixed" (or a list of fixes) gives the alternatives described in pm.py / nrs.py.
"""
from . import atlas, nrs, pm                           # noqa: F401
from .fingerprints import apply_fingerprint, load_published  # noqa: F401

__version__ = "0.1.0"
