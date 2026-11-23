"""Locate the repository's data files (atlas index, fingerprints, example data).

Inside a checkout of CFP-NFP the files next to this package are used directly. When the package is installed
on its own (pip install ...#subdirectory=python), they are downloaded once from the GitHub repository into a
cache folder (CFPNFP_CACHE, default ~/.cache/cfpnfp).
"""
from __future__ import annotations

import os
import urllib.request
from pathlib import Path

REPO_RAW = "https://raw.githubusercontent.com/emergelab/CFP-NFP/master/"
ATLAS_RAW = "https://raw.githubusercontent.com/emergelab/hierarchical-brain-atlas/main/"

# path inside the CFP-NFP repo -> where to fetch it from
FILES = {
    "getPMs/atlasindex.mat": REPO_RAW + "getPMs/atlasindex.mat",
    "FPs/AD_CFP.mat": REPO_RAW + "FPs/AD_CFP.mat",
    "FPs/AD_NFP.mat": REPO_RAW + "FPs/AD_NFP.mat",
    "Example/TS.mat": REPO_RAW + "Example/TS.mat",
    "atlas/A424_Labels_AA-AAc_main_maps.csv": ATLAS_RAW + "data/A424_Labels_AA-AAc_main_maps.csv",
    "atlas/tpl-MNI152NLin2009cAsym_atlas-AA424Parcels2mm_res-02_desc-2mm_dseg.nii.gz":
        ATLAS_RAW + "bids-1.11.0/tpl-MNI152NLin2009cAsym/anat/"
                    "tpl-MNI152NLin2009cAsym_atlas-AA424Parcels2mm_res-02_desc-2mm_dseg.nii.gz",
}


def repo_root() -> Path | None:
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "getPMs" / "funaNRS.m").exists():
            return p
    return None


def path(name: str) -> Path:
    """Local path of a repository data file, downloading it to the cache if needed."""
    root = repo_root()
    if root is not None and (root / name).exists():
        return root / name
    cache = Path(os.environ.get("CFPNFP_CACHE", Path.home() / ".cache" / "cfpnfp"))
    dest = cache / name
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        urllib.request.urlretrieve(FILES[name], tmp)
        tmp.replace(dest)
    return dest
