"""Akiki-Abdallah (AA / AAc) hierarchical network affiliations of the A424 nodes (getPMs/atlasindex.mat)."""
from __future__ import annotations

import numpy as np
import scipy.io as sio

from . import data

# column (1-based, as in Test_getPMs.m) -> number of modules, for the levels used in the paper / example
NAMED_LEVELS = {"AA7": 5, "AA24": 20, "AA50": 43, "AA424": 137}


def _index():
    return sio.loadmat(data.path("getPMs/atlasindex.mat"), squeeze_me=True,
                       struct_as_record=False)["atlasindex"].A424


def affiliations(modularity: str = "AAC") -> np.ndarray:
    """424 x L matrix of module ids (1..m) for every hierarchical level. "AAC" = whole brain (AA),
    "AAcC" = cortical (AAc)."""
    return np.asarray(getattr(_index(), modularity)).astype(int)


def level(name_or_column, modularity: str = "AAC") -> np.ndarray:
    """Module ids of the 424 nodes at one level: a name like "AA24" or a 1-based column number."""
    C = affiliations(modularity)
    col = NAMED_LEVELS[name_or_column] if isinstance(name_or_column, str) else int(name_or_column)
    return C[:, col - 1]


def level_for_modules(m: int, modularity: str = "AAC") -> int:
    """1-based column whose level has exactly m modules (e.g. 58 -> the AA-58 architecture)."""
    C = affiliations(modularity)
    hits = [i + 1 for i in range(C.shape[1]) if C[:, i].max() == m]
    if not hits:
        raise ValueError(f"no level with {m} modules")
    return hits[0]
