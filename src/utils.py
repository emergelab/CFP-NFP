import numpy as np
import h5py


def load_connectome_data(
    filename: str,
    group: str,
    dataset: str,
) -> np.ndarray:
    """Load connectome data from an HDF5 file.

    Parameters
    ----------
    filename : str
        Path to the HDF5 file.
    group : str
        Name of the group in the HDF5 file.
    dataset : str
        Name of the dataset within the group.

    Returns
    -------
    np.ndarray
        The loaded connectome data as a NumPy array.
    """
    with h5py.File(filename, "r") as f:
        data = f[group][dataset][:]

    return data


def compute_features(
    data,
    affiliations: np.ndarray,
    levels: int | None,
    is_corr: bool = False,
    zscore: bool = True,
) -> np.ndarray:
    """_summary_

    Parameters
    ----------
    data : np.array, shape (n_regions, n_features)
        time series or correlations matrix
    affiliations : _type_
        nodal affiliations; % whole-brain (AAC) or cortical (AAcC)
    levels : _type_
        [20, 43]; # [] to run all levels
    is_corr : _type_
        1; # if correlation matrix
    zscore : _type_
        1; # to normalize the connectome

    Returns
    -------
    _type_
        _description_
    """
    # TS = time series or correlations matrix
    # C = nodal affiliations; % whole-brain (AAC) or cortical (AAcC)
    # levels = [20, 43]; # [] to run all levels
    # dCor = 1; # if correlation matrix
    # zscore = 1; # to normalize the connectome

    if data.shape[0] != affiliations.shape[0]:
        raise ValueError("data and affiliations must have the same number of regions")

    n, _ = affiliations.shape

    if levels is None or len(levels) == 0:
        levels = list(range(_))

    n_levels = len(levels)

    # calculate Pearson networks
    if not is_corr:
        corrs = np.corrcoef(data, rowvar=True)
    else:
        corrs = data

    np.fill_diagonal(corrs, 0)  # zero diag
    corrs[corrs > 0.99] = 0.99
    corrs[corrs < -0.99] = -0.99  # clip extreme corr
    # remove Infs and NaNs
    ni = np.isnan(corrs) | np.isinf(corrs)
    if np.any(ni):
        corrs[ni] = 0
    corrs = np.arctanh(corrs)

    # for normalizing matrices
    if zscore == 1:
        mu = np.mean(corrs)
        sigma = np.std(corrs)
        corrs = (corrs - mu) / sigma

    niNRS = np.full((n, n_levels), np.nan)
    neNRS = np.full((n, n_levels), np.nan)  # total number of nodes
    tNRS = np.full((n * (n - 1) // 2 + n, n_levels), np.nan)  # total number of edges

    for idx, level in enumerate(levels):
        affiliation_level = affiliations[:, level - 1]  # matlab is 1-indexed
        level_max = affiliation_level.max()

        if level_max == n:  # FC; faster than the loop
            niNRS[:, idx] = np.mean(corrs, axis=1)  # i.e., nS
            neNRS[:, idx] = 0
            # idx = np.triu(np.ones((m, m)), k=1).astype(bool)
            tNRS[:, idx] = corrs[np.tril_indices_from(corrs)]
        else:
            NRS = np.full((level_max, level_max), np.nan)
            for s in range(1, level_max + 1):
                sindx = affiliation_level == s
                # nodal internal strength
                niNRS[sindx, idx] = np.mean(corrs[sindx][:, sindx], axis=1)
                # nodal external strength
                neNRS[sindx, idx] = np.mean(corrs[sindx][:, ~sindx], axis=1)
                # internetworks
                for t in range(1, level_max + 1):
                    tindx = affiliation_level == t
                    NRS[s - 1, t - 1] = np.mean(corrs[tindx][:, sindx])

            triu_idx = np.tril_indices_from(NRS)
            tNRS[: len(triu_idx[0]), idx] = NRS[triu_idx]

    return tNRS, niNRS, neNRS
