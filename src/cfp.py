import numpy as np
from sklearn.model_selection import KFold
from sklearn.linear_model import Lasso
from scipy.stats import pearsonr


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
    data : _type_
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


def run_predictive_model(X, y, model, k, thresh, spearman, s, v_alpha, lambda_):
    """_summary_

    Parameters
    ----------
    X : np.array, shape (n_subjects, n_features)
        _description_
    y : _type_
        _description_
    model : _type_
        _description_
    k : int
        Number of k-fold cross validation splits.
    thresh : _type_
        _description_
    spearman : _type_
        _description_
    s : _type_
        _description_
    v_alpha : _type_
        _description_
    lambda_ : _type_
        _description_
    """
    # convert to 2D matrix
    # dim = X.shape
    # if len(dim) > 3 and dim[3] > 1:
    #     features = features.reshape(-1, dim[3])
    # elif len(dim) > 2 and dim[2] > 1:
    #     features = features.reshape(-1, dim[2])

    # remove NaN
    nanidx = np.isnan(y)
    y = y[~nanidx]
    features = X[:, ~nanidx]

    # get variables
    num_subj = features.shape[0]
    num_features = features.shape[1]
    if k > num_subj:
        k = num_subj

    # get cross validation indices
    kf = KFold(n_splits=k)

    pos_features = np.full((k, num_features), np.nan)
    neg_features = np.full((k, num_features), np.nan)
    r2w = np.zeros((k, num_features))
    coef_features = np.zeros((num_features, 1))
    behav_pred_both = np.full(num_subj, np.nan)

    for leftout, (train_index, test_index) in enumerate(kf.split(features)):
        # leave out subjects from features and behavior
        test_features = X[test_index]
        train_features = X[train_index]
        train_behav = y[train_index]

        # correlate all features with behavior
        r, p = np.array(
            [pearsonr(train_features[:, i], train_behav) for i in range(num_features)]
        ).T

        # set threshold and define masks
        pos_mask = (r > 0) & (p < thresh)
        neg_mask = (r < 0) & (p < thresh)
        mask = p < thresh

        # update features
        pos_features[leftout] = pos_mask
        neg_features[leftout] = neg_mask

        # get features' weight
        r[np.isnan(r)] = 0
        r2w[leftout] = r**2 * (pos_mask + neg_mask)

        # build and test
        if model in ["cpm", "wcpm"]:
            if model == "cpm":
                # get sum of all features
                train_sumpos = np.sum(train_features[:, pos_mask], axis=0)
                train_sumneg = np.sum(train_features[:, neg_mask], axis=0)
            elif model == "wcpm":
                # get sum of all features
                train_sumpos = np.sum(
                    train_features[:, pos_mask] * r2w[pos_mask, leftout][:, np.newaxis],
                    axis=0,
                )
                train_sumneg = np.sum(
                    train_features[:, neg_mask] * r2w[neg_mask, leftout][:, np.newaxis],
                    axis=0,
                )

            train_sumboth = train_sumpos - train_sumneg

            # build model on TRAIN subjs
            fit_both = np.polyfit(train_sumboth, train_behav, 1)

            # get sum of all features in TEST subjs
            test_sumpos = np.sum(test_features[pos_mask, :], axis=0)
            test_sumneg = np.sum(test_features[neg_mask, :], axis=0)
            test_sumboth = test_sumpos - test_sumneg

            # run model on TEST subjs
            behav_pred_both[test_index] = fit_both[0] * test_sumboth + fit_both[1]

        elif model == "rcpm":
            # build model on TRAIN subs
            if lambda_ is None:
                lasso = Lasso(alpha=v_alpha, fit_intercept=True, max_iter=10000)
                lasso.fit(train_features[mask, :].T, train_behav)
                coef = lasso.coef_
                coef0 = lasso.intercept_
