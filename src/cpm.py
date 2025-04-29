import numpy as np
from sklearn.base import BaseEstimator, RegressorMixin
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from joblib import Parallel, delayed


class ConnectomePredictiveModel(BaseEstimator, RegressorMixin):
    """
    A class for a Connectome Predictive Model (CPM).

    This class is designed to handle the training and prediction of a CPM
    using a given connectome and behavioral data. It includes methods for
    training the model, making predictions, and evaluating the model's performance.
    The model is trained using a connectome matrix and a behavioral vector,
    and it can predict behavioral scores for new connectome data.

    Parameters
    ----------
    model : str, optional
        One of {"linear", "ridge", "lasso"}, by default "linear"
    weighted : bool, optional
        Whether to use weighted summation of features, by default True
    """

    def __init__(
        self,
        model: str = "linear",
        weighted: bool = True,
        p_threshold: float = 0.05,
        cross_validation_k: int = 10,
    ):
        if model not in ["linear", "ridge", "lasso"]:
            raise ValueError("Model must be one of 'linear', 'ridge', or 'lasso'.")
        else:
            self.model = model
        self.weighted = weighted
        self.p_threshold = p_threshold
        self.cross_validation_k = cross_validation_k

    def test(self, X: np.ndarray, y: np.ndarray, reps: int = 1000, worker: int = 1):
        """Test the model on the given data.

        Parameters
        ----------
        X : np.ndarray
            Shape (n_samples, n_regions, n_regions) or (n_samples, n_modalities, n_regions, n_regions)
            Allows for handling of multi-modality connectome data.
        y : np.ndarray
            Shape (n_samples, )
        reps : int, optional
            Number of repetitions for testing, by default 1000
        worker : int, optional
            Number of workers for parallel processing, by default 1

        Returns
        -------

        """
        n = X.shape[0]
        true_statistic = self._compute_statistic(X, y)

        null_distribution = Parallel(n_jobs=worker)(
            delayed(self._compute_statistic)(X, y[np.random.permutation(n)])
            for _ in range(reps)
        )
        null_distribution = np.array(null_distribution)
        p_value = (true_statistic < null_distribution).sum() / (reps + 1)
        return p_value

    def _compute_statistic(self, X: np.ndarray, y: np.ndarray):
        """_summary_

        Parameters
        ----------
        X : np.ndarray
            Shape (n_samples, n_regions, n_regions) or (n_samples, n_modalities, n_regions, n_regions)
            Allows for handling of multi-modality connectome data.
        y : np.ndarray
            Shape (n_samples, )
        reps : int, optional
            Number of repetitions for testing, by default 1000
        """
        if X.ndim not in [3, 4]:
            raise ValueError("X must be a 3D, or 4D array.")

        if self.model == "linear":
            clf = LinearRegression()
        elif self.model == "ridge":
            clf = Ridge()
        elif self.model == "lasso":
            clf = Lasso()

        kf = KFold(n_splits=self.cross_validation_k)
        y_pred = np.zeros(X.shape[0])

        for train_idx, test_idx in kf.split(X):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, _ = y[train_idx], y[test_idx]

            if X.ndim == 3:
                reshape_size = (-1, 1, 1)
            elif X.ndim == 4:
                reshape_size = (-1, 1, 1, 1)

            # compute and select features based on correlation p-vals
            corrs, pvals = pearsonr(X_train, y_train.reshape(reshape_size))
            pos_mask = (corrs > 0) & (pvals < self.p_threshold)
            neg_mask = (corrs < 0) & (pvals < self.p_threshold)

            # sum features based on mask
            # TODO: weighted
            axis = tuple(range(1, X_train.ndim))
            X_train_pos = (X_train * pos_mask).sum(axis=axis) / 2
            X_train_neg = (X_train * neg_mask).sum(axis=axis) / 2
            # X_train_both = X_train_pos - X_train_neg
            X_train_both = np.hstack(
                [X_train_pos.reshape(-1, 1), X_train_neg.reshape(-1, 1)]
            )

            # fit model
            clf.fit(X_train_both, y_train.reshape(-1, 1))

            X_test_pos = (X_test * pos_mask).sum(axis=axis) / 2
            X_test_neg = (X_test * neg_mask).sum(axis=axis) / 2
            # X_test_both = X_test_pos - X_test_neg
            X_test_both = np.hstack(
                [X_test_pos.reshape(-1, 1), X_test_neg.reshape(-1, 1)]
            )

            # predict
            y_pred[test_idx] = clf.predict(X_test_both).squeeze()

        behav_corr, _ = pearsonr(y, y_pred)

        return behav_corr
