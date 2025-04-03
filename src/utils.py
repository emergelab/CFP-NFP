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
