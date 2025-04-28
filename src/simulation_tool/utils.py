import numpy as np
from numpy.typing import NDArray
from scipy.constants import c, e, h

from simulation_tool.typing_ import Array1D


def gaussian(x: NDArray, mean: NDArray, sigma: NDArray, alpha: NDArray) -> NDArray:
    return alpha * np.exp(-((x - mean) ** 2) / (2 * sigma**2))


def convolute_gaussians(
    x: Array1D,
    means: Array1D,
    sigmas: Array1D,
    alphas: Array1D,
) -> Array1D:
    x = x[:, np.newaxis]  # shape (N, 1)
    means = means[np.newaxis, :]  # shape (1, G)
    sigmas = sigmas[np.newaxis, :]  # shape (1, G)
    alphas = alphas[np.newaxis, :]  # shape (1, G)

    gaussians = gaussian(x, means, sigmas, alphas)  # shape (N, G)
    return np.sum(gaussians, axis=1)  # shape (N,)


def get_num_points_for_linspace(start: float, stop: float, step: float) -> int:
    return int(round((stop - start) / step)) + 1


def wavelength_to_energy(wavelength: float) -> float:
    """Converts wavelength in meters to energy in eV."""
    energy = (h / e * c) / wavelength
    return energy


def loguniform(low=0, high=1, size=None):
    return np.exp(np.random.uniform(np.log(low), np.log(high), size))


def uniform(low=0, high=1):
    return np.random.uniform(low, high)


def randn():
    return np.random.randn()
