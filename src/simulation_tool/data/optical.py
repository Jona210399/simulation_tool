from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes
from scipy.constants import nano as NANO

from simulation_tool.constants import M_TO_NM, NM_TO_M
from simulation_tool.typing_ import Array1D
from simulation_tool.utils import (
    convolute_gaussians,
    get_num_points_for_linspace,
    wavelength_to_energy,
)

MIN_WAVELENGTH = 350 * NANO  # 300 * NANO  # UNIT: m
MAX_WAVELENGTH = 900 * NANO  # 800 * NANO  # UNIT: m
WAVE_LENGTH_STEP = NANO  # UNIT: m


@dataclass
class OpticalData:
    wavelengths: Array1D
    n: Array1D
    k: Array1D
    bandgap: float = None

    @classmethod
    def from_file(cls, file_path: str) -> "OpticalData":
        data = pl.read_csv(file_path, separator=" ", truncate_ragged_lines=True)
        wavelengths = data["lambda"].to_numpy()
        n = data["n"].to_numpy()
        k = data["k"].to_numpy()
        return cls(wavelengths, n, k)

    def save(self, save_path: Path):
        save_location = save_path / "nk.txt"
        pl.DataFrame(self.to_saveable_dict()).write_csv(
            file=save_location, separator=" "
        )

    @classmethod
    def new(
        cls,
        min_wavelenght: float,
        max_wavelenght: float,
        wavelenght_step: float,
    ) -> "OpticalData":
        wavelengths = generate_wavelengths(
            min_wavelength=min_wavelenght,
            max_wavelength=max_wavelenght,
            wavelength_step=wavelenght_step,
        )
        n = generate_n(wavelengths)
        k, bandgap = generate_k(wavelengths)
        return cls(
            wavelengths,
            n,
            k,
            bandgap,
        )

    def plot(self, ax: Axes):
        plot_optical_data(
            ax=ax,
            **self.to_plotable_dict(),
        )

    def to_plotable_dict(self) -> dict[str, Array1D]:
        dict_ = asdict(self)
        dict_.pop("bandgap")
        return dict_

    def to_saveable_dict(self) -> dict[str, Array1D]:
        dict_ = self.to_plotable_dict()
        # preserves key order
        return {"lambda": dict_.pop("wavelengths"), **dict_}


def plot_optical_data(
    ax: Axes,
    wavelengths: Array1D,
    n: Array1D,
    k: Array1D,
):
    ax.plot(wavelengths * M_TO_NM, n, "k-", label="n")
    ax.plot(wavelengths * M_TO_NM, k, "b-", label="k")
    ax.legend()
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("n/k [a.u.]")


def generate_n(wavelengths: Array1D) -> Array1D:
    """
    Generates a refractive index n for the given wavelengths.
    Wavelengths are in meters.
    """
    num_gaussians = 10
    n_offset = 1.5
    means = np.random.random(num_gaussians) * 600 + 200
    sigmas = np.random.random(num_gaussians) * 60 + 10
    alphas = np.random.random(num_gaussians) * 0.5

    n = convolute_gaussians(wavelengths * M_TO_NM, means, sigmas, alphas)
    n += n_offset

    return n


def generate_k(wavelengths: Array1D) -> tuple[Array1D, float]:
    """
    Generates a refractive index k for the given wavelengths.
    Wavelengths are in meters.
    Returns the bandgap in eV as well
    """
    num_gaussians = 10
    bandgap = np.random.random() * 200 + 550
    means = np.random.random(num_gaussians) * (bandgap - 200) + 200
    sigmas = np.random.random(num_gaussians) * 60 + 10
    alphas = np.random.random(num_gaussians) * 0.2

    k = convolute_gaussians(wavelengths * M_TO_NM, means, sigmas, alphas)
    bandgap = np.max(means)
    bandgap = wavelength_to_energy(bandgap * NM_TO_M)
    return k, bandgap


def generate_wavelengths(
    min_wavelength: float,
    max_wavelength: float,
    wavelength_step: float,
) -> Array1D:
    """
    Returns a 1D array of wavelengths from MIN_WAVELENGTH to MAX_WAVELENGTH.
    Unit: m
    """
    wavelengths = np.linspace(
        min_wavelength,
        max_wavelength,
        get_num_points_for_linspace(min_wavelength, max_wavelength, wavelength_step),
    )
    return wavelengths
