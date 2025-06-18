from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes
from scipy.constants import c, pi
from scipy.interpolate import interp1d
from scipy.signal import hilbert

from simulation_tool.constants import (
    M_TO_NM,
    N_OFFSET_RANGE,
    NK_FILE_NAME,
    NM_TO_M,
    USE_KK_CONSISTENT_GENERATION,
)
from simulation_tool.typing_ import Array1D, Array2D
from simulation_tool.utils import (
    convolute_gaussians,
    get_num_points_for_linspace,
    wavelenght_to_energy,
)


@dataclass
class OpticalData:
    wavelenghts: Array1D
    n: Array1D
    k: Array1D
    bandgap: float = None

    @classmethod
    def from_file(cls, file_path: Path) -> "OpticalData":
        if file_path.suffix == ".parquet":
            data = pl.read_parquet(source=file_path)
        else:
            data = pl.read_csv(file_path, separator=" ", truncate_ragged_lines=True)

        wavelenghts = data["lambda"].to_numpy()
        n = data["n"].to_numpy()
        k = data["k"].to_numpy()
        return cls(wavelenghts, n, k)

    def save(self, save_dir: Path):
        save_location = save_dir / NK_FILE_NAME
        pl.DataFrame(self.to_saveable_dict()).write_csv(
            file=save_location, separator=" "
        )

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = (save_dir / NK_FILE_NAME).with_suffix(".parquet")
        pl.DataFrame(self.to_saveable_dict()).with_columns(
            [pl.all().cast(dtype)]
        ).write_parquet(file=save_location)

    @classmethod
    def new(
        cls,
        min_wavelenght: float,
        max_wavelenght: float,
        wavelenght_step: float,
    ) -> "OpticalData":
        if USE_KK_CONSISTENT_GENERATION:
            return kramers_kronig_consistent_generation(
                min_wavelenght=min_wavelenght,
                max_wavelenght=max_wavelenght,
                wavelenght_step=wavelenght_step,
            )
        return standard_generation(
            min_wavelenght=min_wavelenght,
            max_wavelenght=max_wavelenght,
            wavelenght_step=wavelenght_step,
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
        return {"lambda": dict_.pop("wavelenghts"), **dict_}

    @staticmethod
    def combine(optical: list["OpticalData"]) -> tuple[Array2D, Array2D, Array2D]:
        n = np.stack([opt.n for opt in optical])
        k = np.stack([opt.k for opt in optical])
        wavelenghts = np.stack([opt.wavelenghts for opt in optical])
        return n, k, wavelenghts


def plot_optical_data(
    ax: Axes,
    wavelenghts: Array1D,
    n: Array1D,
    k: Array1D,
):
    ax.plot(wavelenghts * M_TO_NM, n, "k-", label="n")
    ax.plot(wavelenghts * M_TO_NM, k, "b-", label="k")
    ax.legend()
    ax.set_xlabel("wavelenght [nm]")
    ax.set_ylabel("n/k [a.u.]")


def standard_generation(
    min_wavelenght: float,
    max_wavelenght: float,
    wavelenght_step: float,
) -> OpticalData:
    wavelenghts = _generate_wavelenghts(
        min_wavelenght=min_wavelenght,
        max_wavelenght=max_wavelenght,
        wavelenght_step=wavelenght_step,
    )
    n = _generate_n(wavelenghts)
    k, bandgap = _generate_k(wavelenghts)
    return OpticalData(
        wavelenghts=wavelenghts,
        n=n,
        k=k,
        bandgap=bandgap,
    )


def _generate_n(wavelenghts: Array1D) -> Array1D:
    """
    Generates a refractive index n for the given wavelenghts.
    wavelenghts are in meters.
    """
    num_gaussians = 10
    n_offset = 1.5
    means = np.random.random(num_gaussians) * 600 + 200
    sigmas = np.random.random(num_gaussians) * 60 + 10
    alphas = np.random.random(num_gaussians) * 0.5
    n_offset = np.random.uniform(*N_OFFSET_RANGE)

    n = convolute_gaussians(wavelenghts * M_TO_NM, means, sigmas, alphas)
    n += n_offset

    return n


def _generate_k(wavelenghts: Array1D) -> tuple[Array1D, float]:
    """
    Generates a refractive index k for the given wavelenghts.
    wavelenghts are in meters.
    Returns the bandgap in eV as well
    """
    num_gaussians = 10
    bandgap = np.random.random() * 200 + 550
    means = np.random.random(num_gaussians) * (bandgap - 200) + 200
    sigmas = np.random.random(num_gaussians) * 60 + 10
    alphas = np.random.random(num_gaussians) * 0.2

    k = convolute_gaussians(wavelenghts * M_TO_NM, means, sigmas, alphas)
    bandgap = np.max(means)
    bandgap = wavelenght_to_energy(bandgap * NM_TO_M)
    return k, bandgap


def _generate_wavelenghts(
    min_wavelenght: float,
    max_wavelenght: float,
    wavelenght_step: float,
) -> Array1D:
    """
    Returns a 1D array of wavelenghts from min_wavelenght to max_wavelenght.
    Unit: m
    """
    wavelenghts = np.linspace(
        min_wavelenght,
        max_wavelenght,
        get_num_points_for_linspace(min_wavelenght, max_wavelenght, wavelenght_step),
    )
    return wavelenghts


def kramers_kronig_consistent_generation(
    min_wavelenght: float,
    max_wavelenght: float,
    wavelenght_step: float,
) -> OpticalData:
    wavelenghts = _generate_wavelenghts(
        min_wavelenght=min_wavelenght,
        max_wavelenght=max_wavelenght,
        wavelenght_step=wavelenght_step,
    )
    k, bandgap = _generate_k(wavelenghts)

    n_offset = np.random.uniform(*N_OFFSET_RANGE)
    n = _kk_hilbert_n_from_k(wavelenghts, k, n_offset=n_offset)
    return OpticalData(
        wavelenghts=wavelenghts,
        n=n,
        k=k,
        bandgap=bandgap,
    )


def _wavelenghts_to_omega(wavelenghts: Array1D) -> Array1D:
    return 2 * pi * c / wavelenghts


def _kk_hilbert_n_from_k(
    wavelenghts: Array1D,
    k: Array1D,
    n_offset: float = 1.5,
) -> Array1D:
    """
    Compute n(λ) from k(λ) using the Kramers–Kronig relation via Hilbert transform.

    Assumptions:
    - wavelenghts is uniformly spaced in λ
    - k(λ) is defined on the same λ grid
    - Output n(λ) is consistent with the causality-based KK relation in frequency space
    """
    # Convert to angular frequency (ω), which is *non-uniform* even if λ is uniform
    omega = _wavelenghts_to_omega(wavelenghts)
    sort_idx = np.argsort(omega)
    omega_sorted = omega[sort_idx]
    k_sorted = k[sort_idx]

    # Interpolate k(ω) onto a *uniformly spaced ω grid*
    num_points = len(omega_sorted)
    omega_uniform = np.linspace(omega_sorted[0], omega_sorted[-1], num_points)
    k_interp = interp1d(omega_sorted, k_sorted, kind="cubic", fill_value="extrapolate")
    k_uniform = k_interp(omega_uniform)

    # Apply Hilbert transform in ω domain
    hilbert_k = np.imag(hilbert(k_uniform))  # Hilbert transform handles principal value
    n_uniform = 1 + (2 / pi) * hilbert_k + n_offset

    # Interpolate result back to original (unsorted) ω values
    omega_unsorted = _wavelenghts_to_omega(wavelenghts)
    n_interp = interp1d(
        omega_uniform, n_uniform, kind="cubic", fill_value="extrapolate"
    )
    n = n_interp(omega_unsorted)

    return n
