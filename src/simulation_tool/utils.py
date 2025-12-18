import shutil
import zipfile
from pathlib import Path
from typing import Protocol, TypeVar

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from numpy.typing import NDArray
from scipy.constants import c, e, h

from simulation_tool.typing_ import Array1D

T = TypeVar("T")


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


def wavelenght_to_energy(wavelength: float) -> float:
    """Converts wavelength in meters to energy in eV."""
    energy = (h / e * c) / wavelength
    return energy


def add_prefix_to_keys(
    dictionary: dict[str, T],
    prefix: str,
) -> dict[str, T]:
    """Adds a prefix to all keys in a dictionary."""
    return {f"{prefix}{key}": value for key, value in dictionary.items()}


def get_sorted_run_dirs(
    simulation_dir: Path,
    run_dir_prefix: str = "run",
) -> list[Path]:
    return sorted(
        (
            p
            for p in simulation_dir.iterdir()
            if p.is_dir() and p.name.startswith(f"{run_dir_prefix}_")
        ),
        key=lambda p: int(p.name.split("_")[-1]),
    )


class Plotable(Protocol):
    def plot(self, ax: Axes, **kwargs): ...


def save_figure(plotable: Plotable, save_path: Path, dpi: int, **kwargs) -> None:
    fig, ax = plt.subplots()
    plotable.plot(ax=ax, **kwargs)
    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi)
    plt.close(fig)


def zip_folder(path: Path, *, remove_origin: bool = False) -> tuple[Path, str]:
    """
    Create a ZIP archive from a folder.
    The archive is created next to the source folder and named '<folder>.zip'.
    """
    if not path.exists():
        raise ValueError(f"Path does not exist: {path}")
    if not path.is_dir():
        raise ValueError(f"Path is not a directory: {path}")

    zip_path = path.with_suffix(".zip")

    shutil.make_archive(
        base_name=str(zip_path.with_suffix("")),
        format="zip",
        root_dir=path.parent,
        base_dir=path.name,
    )

    if not zip_path.exists():
        return zip_path, f"{zip_path.name} archive was not created."

    with zipfile.ZipFile(zip_path, "r") as zf:
        if not zf.namelist():
            zip_path.unlink(missing_ok=True)
            return zip_path, f"Created {zip_path.name} archive is empty."

        bad_file = zf.testzip()
        if bad_file is not None:
            zip_path.unlink(missing_ok=True)
            return zip_path, f"Corrupt file in {zip_path.name} archive: {bad_file}"

    if remove_origin:
        shutil.rmtree(path)

    return zip_path, f"Successfully created {zip_path.name} archive."
