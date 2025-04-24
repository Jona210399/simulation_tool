from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D, PathLike


@dataclass
class UVVisData:
    wavelengths: Array1D
    absorption: Array1D
    reflection: Array1D
    transmission: Array1D

    @classmethod
    def from_files(
        cls,
        uvvis_file: PathLike,
        rt_file: PathLike,
    ) -> "UVVisData":
        uvvis_data = np.loadtxt(uvvis_file)
        rt_data = np.loadtxt(rt_file)
        wavelengths = uvvis_data[:, 0] * M_TO_NM
        absorption = uvvis_data[:, 1]
        reflection = rt_data[:, 1]
        absorption = 1.0 - np.exp(-absorption)
        transmission = 1.0 - reflection - absorption

        return cls(wavelengths, absorption, reflection, transmission)


def plot_uvvis(
    wavelengths: Array1D,
    absorption: Array1D,
    reflection: Array1D,
    transmission: Array1D,
    dpi: int,
    save_path: Path = Path("."),
):
    plt.figure()
    plt.plot(
        wavelengths,
        absorption,
        "-",
        color="C0",
        label="absorption",
    )
    plt.plot(
        wavelengths,
        reflection,
        "-",
        color="C1",
        label="reflection",
    )
    plt.plot(
        wavelengths,
        transmission,
        "-",
        color="C2",
        label="transmission",
    )
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Absorption [a.u.]")
    plt.legend(loc="best")
    plt.savefig(f"{save_path}/uvvis_absorption_spectrum.png", dpi=dpi)
    plt.close()
