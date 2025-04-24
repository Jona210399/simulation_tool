from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D


@dataclass
class EQEData:
    wavelengths: Array1D
    eqe: Array1D

    @classmethod
    def from_file(cls, file_path: Path) -> "EQEData":
        data = pl.read_csv(
            file_path,
            separator=" ",
        )
        wavelengths = data["lambda"].to_numpy() * M_TO_NM
        eqe = data["EQE"].to_numpy()
        return cls(wavelengths, eqe)


def plot_EQE(
    eqe: Array1D,
    wavelengths: Array1D,
    dpi: int,
    save_path: Path = Path("."),
):
    plt.figure()
    plt.plot(wavelengths, eqe)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("EQE")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("EQE [%%]")
    plt.savefig(f"{save_path}/EQE.png", dpi=dpi)
    plt.close()
