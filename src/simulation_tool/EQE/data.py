from dataclasses import dataclass
from pathlib import Path

import polars as pl
from matplotlib.axes import Axes

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

        return cls.from_dataframe(data)

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame) -> "EQEData":
        wavelengths = df["lambda"].to_numpy()
        eqe = df["EQE"].to_numpy()
        return cls(wavelengths, eqe)

    def plot(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_EQE(
            ax=ax,
            eqe=self.eqe * 100,
            wavelengths=self.wavelengths * M_TO_NM,
            linewidth=linewidth,
            alpha=alpha,
        )


def plot_EQE(
    ax: Axes,
    eqe: Array1D,
    wavelengths: Array1D,
    linewidth: float = 1.0,
    alpha: float = 1.0,
):
    ax.plot(wavelengths, eqe, linewidth=linewidth, alpha=alpha)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("EQE [%]")
