from dataclasses import dataclass

import numpy as np
from matplotlib.axes import Axes

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
        wavelengths = uvvis_data[:, 0]
        absorption = uvvis_data[:, 1]
        reflection = rt_data[:, 1]
        absorption = 1.0 - np.exp(-absorption)
        transmission = 1.0 - reflection - absorption

        return cls(wavelengths, absorption, reflection, transmission)

    def plot(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_uvvis(
            ax=ax,
            wavelengths=self.wavelengths * M_TO_NM,
            absorption=self.absorption * 100,
            reflection=self.reflection * 100,
            transmission=self.transmission * 100,
            linewidth=linewidth,
            alpha=alpha,
        )

    def plot_absorption(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_uvvis(
            ax=ax,
            wavelengths=self.wavelengths * M_TO_NM,
            absorption=self.absorption * 100,
            reflection=None,
            transmission=None,
            linewidth=linewidth,
            alpha=alpha,
        )

    def plot_reflection(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_uvvis(
            ax=ax,
            wavelengths=self.wavelengths * M_TO_NM,
            absorption=None,
            reflection=self.reflection * 100,
            transmission=None,
            linewidth=linewidth,
            alpha=alpha,
        )

    def plot_transmission(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_uvvis(
            ax=ax,
            wavelengths=self.wavelengths * M_TO_NM,
            absorption=None,
            reflection=None,
            transmission=self.transmission * 100,
            linewidth=linewidth,
            alpha=alpha,
        )


def plot_uvvis(
    ax: Axes,
    wavelengths: Array1D,
    absorption: Array1D | None,
    reflection: Array1D | None,
    transmission: Array1D | None,
    linewidth: float = 1.0,
    alpha: float = 1.0,
):
    ylabel = ""

    if absorption is not None:
        ax.plot(
            wavelengths,
            absorption,
            "-",
            label="absorption",
            linewidth=linewidth,
            alpha=alpha,
        )
        ylabel += "absorption"

    if reflection is not None:
        ax.plot(
            wavelengths,
            reflection,
            "-",
            label="reflection",
            linewidth=linewidth,
            alpha=alpha,
        )
        ylabel += "reflection" if ylabel == "" else " / reflection"

    if transmission is not None:
        ax.plot(
            wavelengths,
            transmission,
            "-",
            label="transmission",
            linewidth=linewidth,
            alpha=alpha,
        )
        ylabel += "transmission" if ylabel == "" else " / transmission"

    ylabel += " [%]"
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel(ylabel)

    if sum(x is not None for x in (absorption, reflection, transmission)) > 1:
        ax.legend(loc="best")
