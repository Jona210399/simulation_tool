from dataclasses import dataclass

import numpy as np
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM, NM_TO_M
from simulation_tool.typing_ import Array1D


@dataclass
class AbsorptionCoefficientData:
    alpha: Array1D
    x: Array1D

    @classmethod
    def from_file(cls, file_path: str) -> "AbsorptionCoefficientData":
        data = np.loadtxt(file_path)
        x = data[:, 0]
        absorption_coefficient = data[:, 1]
        return cls(absorption_coefficient, x)

    def plot(
        self,
        ax: Axes,
    ):
        plot_absorption_coefficient(
            ax=ax,
            x=self.x * M_TO_NM,
            alpha=self.alpha * NM_TO_M,
        )


def plot_absorption_coefficient(
    ax: Axes,
    x: Array1D,
    alpha: Array1D,
):
    ax.plot(x, alpha, "k-")
    ax.set_xlabel("Position [nm]")
    ax.set_ylabel("Absorption coefficient [1/nm]")
