from dataclasses import dataclass

import numpy as np
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D


@dataclass
class ElectricFieldData:
    electric_field: Array1D
    x: Array1D

    @classmethod
    def from_file(cls, file_path: str) -> "ElectricFieldData":
        data = np.loadtxt(file_path)
        x = data[:, 0]
        electric_field = data[:, 1]
        return cls(electric_field, x)

    def plot(
        self,
        ax: Axes,
    ):
        plot_electric_field(
            ax=ax,
            electric_field=self.electric_field,
            x=self.x * M_TO_NM,
        )


def plot_electric_field(
    ax: Axes,
    electric_field: Array1D,
    x: Array1D,
):
    ax.plot(x, electric_field, "k-")
    ax.set_xlabel("Position [nm]")
    ax.set_ylabel("Squared electrical field [a.u.]")
