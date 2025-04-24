from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D


@dataclass
class ElectricFieldData:
    electric_field: Array1D
    x: Array1D

    @classmethod
    def from_file(cls, file_path: str) -> "ElectricFieldData":
        data = np.loadtxt(file_path)
        x = data[:, 0] * M_TO_NM
        electric_field = data[:, 1]
        return cls(electric_field, x)


def plot_electric_field(
    electric_field: Array1D,
    x: Array1D,
    dpi: int,
    save_path: Path = Path("."),
):
    plt.figure()
    plt.plot(x, electric_field, "k-")
    plt.xlabel("Position [nm]")
    plt.ylabel("Squared electrical field [a.u.]")
    plt.savefig(f"{save_path}/E_of_x.png", dpi=dpi)
    plt.close()
