from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

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
        dpi: int,
        save_path: Path = Path("."),
    ):
        plot_absorption_coefficient(
            self.x * M_TO_NM,
            self.alpha * NM_TO_M,
            dpi,
            save_path,
        )


def plot_absorption_coefficient(
    x: Array1D,
    alpha: Array1D,
    dpi: int,
    save_path: Path = Path("."),
):
    plt.figure()
    plt.plot(x, alpha, "k-")
    plt.xlabel("Position [nm]")
    plt.ylabel("Absorption coefficient [1/nm]")
    plt.savefig(f"{save_path}/alpha_of_x.png", dpi=dpi)
    plt.close()
