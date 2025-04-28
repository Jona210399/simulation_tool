from dataclasses import dataclass
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from simulation_tool.constants import M_TO_NM
from simulation_tool.templates.layer import Layer


@dataclass
class BandDiagramData:
    W_L: float
    W_R: float
    L_TCO: float
    L_BE: float
    Ls: list[float]
    E_cs: list[float]
    E_vs: list[float]

    @classmethod
    def from_layers(
        cls,
        W_L: float,
        W_R: float,
        L_TCO: float,
        L_BE: float,
        layers: list[Layer],
    ):
        Ls = [layer.L for layer in layers]
        E_cs = [layer.E_c for layer in layers]
        E_vs = [layer.E_v for layer in layers]
        return cls(W_L, W_R, L_TCO, L_BE, Ls, E_cs, E_vs)

    def plot(self, dpi: int, adjust_y_axis: bool, save_path: Path = Path(".")):
        plot_band_diagram(
            self,
            dpi=dpi,
            adjust_y_axis=adjust_y_axis,
            save_path=save_path,
        )


def plot_band_diagram(
    data: BandDiagramData,
    dpi: int,
    adjust_y_axis: bool,
    save_path: Path = Path("."),
):
    _, ax = plt.subplots(figsize=(10, 6))

    widths = np.array([data.L_TCO] + data.Ls + [data.L_BE]) * M_TO_NM
    E_cs = [data.W_L] + data.E_cs + [data.W_R]
    E_vs = [None] + data.E_vs + [None]
    names = [None] + [f"Layer {i + 1}" for i in range(len(data.Ls))] + [None]
    colors = ["C0"] + [f"C{i + 1}" for i in range(len(data.Ls))] + ["C0"]

    x_position = 0.0

    for width, E_c, E_v, name, color in zip(
        widths,
        E_cs,
        E_vs,
        names,
        colors,
    ):
        conduction_band = patches.Rectangle(
            (x_position, 0),
            width,
            E_c,
            linewidth=1,
            edgecolor="black",
            facecolor=color,
            alpha=0.3,
        )

        ax.add_patch(conduction_band)
        ax.plot(
            [x_position, x_position + width], [E_c, E_c], color=color, linestyle="--"
        )

        if E_v is None:
            x_position += width
            continue

        valence_band = patches.Rectangle(
            (x_position, E_v),
            width,
            50.0 - E_v,
            linewidth=1,
            edgecolor="black",
            facecolor=color,
            alpha=0.3,
        )
        ax.add_patch(valence_band)
        ax.plot(
            [x_position, x_position + width],
            [E_v, E_v],
            color=color,
            linestyle="--",
        )

        ax.text(
            x_position + width / 2,
            E_c + (E_v - E_c) / 2,
            name,
            ha="center",
            va="center",
            fontsize=10,
            rotation=90,
        )
        x_position += width

    y_min = min(E_cs) - 0.5 if adjust_y_axis else 1.0
    y_max = max(data.E_vs) + 0.5 if adjust_y_axis else 7.0

    plt.xlabel("Distance [nm]")
    plt.ylabel("Energy [eV]")
    plt.title("Solar Cell Band Diagram")
    plt.xlim(-0.1 * x_position, x_position + 0.1 * x_position)
    plt.ylim(y_min, y_max)
    plt.tight_layout()
    plt.savefig(f"{save_path}/band_diagram.png", dpi=dpi)
    plt.close()
