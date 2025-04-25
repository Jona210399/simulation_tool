from dataclasses import dataclass
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt

from simulation_tool.constants import M_TO_NM


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
    def from_settings(
        cls,
        settings: dict[str, str | int | float],
        num_layers: int = 3,
    ) -> "BandDiagramData":
        W_L = settings["W_L"]
        W_R = settings["W_R"]
        L_TCO = settings["L_TCO"] * M_TO_NM
        L_BE = settings["L_BE"] * M_TO_NM
        Ls = [settings[f"l{index + 1}.L"] * M_TO_NM for index in range(num_layers)]
        E_cs = [settings[f"l{index + 1}.E_c"] for index in range(num_layers)]
        E_vs = [settings[f"l{index + 1}.E_v"] for index in range(num_layers)]
        return cls(W_L, W_R, L_TCO, L_BE, Ls, E_cs, E_vs)

    def plot(self, dpi: int, save_path: Path = Path(".")):
        plot_band_diagram(self, dpi, save_path)


def plot_band_diagram(
    data: BandDiagramData,
    dpi: int,
    save_path: Path = Path("."),
):
    fig, ax = plt.subplots(figsize=(10, 6))
    # layer 0

    L_TCO = data.L_TCO
    L_BE = data.L_BE
    W_L = data.W_L
    W_R = data.W_R
    Ls = data.Ls
    E_cs = data.E_cs
    E_vs = data.E_vs

    x_position = 0.0
    width = L_TCO
    height = W_L
    color = "C0"
    rect = patches.Rectangle(
        (x_position, 0),
        width,
        height,
        linewidth=1,
        edgecolor="black",
        facecolor=color,
        alpha=0.3,
    )
    ax.add_patch(rect)
    ax.plot(
        [x_position, x_position + L_TCO], [height, height], color=color, linestyle="--"
    )
    x_position += L_TCO
    # layers 1, 2, and 3
    for i in range(3):
        L = Ls[i]
        E_c = E_cs[i]
        E_v = E_vs[i]
        name = f"Layer {i + 1}"
        width = L
        height1 = E_c
        height2 = 10.0 - E_v
        color = f"C{i + 1}"
        # Create a rectangle representing the band gap of the layer.
        rect1 = patches.Rectangle(
            (x_position, 0),
            width,
            height1,
            linewidth=1,
            edgecolor="black",
            facecolor=color,
            alpha=0.3,
        )
        rect2 = patches.Rectangle(
            (x_position, E_v),
            width,
            height2,
            linewidth=1,
            edgecolor="black",
            facecolor=color,
            alpha=0.3,
        )
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        # Plot dashed lines for the conduction and valence band edges.
        ax.plot([x_position, x_position + L], [E_c, E_c], color=color, linestyle="--")
        ax.plot([x_position, x_position + L], [E_v, E_v], color=color, linestyle="--")
        # Annotate the layer name in the middle of the box.
        ax.text(
            x_position + width / 2,
            E_c + height / 2,
            name,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
        # Update the x-position for the next layer.
        x_position += L
    # layer 4
    width = L_BE
    height = W_R
    color = "C0"
    rect = patches.Rectangle(
        (x_position, 0),
        width,
        height,
        linewidth=1,
        edgecolor="black",
        facecolor=color,
        alpha=0.3,
    )
    ax.add_patch(rect)
    ax.plot(
        [x_position, x_position + L_BE], [height, height], color=color, linestyle="--"
    )
    x_position += L_BE
    total_thickness = x_position

    # Set plot labels and title.
    ax.set_xlabel("Distance [nm]")
    ax.set_ylabel("Energy [eV]")
    ax.set_title("Solar Cell Band Diagram")
    ax.set_xlim(-0.1 * total_thickness, total_thickness + 0.1 * total_thickness)
    # Determine the overall energy range from all layers and contacts.
    # energy_min = min(min(E_cs), W_L, W_R) - 0.5
    # energy_max = max(max(E_vs), W_L, W_R) + 0.5
    ax.set_ylim(1.0, 7.0)

    plt.tight_layout()
    plt.savefig(f"{save_path}/band_diagram.png", dpi=dpi)
    plt.close()
