from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from simulation_tool.constants import LINE_ALPHA, LINE_WIDTH, M_TO_NM
from simulation_tool.data import UVVisData
from simulation_tool.data.optical import OpticalData
from simulation_tool.EQE.data import EQEData
from simulation_tool.jV.data import JVData
from simulation_tool.typing_ import Array1D


def plot_device_characteristics_distributions(
    pces: Array1D,
    jscs: Array1D,
    vocs: Array1D,
    ffs: Array1D,
    save_path: Path,
    dpi: int,
):
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)

    plt.hist(pces, bins=20, color="C0", edgecolor="black")
    plt.xlabel("PCE [%]")
    plt.ylabel("Counts")
    plt.title("PCE distribution")

    plt.subplot(2, 2, 2)
    plt.hist(jscs, bins=20, color="C1", edgecolor="black")
    plt.xlabel("Jsc [A/m²]")
    plt.ylabel("Counts")
    plt.title("Jsc distribution")

    plt.subplot(2, 2, 3)
    plt.hist(vocs, bins=20, color="C2", edgecolor="black")
    plt.xlabel("Voc [V]")
    plt.ylabel("Counts")
    plt.title("Voc distribution")

    plt.subplot(2, 2, 4)
    plt.hist(ffs * 100, bins=20, color="C3", edgecolor="black")
    plt.xlabel("Fill Factor [%]")
    plt.ylabel("Counts")
    plt.title("FF distribution")

    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_jVs(
    jVs: list[JVData],
    save_path: Path,
    dpi: int,
    label: bool = False,
    ylim: tuple[float, float] | None = (-200, 200),
):
    _, ax = plt.subplots(figsize=(8, 6))

    for jV in jVs:
        jV.plot(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA, label=label)

    plt.ylim(*ylim) if ylim else None
    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_EQEs(eqes: list[EQEData], save_path: Path, dpi: int = 300):
    _, ax = plt.subplots(figsize=(8, 6))

    for eqe in eqes:
        eqe.plot(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_uvvis(
    uvvis: list[UVVisData],
    save_dir: Path,
    dpi: int,
):
    _, ax = plt.subplots(figsize=(8, 6))

    for uvv in uvvis:
        uvv.plot_absorption(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_dir / "absorption.png", dpi=dpi)
    plt.close()

    _, ax = plt.subplots(figsize=(8, 6))
    for uvv in uvvis:
        uvv.plot_reflection(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_dir / "reflection.png", dpi=dpi)
    plt.close()


def plot_optical_data(
    opticals: list[OpticalData],
    save_path: Path,
    dpi: int = 300,
):
    _, ax = plt.subplots(figsize=(8, 6))
    ax: Axes = ax

    for optical in opticals:
        n = ax.plot(
            optical.wavelenghts * M_TO_NM,
            optical.n,
            linewidth=LINE_WIDTH,
            alpha=LINE_ALPHA,
            linestyle="-",
        )
        color = n[0].get_color()

        ax.plot(
            optical.wavelenghts * M_TO_NM,
            optical.k,
            linewidth=LINE_WIDTH,
            alpha=LINE_ALPHA,
            color=color,
            linestyle="--",
        )

        ax.set_xlabel("wavelenght [nm]")
        ax.set_ylabel("n/k [a.u.]")

    legend_elements = [
        Line2D([0], [0], color="black", lw=LINE_WIDTH, linestyle="-", label="n"),
        Line2D([0], [0], color="black", lw=LINE_WIDTH, linestyle="--", label="k"),
    ]
    ax.legend(handles=legend_elements)

    plt.savefig(save_path, dpi=dpi)
    plt.close()
