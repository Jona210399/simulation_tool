from dataclasses import dataclass

import numpy as np
import polars as pl
from matplotlib.axes import Axes

from simulation_tool.typing_ import Array1D, PathLike


@dataclass
class JVData:
    v_ext: Array1D
    j_ext: Array1D
    voc: float
    jsc: float
    ff: float

    @classmethod
    def from_files(
        cls,
        device_characteristics_file: PathLike,
        jv_file: PathLike,
    ) -> "JVData":
        if device_characteristics_file.exists():
            device_characteristics = pl.read_csv(
                device_characteristics_file,
                separator=" ",
            )
            jsc, voc, ff = (
                device_characteristics["Jsc"].item(),
                device_characteristics["Voc"].item(),
                device_characteristics["FF"].item(),
            )

        else:
            voc = np.nan
            jsc = np.nan
            ff = np.nan

        jv_data = pl.read_csv(jv_file, separator=" ", truncate_ragged_lines=True)

        return cls(
            jv_data["Vext"].to_numpy(),
            jv_data["Jext"].to_numpy(),
            voc,
            jsc,
            ff,
        )

    def calculate_pce(self, pin: float = 1000.0) -> float:
        """
        Calculate the power conversion efficiency (PCE) in % based on the
        current density (Jsc), open-circuit voltage (Voc), fill factor (FF),
        and the incident light power density (Pin) in (W/m²).
        """
        return -1.0 * (self.jsc * self.voc * self.ff) / pin * 100.0

    def plot(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        label: bool = True,
        alpha: float = 1.0,
    ):
        plot_jV(
            ax=ax,
            v_ext=self.v_ext,
            j_ext=self.j_ext,
            voc=self.voc,
            jsc=self.jsc,
            ff=self.ff,
            linewidth=linewidth,
            label=label,
            alpha=alpha,
        )


def plot_jV(
    ax: Axes,
    v_ext: Array1D,
    j_ext: Array1D,
    voc: float,
    jsc: float,
    ff: float,
    linewidth: float = 1.0,
    label: bool = True,
    alpha: float = 1.0,
):
    label = f"Voc={voc:.2f}, Jsc={jsc:.1f}, FF={ff * 100.0:.1f}%" if label else None
    ax.plot(
        v_ext,
        j_ext,
        label=label,
        linewidth=linewidth,
        alpha=alpha,
    )
    ax.set_xlabel("V$_{ext}$ [V]")
    ax.set_ylabel("J$_{ext}$ [A m$^{-2}$]")

    ax.axhline(0, color="k", linestyle="--")
    ax.axvline(0, color="k", linestyle="--")
    ax.set_ylim([-200, 200])

    if label:
        ax.legend(loc="upper left")
