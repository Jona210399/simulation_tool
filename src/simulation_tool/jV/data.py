from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

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
        device_characteristics = pl.read_csv(
            device_characteristics_file,
            separator=" ",
        )

        jv_data = pl.read_csv(jv_file, separator=" ", truncate_ragged_lines=True)

        jsc, voc, ff = (
            device_characteristics["Jsc"].item(),
            device_characteristics["Voc"].item(),
            device_characteristics["FF"].item(),
        )

        return cls(
            jv_data["Vext"].to_numpy(),
            jv_data["Jext"].to_numpy(),
            voc,
            jsc,
            ff,
        )

    def plot(
        self,
        dpi: int,
        save_path: Path = Path("."),
    ):
        plot_jV(
            self.v_ext,
            self.j_ext,
            self.voc,
            self.jsc,
            self.ff,
            dpi,
            save_path,
        )


def plot_jV(
    v_ext: Array1D,
    j_ext: Array1D,
    voc: float,
    jsc: float,
    ff: float,
    dpi: int,
    save_path: Path = Path("."),
):
    plt.figure()
    plt.plot(
        v_ext,
        j_ext,
        label=f"Voc={voc:.2f}, Jsc={jsc:.1f}, FF={ff * 100.0:.1f}%",
    )
    plt.xlabel("V$_{ext}$ [V]")
    plt.ylabel("J$_{ext}$ [A m$^{-2}$]")
    plt.xlim([-0.5, 1.7])
    plt.ylim([-200, 100])
    plt.plot([-10, 10], [0, 0], "k--")
    plt.plot([0, 0], [-1000.0, 1000.0], "k--")
    plt.legend(loc="upper left")
    plt.xlabel("Voltage [V]")
    plt.ylabel("Current density [A/m^2]")
    plt.savefig(f"{save_path}/jV.png", dpi=dpi)
    plt.close()
