from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes

from simulation_tool.constants import JV_SIM_OUTPUT_FILE_NAME
from simulation_tool.typing_ import Array1D, Array2D


@dataclass
class JVData:
    Vext: Array1D
    Jext: Array1D
    voc: float
    jsc: float
    ff: float

    @classmethod
    def from_files(
        cls,
        device_characteristics_file: Path,
        jv_file: Path,
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

        if jv_file.suffix == ".parquet":
            jv_data = pl.read_parquet(source=jv_file)
            vext = jv_data["Vext"].to_numpy()
            jext = jv_data["Jext"].to_numpy()
        else:
            jv_data = np.loadtxt(jv_file, skiprows=1, usecols=(0, 1))
            vext = jv_data[:, 0]
            jext = jv_data[:, 1]

        return cls(
            Vext=vext,
            Jext=jext,
            voc=voc,
            jsc=jsc,
            ff=ff,
        )

    def to_saveable_dict(self) -> dict[str, Array1D]:
        dict_ = asdict(self)
        keys_to_pop = ["voc", "jsc", "ff"]
        for key in keys_to_pop:
            dict_.pop(key)
        return dict_

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = (save_dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet")
        pl.DataFrame(self.to_saveable_dict()).with_columns(
            [pl.all().cast(dtype)]
        ).write_parquet(file=save_location)

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
            v_ext=self.Vext,
            j_ext=self.Jext,
            voc=self.voc,
            jsc=self.jsc,
            ff=self.ff,
            linewidth=linewidth,
            label=label,
            alpha=alpha,
        )

    @staticmethod
    def combine(
        jVs: list["JVData"],
    ) -> tuple[Array2D, Array2D, Array1D, Array1D, Array1D, Array1D]:
        v = np.stack([jv.Vext for jv in jVs])
        j = np.stack([jv.Jext for jv in jVs])
        ff = np.array([jv.ff for jv in jVs])
        voc = np.array([jv.voc for jv in jVs])
        jsc = np.array([jv.jsc for jv in jVs])
        pce = np.array([jv.calculate_pce() for jv in jVs])

        return v, j, ff, voc, jsc, pce


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

    if label:
        ax.legend(loc="upper left")
