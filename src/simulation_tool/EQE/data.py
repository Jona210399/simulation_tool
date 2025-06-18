from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D, Array2D


@dataclass
class EQEData:
    wavelenghts: Array1D
    EQE: Array1D

    @classmethod
    def from_file(cls, file_path: Path) -> "EQEData":
        if file_path.suffix == ".parquet":
            data = pl.read_parquet(source=file_path)
        else:
            data = pl.read_csv(
                file_path,
                separator=" ",
            )

        return cls.from_dataframe(data)

    @classmethod
    def from_dataframe(cls, df: pl.DataFrame) -> "EQEData":
        wavelengths = df["lambda"].to_numpy()
        eqe = df["EQE"].to_numpy()
        return cls(wavelengths, eqe)

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = save_dir / "EQE.parquet"
        pl.DataFrame(self.to_saveable_dict()).with_columns(
            [pl.all().cast(dtype)]
        ).write_parquet(file=save_location)

    def to_saveable_dict(self) -> dict[str, Array1D]:
        dict_ = asdict(self)
        # preserves key order
        return {"lambda": dict_.pop("wavelenghts"), **dict_}

    def plot(
        self,
        ax: Axes,
        linewidth: float = 1.0,
        alpha: float = 1.0,
    ):
        plot_EQE(
            ax=ax,
            eqe=self.EQE * 100,
            wavelengths=self.wavelenghts * M_TO_NM,
            linewidth=linewidth,
            alpha=alpha,
        )

    @staticmethod
    def combine(eqes: list["EQEData"]) -> tuple[Array2D, Array2D]:
        eqe = np.stack([eqe.EQE for eqe in eqes])
        wavelenghts = np.stack([eqe.wavelenghts for eqe in eqes])

        return eqe, wavelenghts


def plot_EQE(
    ax: Axes,
    eqe: Array1D,
    wavelengths: Array1D,
    linewidth: float = 1.0,
    alpha: float = 1.0,
):
    ax.plot(wavelengths, eqe, linewidth=linewidth, alpha=alpha)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel("EQE [%]")
