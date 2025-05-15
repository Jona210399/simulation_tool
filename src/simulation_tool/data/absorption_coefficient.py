from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM, NM_TO_M
from simulation_tool.templates.simss import SimssOutputFiles
from simulation_tool.typing_ import Array1D


@dataclass
class AbsorptionCoefficientData:
    x: Array1D
    alpha: Array1D

    @classmethod
    def from_file(cls, file_path: Path) -> "AbsorptionCoefficientData":
        if file_path.suffix == ".parquet":
            data = pl.read_parquet(source=file_path)
        else:
            data = np.loadtxt(file_path)

        return cls(x=data[:, 0], alpha=data[:, 1])

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = save_dir / SimssOutputFiles.ALPHA_OF_X.with_suffix(".parquet")
        pl.DataFrame(asdict(self)).with_columns([pl.all().cast(dtype)]).write_parquet(
            file=save_location
        )

    def plot(
        self,
        ax: Axes,
    ):
        plot_absorption_coefficient(
            ax=ax,
            x=self.x * M_TO_NM,
            alpha=self.alpha * NM_TO_M,
        )


def plot_absorption_coefficient(
    ax: Axes,
    x: Array1D,
    alpha: Array1D,
):
    ax.plot(x, alpha, "k-")
    ax.set_xlabel("Position [nm]")
    ax.set_ylabel("Absorption coefficient [1/nm]")
