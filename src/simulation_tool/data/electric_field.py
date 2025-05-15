from dataclasses import asdict, dataclass
from pathlib import Path

import polars as pl
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM
from simulation_tool.templates.simss import SimssOutputFiles
from simulation_tool.typing_ import Array1D


@dataclass
class ElectricFieldData:
    x: Array1D
    electric_field: Array1D

    @classmethod
    def from_file(cls, file_path: Path) -> "ElectricFieldData":
        if file_path.suffix == ".parquet":
            data = pl.read_parquet(source=file_path)
        else:
            data = pl.read_csv(
                file_path,
                separator=" ",
                has_header=False,
            )

        return cls(
            x=data[:, 0].to_numpy(),
            electric_field=data[:, 1].to_numpy(),
        )

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = save_dir / SimssOutputFiles.E_OF_X.with_suffix(".parquet")
        pl.DataFrame(asdict(self)).with_columns([pl.all().cast(dtype)]).write_parquet(
            file=save_location
        )

    def plot(
        self,
        ax: Axes,
    ):
        plot_electric_field(
            ax=ax,
            electric_field=self.electric_field,
            x=self.x * M_TO_NM,
        )


def plot_electric_field(
    ax: Axes,
    electric_field: Array1D,
    x: Array1D,
):
    ax.plot(x, electric_field, "k-")
    ax.set_xlabel("Position [nm]")
    ax.set_ylabel("Squared electrical field [a.u.]")
