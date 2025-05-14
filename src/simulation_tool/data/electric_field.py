from dataclasses import asdict, dataclass
from pathlib import Path

import polars as pl
from matplotlib.axes import Axes

from simulation_tool.constants import M_TO_NM
from simulation_tool.typing_ import Array1D


@dataclass
class ElectricFieldData:
    electric_field: Array1D
    x: Array1D

    @classmethod
    def from_file(cls, file_path: Path) -> "ElectricFieldData":
        if file_path.suffix == ".parquet":
            data = pl.read_parquet(source=file_path)
            return cls(
                electric_field=data[:, 0].to_numpy(),
                x=data[:, 1].to_numpy(),
            )
        data = pl.read_csv(
            file_path,
            separator=" ",
            has_header=False,
        )
        return cls(
            electric_field=data[:, 1].to_numpy(),
            x=data[:, 0].to_numpy(),
        )

    def to_parquet(
        self,
        save_dir: Path,
        dtype: pl.DataType = pl.Float32,
    ):
        save_location = save_dir / "E_of_x.parquet"
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
