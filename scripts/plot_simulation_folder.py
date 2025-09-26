from pathlib import Path

from simulation_tool.data import AbsorptionCoefficientData, OpticalData, UVVisData
from simulation_tool.EQE.data import EQEData
from simulation_tool.jV.data import JVData
from simulation_tool.utils import save_figure


def plot_single_simulation_folder(simulation_folder: Path) -> None:
    data_folder = simulation_folder / "data"
    plot_folder = simulation_folder / "plots"
    plot_folder.mkdir(exist_ok=True)

    uvvis = UVVisData.from_file(data_folder / "uvvis.parquet")
    save_figure(uvvis, plot_folder / f"{uvvis.__class__.__name__}.png", dpi=300)

    optical = OpticalData.from_file(data_folder / "nk.parquet")
    save_figure(optical, plot_folder / f"{optical.__class__.__name__}.png", dpi=300)

    absorption = AbsorptionCoefficientData.from_file(data_folder / "alpha_of_x.parquet")
    save_figure(
        absorption, plot_folder / f"{absorption.__class__.__name__}.png", dpi=300
    )

    jv = JVData.from_files(
        device_characteristics_file=data_folder / "scPars.txt",
        jv_file=data_folder / "jV.parquet",
    )

    save_figure(jv, plot_folder / f"{jv.__class__.__name__}.png", dpi=300)

    eqe = EQEData.from_file(data_folder / "EQE.parquet")
    save_figure(eqe, plot_folder / f"{eqe.__class__.__name__}.png", dpi=300)

    return


if __name__ == "__main__":
    plot_single_simulation_folder(
        Path(
            "/p/project1/solai/oestreicher1/repos/simsalabim/simulation_tool/simulation_runs/697856/run_0"
        )
    )
