import shutil
from dataclasses import dataclass
from pathlib import Path

from simulation_tool.data import AbsorptionCoefficientData, ElectricFieldData, UVVisData
from simulation_tool.exceptions import DeviceParametersIncompleteError, SimulationError
from simulation_tool.jV.data import JVData
from simulation_tool.simsalabim.run import clean_up_simulation_output, run_simulation
from simulation_tool.templates.simss import SimssOutputFiles, UserInterface
from simulation_tool.utils import Plotable, save_figure


@dataclass
class JVSimulationOutput:
    uvvis: UVVisData
    electric_field: ElectricFieldData
    absorption_coefficient: AbsorptionCoefficientData
    jv: JVData


def run_jV_simulation(
    session_path: Path,
    path_to_executable: Path,
    simss_ui: UserInterface,
    setup_file: Path,
    plot_output: bool,
    plot_dpi: int,
) -> JVSimulationOutput | SimulationError | DeviceParametersIncompleteError:
    result = run_simulation(
        path_to_executable=path_to_executable,
        cmd_pars=[{"par": "dev_par_file", "val": str(setup_file)}],
        session_path=session_path,
        simss_ui=simss_ui,
    )

    if isinstance(result, SimulationError):
        return result

    uvvis = UVVisData.from_files(
        uvvis_file=session_path / SimssOutputFiles.ABSORPTION_SPECTRUM,
        rt_file=session_path / SimssOutputFiles.REFLECTION_TRANSMISSION_SPECTRUM,
    )

    electric_field = ElectricFieldData.from_file(session_path / SimssOutputFiles.E_OF_X)

    absorption_coefficient = AbsorptionCoefficientData.from_file(
        session_path / SimssOutputFiles.ALPHA_OF_X,
    )

    jv = JVData.from_files(
        device_characteristics_file=session_path / simss_ui.scParsFile,
        jv_file=session_path / simss_ui.JVFile,
    )

    output = [uvvis, electric_field, absorption_coefficient, jv]

    if plot_output:
        _plot_output(
            session_path=session_path,
            plotables=output,
            plot_dpi=plot_dpi,
        )

    _preserve_output(
        session_path=session_path, simulation_output=output, simss_ui=simss_ui
    )

    clean_up_simulation_output(session_path=session_path, simss_ui=simss_ui)

    if isinstance(result, DeviceParametersIncompleteError):
        return result

    return JVSimulationOutput(
        uvvis=uvvis,
        electric_field=electric_field,
        absorption_coefficient=absorption_coefficient,
        jv=jv,
    )


def _plot_output(
    session_path: Path,
    plotables: list[Plotable],
    plot_dpi: int,
) -> None:
    for plotable in plotables:
        save_figure(
            plotable,
            save_path=session_path / f"{plotable.__class__.__name__}.png",
            dpi=plot_dpi,
        )


def _preserve_output(
    session_path: Path,
    simulation_output: tuple[
        UVVisData, ElectricFieldData, AbsorptionCoefficientData, JVData
    ],
    simss_ui: UserInterface,
):
    for output in simulation_output:
        output.to_parquet(
            save_dir=session_path,
        )

    if (path := session_path / simss_ui.scParsFile).exists():
        if path.suffix == ".txt":
            raise ValueError(
                f"File {path} already has the .txt extension. Cannot rename it again."
            )
        shutil.move(path, path.with_suffix(".txt"))
