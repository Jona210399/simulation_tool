import shutil
from pathlib import Path

from pySIMsalabim.experiments.JV_steady_state import run_SS_JV

from simulation_tool.data import AbsorptionCoefficientData, ElectricFieldData, UVVisData
from simulation_tool.exceptions import SimulationError
from simulation_tool.jV.data import JVData


def run_jV_simulation(
    session_path: Path,
    simss_device_parameters: Path,
) -> None | SimulationError:
    result, message = run_SS_JV(
        str(simss_device_parameters),
        session_path,
        G_fracs=None,
    )

    if result.returncode != 0:
        return SimulationError(
            simulation_type="jV",
            return_value=result.returncode,
            message=message,
        )


def create_jV_simulation_plots(
    session_path: Path,
    dpi: int,
):
    JVData.from_files(
        device_characteristics_file=session_path / "scPars.txt",
        jv_file=session_path / "JV.dat",
    ).plot(
        save_path=session_path,
        dpi=dpi,
    )

    UVVisData.from_files(
        f"{session_path}/AbsorptionSpectrum.txt",
        f"{session_path}/reflection_transmission_spectrum.txt",
    ).plot(
        dpi=dpi,
        save_path=session_path,
    )

    ElectricFieldData.from_file(f"{session_path}/E_of_x.txt").plot(
        dpi=dpi,
        save_path=session_path,
    )

    AbsorptionCoefficientData.from_file(f"{session_path}/alpha_of_x.txt").plot(
        dpi=dpi,
        save_path=session_path,
    )


def preserve_jV_simulation_output(
    session_path: Path,
):
    shutil.move(
        session_path / "scPars.txt",
        session_path / "device_characteristics.txt",
    )
    shutil.move(
        session_path / "JV.dat",
        session_path / "jV.txt",
    )
