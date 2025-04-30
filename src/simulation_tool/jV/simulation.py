import shutil
from pathlib import Path

from pySIMsalabim.utils.general import run_simulation

from simulation_tool.data import AbsorptionCoefficientData, ElectricFieldData, UVVisData
from simulation_tool.exceptions import SimulationError
from simulation_tool.jV.data import JVData
from simulation_tool.utils import save_figure


def run_jV_simulation(
    session_path: Path,
    simss_device_parameters: Path,
) -> None | SimulationError:
    result, message = run_simulation(
        sim_type="simss",
        cmd_pars=[{"par": "dev_par_file", "val": str(simss_device_parameters)}],
        session_path=str(session_path),
        run_mode=True,
    )

    if result.returncode != 0:
        return SimulationError(
            simulation_type="jV",
            return_value=result.returncode,
            message=message,
        )

    files_to_check = [
        "AbsorptionSpectrum.txt",
        "reflection_transmission_spectrum.txt",
        "E_of_x.txt",
        "alpha_of_x.txt",
        "JV.dat",
        "scPars.dat",
    ]
    for file in files_to_check:
        if not (session_path / file).exists():
            return SimulationError(
                simulation_type="jV",
                return_value=1,
                message=f"Simulation did not produce the expected output files. Missing file: {file}",
            )


def create_jV_simulation_plots(
    session_path: Path,
    dpi: int,
):
    save_figure(
        UVVisData.from_files(
            uvvis_file=session_path / "AbsorptionSpectrum.txt",
            rt_file=session_path / "reflection_transmission_spectrum.txt",
        ),
        save_path=session_path / "uvvis.png",
        dpi=dpi,
    )

    save_figure(
        ElectricFieldData.from_file(f"{session_path}/E_of_x.txt"),
        save_path=session_path / "E_of_x.png",
        dpi=dpi,
    )

    save_figure(
        AbsorptionCoefficientData.from_file(f"{session_path}/alpha_of_x.txt"),
        save_path=session_path / "alpha_of_x.png",
        dpi=dpi,
    )

    save_figure(
        JVData.from_files(
            device_characteristics_file=session_path / "scPars.dat",
            jv_file=session_path / "JV.dat",
        ),
        save_path=session_path / "jV.png",
        dpi=dpi,
    )


def preserve_jV_simulation_output(
    session_path: Path,
):
    shutil.move(
        session_path / "scPars.dat",
        session_path / "device_characteristics.txt",
    )
    shutil.move(
        session_path / "JV.dat",
        session_path / "jV.txt",
    )
