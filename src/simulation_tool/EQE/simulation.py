import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

from pySIMsalabim.experiments.EQE import run_EQE

from simulation_tool.EQE.data import EQEData
from simulation_tool.exceptions import SimulationError


@dataclass
class EQEParameters:
    spectrum: str
    lambda_min: float
    lambda_max: float
    lambda_step: float
    Vext: float


def run_EQE_simulation(
    session_path: Path,
    simss_device_parameters: Path,
    eqe_parameters: EQEParameters,
) -> None | SimulationError:
    return_value, message = run_EQE(
        simss_device_parameters=str(simss_device_parameters),
        session_path=session_path,
        run_mode=False,
        **asdict(eqe_parameters),
    )
    if return_value != 0:
        return SimulationError(
            simulation_type="EQE",
            return_value=return_value,
            message=message,
        )


def create_EQE_simulation_plots(
    session_path: Path,
    dpi: int,
):
    EQEData.from_file(session_path / "EQE.dat").plot(
        dpi=dpi,
        save_path=session_path,
    )


def preserve_EQE_simulation_output(
    session_path: Path,
):
    shutil.move(
        session_path / "EQE.dat",
        session_path / "EQE.txt",
    )
