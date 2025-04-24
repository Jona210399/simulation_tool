import os
import shutil
from dataclasses import asdict
from pathlib import Path

import joblib
import numpy as np
import pySIMsalabim as sim
from pySIMsalabim.experiments.EQE import run_EQE
from pySIMsalabim.experiments.JV_steady_state import run_SS_JV

import simulation_tool.templates as templates
import simulation_tool.templates.to_text as to_text
from simulation_tool.absorption_coefficient import (
    AbsorptionCoefficientData,
    plot_absorption_coefficient,
)
from simulation_tool.band_diagramm import BandDiagramData, plot_band_diagram
from simulation_tool.electric_field import ElectricFieldData, plot_electric_field
from simulation_tool.EQE import EQEData, plot_EQE
from simulation_tool.jV import JVData, plot_jV
from simulation_tool.optical import OpticalData
from simulation_tool.session import set_up_session
from simulation_tool.templates.random import sample_random_parameters
from simulation_tool.uvvis import UVVisData, plot_uvvis

PATH_TO_SIMSS = Path(
    "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/SimSS"
)
DPI = 60
G_FRAC = 1.0
SPECTRUM = "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/Data/AM15G.txt"

LAMBDA_MIN = 350
LAMBDA_MAX = 900
LAMBDA_STEP = 10
VEXT = 0


def create_settings(
    session_path: Path,
    path_to_simss: Path,
    bandgap: float,
    manual_randomization: bool,
) -> dict[str, str | int | float]:
    if manual_randomization:
        simss_config = templates.get_simss_config(session_path, path_to_simss)
        layer1_config = templates.get_layer1_config(path_to_simss)
        layer2_config = templates.get_layer2_config(session_path, bandgap)
        layer3_config = templates.get_layer3_config(path_to_simss)

        settings = {}

        configs = [simss_config, layer1_config, layer2_config, layer3_config]
        for config in configs:
            settings.update(config)
    else:
        settings = sample_random_parameters(
            session_path,
            path_to_simss,
            bandgap,
        )
    simss_config_text = to_text.get_simss_config(settings)
    layer1_config_text = to_text.get_layer_config(settings, "l1")
    layer2_config_text = to_text.get_layer_config(settings, "l2")
    layer3_config_text = to_text.get_layer_config(settings, "l3")

    with open(session_path / "simulation_setup_custom.txt", "w") as outfile:
        outfile.write(simss_config_text)

    with open(session_path / "L1_parameters.txt", "w") as outfile:
        outfile.write(layer1_config_text)

    with open(session_path / "L2_parameters.txt", "w") as outfile:
        outfile.write(layer2_config_text)

    with open(session_path / "L3_parameters.txt", "w") as outfile:
        outfile.write(layer3_config_text)

    return settings


def prepare_simulation(session_path: Path, manual_randomization: bool):
    optical_data = OpticalData.new()
    optical_data.plot(dpi=DPI, save_path=session_path)
    optical_data.save(session_path)

    settings = create_settings(
        session_path,
        PATH_TO_SIMSS,
        optical_data.bandgap,
        manual_randomization=manual_randomization,
    )

    plot_band_diagram(
        BandDiagramData.from_settings(settings),
        save_path=session_path,
        dpi=DPI,
    )


class SimulationError(Exception):
    def __init__(
        self,
        simulation_type: str,
        return_value: int,
        message: str,
    ):
        error_msg = f"{simulation_type} simulation failed with return value {return_value} and message: {message}"

        self.simulation_type = simulation_type
        self.return_value = return_value
        self.message = message

        super().__init__(error_msg)


def run_jV_simulation(
    session_path: Path,
    simss_device_parameters: Path,
) -> None | SimulationError:
    return_value, message = run_SS_JV(
        str(simss_device_parameters),
        session_path,
        JV_file_name="JV.dat",
        varFile="Var.dat",
        G_fracs=[G_FRAC],
    )

    if return_value != 0:
        return SimulationError(
            simulation_type="jV",
            return_value=return_value,
            message=message,
        )


def process_jV_simulation_output(session_path: Path) -> tuple[Path, Path]:
    shutil.move(
        f"{session_path}/log_Gfrac_{G_FRAC}.txt",
        f"{session_path}/log_jV.txt",
    )

    os.remove(f"{session_path}/Var_Gfrac_{G_FRAC}.dat")

    device_characteristics_file = shutil.move(
        f"{session_path}/scPars_Gfrac_{G_FRAC}.txt",
        f"{session_path}/device_characteristics.txt",
    )
    jv_file = shutil.move(
        f"{session_path}/JV_Gfrac_{G_FRAC}.dat",
        f"{session_path}/jV.txt",
    )

    return device_characteristics_file, jv_file


def plot_jV_simulation_output(
    session_path: Path,
    device_characteristics_file: Path,
    jv_file: Path,
):
    plot_jV(
        **asdict(
            JVData.from_files(
                device_characteristics_file=device_characteristics_file,
                jv_file=jv_file,
            )
        ),
        save_path=session_path,
        dpi=DPI,
    )

    plot_uvvis(
        **asdict(
            UVVisData.from_files(
                f"{session_path}/AbsorptionSpectrum.txt",
                f"{session_path}/reflection_transmission_spectrum.txt",
            )
        ),
        dpi=DPI,
        save_path=session_path,
    )

    plot_electric_field(
        **asdict(ElectricFieldData.from_file(f"{session_path}/E_of_x.txt")),
        dpi=DPI,
        save_path=session_path,
    )

    plot_absorption_coefficient(
        **asdict(AbsorptionCoefficientData.from_file(f"{session_path}/alpha_of_x.txt")),
        dpi=DPI,
        save_path=session_path,
    )


def run_EQE_simulation(
    session_path: Path,
    simss_device_parameters: Path,
) -> None | SimulationError:
    EQE_outfile_name = "EQE.dat"
    spectrum = SPECTRUM
    return_value, message = run_EQE(
        str(simss_device_parameters),
        session_path,
        spectrum,
        LAMBDA_MIN,
        LAMBDA_MAX,
        LAMBDA_STEP,
        VEXT,
        outfile_name=EQE_outfile_name,
        JV_file_name="JV.dat",
        run_mode=False,
        parallel=False,
    )
    if return_value != 0:
        return SimulationError(
            simulation_type="EQE",
            return_value=return_value,
            message=message,
        )


def run_simulation(session_path: Path, do_EQE: bool):
    simss_device_parameters = session_path / "simulation_setup_custom.txt"

    simulation_result = run_jV_simulation(
        session_path,
        simss_device_parameters=simss_device_parameters,
    )
    if isinstance(simulation_result, SimulationError):
        shutil.rmtree(session_path)
        return f"ERROR: {simulation_result}"

    device_characteristics_file, jv_file = process_jV_simulation_output(session_path)

    plot_jV_simulation_output(
        session_path=session_path,
        device_characteristics_file=device_characteristics_file,
        jv_file=jv_file,
    )

    if not do_EQE:
        os.system(f"rm {session_path}/simss")
        return "DONE"

    simulation_result = run_EQE_simulation(
        session_path,
        simss_device_parameters=simss_device_parameters,
    )

    if isinstance(simulation_result, SimulationError):
        shutil.rmtree(session_path)
        return f"ERROR: {simulation_result}"

    plot_EQE(
        **asdict(EQEData.from_file(f"{session_path}/EQE.dat")),
        dpi=DPI,
        save_path=session_path,
    )

    sim.delete_folders("tmp", session_path)
    os.remove(session_path / "simss")
    return "DONE"


def run(
    simulation_dir: Path,
    process_id: int,
    manual_randomization: bool,
    do_EQE: bool,
):
    np.random.seed(process_id)
    session_path = simulation_dir / f"run_{process_id}"
    set_up_session(PATH_TO_SIMSS, session_path)
    prepare_simulation(session_path, manual_randomization)
    return run_simulation(session_path, do_EQE)


def run_parallel(
    simulation_dir: Path,
    num_workers: int = 4,
    num_todo: int = 4,
    manual_randomization: bool = False,
    do_EQE: bool = True,
):
    parallel = joblib.Parallel(num_workers, return_as="generator")

    jobs = [
        joblib.delayed(run)(simulation_dir, process_id, manual_randomization, do_EQE)
        for process_id in range(num_todo)
    ]

    results = list(parallel(jobs))

    with open(simulation_dir / "results.txt", "w") as outfile:
        for idx, x in enumerate(results):
            outfile.write(f"{idx} {x}\n")


def main():
    simulation_dir = Path.cwd() / "simulation_runs"
    run_parallel(simulation_dir, num_workers=1, num_todo=1)


if __name__ == "__main__":
    main()
