from pathlib import Path

import joblib
import numpy as np

import simulation_tool.templates.to_text as to_text
from simulation_tool import SimulationError
from simulation_tool.data import BandDiagramData, OpticalData
from simulation_tool.EQE import (
    EQEParameters,
    create_EQE_simulation_plots,
    preserve_EQE_simulation_output,
    run_EQE_simulation,
)
from simulation_tool.jV import (
    create_jV_simulation_plots,
    preserve_jV_simulation_output,
    run_jV_simulation,
)
from simulation_tool.session import clean_session, delete_session, set_up_session
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.random import sample_random_parameters
from simulation_tool.templates.simss import SimssConfig

PATH_TO_SIMSS = Path(
    "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/SimSS"
)
DPI = 60
SPECTRUM = "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/Data/AM15G.txt"

EQE_PARAMETERS = EQEParameters(
    spectrum=SPECTRUM,
    lambda_min=350,
    lambda_max=900,
    lambda_step=10,
    Vext=0,
)


def create_settings(
    session_path: Path,
    path_to_simss: Path,
    bandgap: float,
    manual_randomization: bool,
) -> dict[str, str | int | float]:
    if manual_randomization:
        simss_config = SimssConfig.from_session(
            session_path=session_path, path_to_simss=path_to_simss
        ).to_dict()
        layer1 = Layer.get_default_layer1(path_to_simss).to_dict(prefix="l1")
        layer2 = Layer.get_default_layer2(session_path, bandgap).to_dict(prefix="l2")
        layer3 = Layer.get_default_layer3(path_to_simss).to_dict(prefix="l3")

        settings = {}

        configs = [simss_config, layer1, layer2, layer3]
        for config in configs:
            settings.update(config)
    else:
        settings = sample_random_parameters(
            session_path,
            path_to_simss,
            bandgap,
        )
    simss_config = to_text.get_simss_config(settings)
    layer1 = to_text.get_layer_config(settings, "l1")
    layer2 = to_text.get_layer_config(settings, "l2")
    layer3 = to_text.get_layer_config(settings, "l3")

    with open(session_path / "simulation_setup.txt", "w") as outfile:
        outfile.write(simss_config)

    with open(session_path / "L1_parameters.txt", "w") as outfile:
        outfile.write(layer1)

    with open(session_path / "L2_parameters.txt", "w") as outfile:
        outfile.write(layer2)

    with open(session_path / "L3_parameters.txt", "w") as outfile:
        outfile.write(layer3)

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
    BandDiagramData.from_settings(settings).plot(
        dpi=DPI,
        save_path=session_path,
    )


def run_simulations(session_path: Path):
    simss_device_parameters = session_path / "simulation_setup.txt"

    simulation_result = run_jV_simulation(
        session_path,
        simss_device_parameters=simss_device_parameters,
    )
    if isinstance(simulation_result, SimulationError):
        delete_session(session_path)
        return f"ERROR: {simulation_result}"

    create_jV_simulation_plots(session_path=session_path, dpi=DPI)
    preserve_jV_simulation_output(session_path=session_path)

    simulation_result = run_EQE_simulation(
        session_path,
        simss_device_parameters=simss_device_parameters,
        eqe_parameters=EQE_PARAMETERS,
    )

    if isinstance(simulation_result, SimulationError):
        delete_session(session_path)
        return f"ERROR: {simulation_result}"

    create_EQE_simulation_plots(session_path=session_path, dpi=DPI)
    preserve_EQE_simulation_output(session_path=session_path)

    return "DONE"


def run(
    simulation_dir: Path,
    process_id: int,
    manual_randomization: bool,
):
    np.random.seed(process_id)
    session_path = simulation_dir / f"run_{process_id}"
    set_up_session(PATH_TO_SIMSS, session_path)
    prepare_simulation(session_path, manual_randomization)
    result = run_simulations(session_path)
    clean_session(session_path)
    return result


def run_parallel(
    simulation_dir: Path,
    num_workers: int = 4,
    num_todo: int = 4,
    manual_randomization: bool = False,
):
    parallel = joblib.Parallel(num_workers, return_as="generator")

    jobs = [
        joblib.delayed(run)(simulation_dir, process_id, manual_randomization)
        for process_id in range(num_todo)
    ]

    results = list(parallel(jobs))

    with open(simulation_dir / "results.txt", "w") as outfile:
        for idx, x in enumerate(results):
            outfile.write(f"{idx} {x}\n")


def main():
    simulation_dir = Path.cwd() / "simulation_runs"
    run_parallel(simulation_dir, num_workers=1, num_todo=1, manual_randomization=False)


if __name__ == "__main__":
    main()
