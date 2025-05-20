import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
from tqdm import tqdm

import simulation_tool.templates.to_text as to_text
from simulation_tool import DeviceParametersIncompleteError, SimulationError
from simulation_tool.constants import (
    DPI,
    LAMBDA_STEP,
    ORIGINAL_RANDOMIZATION,
    PATH_TO_SIMSS_EXECUTABLE,
    PLOTTING_ENABLED,
    RUN_DIR_PREFIX,
    WAVE_LENGTH_STEP,
)
from simulation_tool.data import BandDiagramData, OpticalData
from simulation_tool.EQE import (
    EQEParameters,
    run_EQE_simulation,
)
from simulation_tool.jV import (
    run_jV_simulation,
)
from simulation_tool.randomization.kf_adapted import (
    randomize_device as randomize_device_kf_adapted,  # noqa: F401
)
from simulation_tool.randomization.original import (
    randomize_device as randomize_device_original,
)
from simulation_tool.session import clean_session, delete_session, set_up_session
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig, UserInterface
from simulation_tool.utils import save_figure

randomize_device = (
    randomize_device_original if ORIGINAL_RANDOMIZATION else randomize_device_kf_adapted
)


def create_layer_stack(
    session_path: Path,
    bandgap: float,
) -> tuple[Layer, Layer, Layer]:
    layer1 = Layer.get_default_layer1()
    layer2 = Layer.get_default_layer2(
        session_path=session_path,
        bandgap=bandgap,
    )
    layer3 = Layer.get_default_layer3()

    return layer1, layer2, layer3


def write_simulation_input_files(
    simss_config: SimssConfig,
    layers: tuple[Layer, Layer, Layer],
):
    with open(simss_config.setup_file, "w") as outfile:
        outfile.write(to_text.simms_config_to_text(simss_config=simss_config))

    simss_config.to_json_file(simss_config.setup_file.with_suffix(".json"))

    layer_files = [
        simss_config.layers.l1,
        simss_config.layers.l2,
        simss_config.layers.l3,
    ]

    for layer_file, layer in zip(layer_files, layers):
        with open(layer_file, "w") as outfile:
            outfile.write(to_text.layer_to_text(layer=layer))

        layer.to_json_file(layer_file.with_suffix(".json"))


def prepare_simulation(
    session_path: Path,
    randomized: bool,
) -> SimssConfig:
    simss_config = SimssConfig.from_session(session_path=session_path)

    optical_data = OpticalData.new(
        min_wavelenght=simss_config.optics.lambda_min,
        max_wavelenght=simss_config.optics.lambda_max,
        wavelenght_step=WAVE_LENGTH_STEP,
    )

    optical_data.save(session_path)
    optical_data.to_parquet(session_path)

    layers = create_layer_stack(
        session_path=session_path,
        bandgap=optical_data.bandgap,
    )

    if randomized:
        layers, simss_config = randomize_device(
            *layers, simss_config=simss_config, bandgap=optical_data.bandgap
        )

    write_simulation_input_files(simss_config, layers)

    if PLOTTING_ENABLED:
        save_figure(optical_data, save_path=session_path / "nk.png", dpi=DPI)

        save_figure(
            BandDiagramData.from_layers(
                W_L=simss_config.contacts.W_L,
                W_R=simss_config.contacts.W_R,
                L_TCO=simss_config.optics.L_TCO,
                L_BE=simss_config.optics.L_BE,
                layers=layers,
            ),
            save_path=session_path / "band_diagram.png",
            dpi=DPI,
        )

    return simss_config


def run_simulations(
    session_path: Path,
    path_to_executable: Path,
    simss_ui: UserInterface,
    setup_file: Path,
    eqe_parameters: EQEParameters,
) -> str:
    simulation_result = run_jV_simulation(
        session_path=session_path,
        path_to_executable=path_to_executable,
        simss_ui=simss_ui,
        setup_file=setup_file,
        plot_output=PLOTTING_ENABLED,
        plot_dpi=DPI,
    )
    if isinstance(
        simulation_result, (SimulationError, DeviceParametersIncompleteError)
    ):
        delete_session(session_path)
        return f"ERROR {simulation_result}"

    simulation_result = run_EQE_simulation(
        session_path=session_path,
        setup_file=setup_file,
        path_to_executable=path_to_executable,
        simss_ui=simss_ui,
        eqe_parameters=eqe_parameters,
        plot_output=PLOTTING_ENABLED,
        plot_dpi=DPI,
    )

    if isinstance(
        simulation_result, (SimulationError, DeviceParametersIncompleteError)
    ):
        delete_session(session_path)
        return f"ERROR {simulation_result}"

    return "DONE"


def run(
    simulation_dir: Path,
    process_id: int,
    randomized: bool,
    random_seed: int,
) -> str:
    np.random.seed(random_seed)
    session_path = simulation_dir / f"{RUN_DIR_PREFIX}_{process_id}"
    set_up_session(
        path_to_executable=PATH_TO_SIMSS_EXECUTABLE,
        session_path=session_path,
    )

    smiss_config = prepare_simulation(
        session_path=session_path,
        randomized=randomized,
    )

    eqe_parameters = EQEParameters(
        spectrum=smiss_config.optics.spectrum,
        lambda_min=smiss_config.optics.lambda_min,
        lambda_max=smiss_config.optics.lambda_max,
        lambda_step=LAMBDA_STEP,
    )

    result = run_simulations(
        session_path=session_path,
        path_to_executable=PATH_TO_SIMSS_EXECUTABLE,
        simss_ui=smiss_config.ui,
        setup_file=smiss_config.setup_file,
        eqe_parameters=eqe_parameters,
    )

    clean_session(
        session_path=session_path,
        simss_config=smiss_config,
    )

    return result


def run_parallel(
    simulation_dir: Path,
    num_workers: int = 4,
    num_todo: int = 4,
    randomized: bool = False,
):
    print(f"Number of simulations to run: {num_todo}")
    print(f"Number of workers: {num_workers}")

    random_seed = int(os.getenv("SLURM_JOB_ID", 0))

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(
                run,
                simulation_dir,
                process_id,
                randomized,
                random_seed + process_id,
            )
            for process_id in range(num_todo)
        ]
        error_count = 0
        with tqdm(total=num_todo) as pbar:
            for future in as_completed(futures):
                pbar.update(1)
                result = future.result()
                if result.startswith("ERROR"):
                    error_count += 1

    print(f"Number of errors: {error_count}")
    print(f"Number of successful runs: {num_todo - error_count}")
    print(
        f"Percentage of successful runs: {(num_todo - error_count) / num_todo * 100:.2f}%"
    )


def main():
    identifier = os.getenv("SLURM_JOB_ID", "")
    simulation_dir = Path.cwd() / "simulation_runs" / identifier
    num_workers = int(os.getenv("SLURM_CPUS_PER_TASK", 1))
    num_todo = 10_000

    start = datetime.now()
    print(f"Simulation started at: {start}")
    run_parallel(
        simulation_dir,
        num_workers=num_workers,
        num_todo=num_todo,
        randomized=True,
    )
    end = datetime.now()
    print(f"Simulation ended at: {end}")
    print(f"Total time: {end - start}")


if __name__ == "__main__":
    main()
