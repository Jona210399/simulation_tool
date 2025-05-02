from dataclasses import asdict
from pathlib import Path

import joblib
import numpy as np

import simulation_tool.templates.to_text as to_text
from simulation_tool import SimulationError
from simulation_tool.constants import NANO, RUN_DIR_PREFIX
from simulation_tool.data import BandDiagramData, OpticalData
from simulation_tool.EQE import (
    EQEParameters,
    create_EQE_simulation_plots,
    run_EQE_simulation,
)
from simulation_tool.jV import (
    create_jV_simulation_plots,
    preserve_jV_simulation_output,
    run_jV_simulation,
)
from simulation_tool.session import clean_session, delete_session, set_up_session
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.randomization import (
    randomize_layer1,
    randomize_layer2,
    randomize_layer3,
    randomize_simss_config,
)
from simulation_tool.templates.simss import SimssConfig
from simulation_tool.utils import save_figure

PATH_TO_SIMSS = Path(
    "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/SimSS"
)
DPI = 60
ADJUST_BAND_DIAGRAM_Y_AXIS = True
WAVE_LENGTH_STEP = NANO  # Used for optical data generation
LAMBDA_STEP = 10 * NANO  # Used for EQE simulation


def create_layer_stack(
    session_path: Path,
    path_to_simss: Path,
    bandgap: float,
) -> list[Layer]:
    layer1 = Layer.get_default_layer1(path_to_simss=path_to_simss)
    layer2 = Layer.get_default_layer2(
        session_path=session_path,
        bandgap=bandgap,
    )
    layer3 = Layer.get_default_layer3(path_to_simss=path_to_simss)

    return [layer1, layer2, layer3]


def write_simulation_input_files(
    simss_config: SimssConfig,
    layers: list[Layer],
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
    path_to_simss: Path,
    randomized: bool,
) -> SimssConfig:
    simss_config = SimssConfig.from_session(
        session_path=session_path,
        path_to_simss=path_to_simss,
    )

    optical_data = OpticalData.new(
        min_wavelenght=simss_config.optics.lambda_min,
        max_wavelenght=simss_config.optics.lambda_max,
        wavelenght_step=WAVE_LENGTH_STEP,
    )

    save_figure(optical_data, save_path=session_path / "nk.png", dpi=DPI)
    optical_data.save(session_path)

    layers = create_layer_stack(
        session_path=session_path,
        path_to_simss=path_to_simss,
        bandgap=optical_data.bandgap,
    )

    layer1, layer2, layer3 = layers

    if randomized:
        layer2 = randomize_layer2(layer=layer2, bandgap=optical_data.bandgap)
        layer1 = randomize_layer1(
            layer=layer1,
            layer2_E_c=layer2.E_c,
            layer2_E_v=layer2.E_v,
        )
        layer3 = randomize_layer3(
            layer=layer3,
            layer2_E_c=layer2.E_c,
            layer2_E_v=layer2.E_v,
        )
        simss_config = randomize_simss_config(
            simss_config=simss_config,
            layer1_E_c=layer1.E_c,
            layer3_E_v=layer3.E_v,
        )

    write_simulation_input_files(simss_config, layers)

    BandDiagramData.from_layers(
        W_L=simss_config.contacts.W_L,
        W_R=simss_config.contacts.W_R,
        L_TCO=simss_config.optics.L_TCO,
        L_BE=simss_config.optics.L_BE,
        layers=layers,
    ).plot(
        dpi=DPI,
        adjust_y_axis=ADJUST_BAND_DIAGRAM_Y_AXIS,
        save_path=session_path,
    )

    return simss_config


def run_simulations(
    session_path: Path,
    path_to_executable: Path,
    setup_file: Path,
    eqe_parameters: EQEParameters,
) -> str:
    simulation_result = run_jV_simulation(
        session_path=session_path,
        setup_file=setup_file,
        path_to_executable=path_to_executable,
    )
    if isinstance(simulation_result, SimulationError):
        delete_session(session_path)
        return f"ERROR {simulation_result}"

    create_jV_simulation_plots(
        session_path=session_path,
        dpi=DPI,
    )
    preserve_jV_simulation_output(
        session_path=session_path,
    )

    simulation_result = run_EQE_simulation(
        session_path=session_path,
        setup_file=setup_file,
        path_to_executable=path_to_executable,
        **asdict(eqe_parameters),
    )

    if isinstance(simulation_result, SimulationError):
        delete_session(session_path)
        return f"ERROR {simulation_result}"

    create_EQE_simulation_plots(session_path=session_path, dpi=DPI)

    return "DONE"


def run(
    simulation_dir: Path,
    process_id: int,
    randomized: bool,
):
    np.random.seed(process_id)
    session_path = simulation_dir / f"{RUN_DIR_PREFIX}_{process_id}"
    set_up_session(
        path_to_executable=PATH_TO_SIMSS / "simss",
        session_path=session_path,
    )

    smiss_config = prepare_simulation(
        session_path=session_path,
        path_to_simss=PATH_TO_SIMSS,
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
        path_to_executable=PATH_TO_SIMSS / "simss",
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
    parallel = joblib.Parallel(num_workers, return_as="generator")

    jobs = [
        joblib.delayed(run)(simulation_dir, process_id, randomized)
        for process_id in range(num_todo)
    ]

    results = list(parallel(jobs))

    with open(simulation_dir / "results.txt", "w") as outfile:
        for idx, x in enumerate(results):
            outfile.write(f"{idx} {x}\n")


def main():
    simulation_dir = Path.cwd() / "simulation_runs"
    run_parallel(simulation_dir, num_workers=1, num_todo=1, randomized=False)


if __name__ == "__main__":
    main()
