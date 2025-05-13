import shutil
from pathlib import Path

import pySIMsalabim as sim

from simulation_tool.templates.simss import SimssConfig


def set_up_session(
    path_to_executable: Path,
    session_path: Path,
):
    if not path_to_executable.exists():
        raise FileNotFoundError(
            f"Path to  xecutable: {path_to_executable} does not exist."
        )

    session_path.mkdir(parents=True, exist_ok=True)


def collect_session_plots(session_path: Path):
    plots_dir = session_path / "plots"
    plots_dir.mkdir(exist_ok=True)

    for png_file in session_path.glob("*.png"):
        destination = plots_dir / png_file.name
        shutil.move(png_file, destination)


def collect_session_data(session_path: Path):
    data_dir = session_path / "data"
    data_dir.mkdir(exist_ok=True)

    for data_file in session_path.glob("*.txt"):
        destination = data_dir / data_file.name
        shutil.move(data_file, destination)


def delete_session(session_path: Path):
    if session_path.exists():
        shutil.rmtree(session_path)


def clean_session(
    session_path: Path,
    simss_config: SimssConfig,
):
    if not session_path.exists():
        return

    (session_path / simss_config.ui.logFile).unlink(missing_ok=True)

    collect_session_plots(session_path)
    collect_session_data(session_path)

    (session_path / simss_config.ui.varFile).unlink(missing_ok=True)
    simss_config.setup_file.unlink(missing_ok=True)
    (session_path / simss_config.ui.JVFile).unlink(missing_ok=True)
    (session_path / simss_config.ui.scParsFile).unlink(missing_ok=True)

    layer_files = [
        simss_config.layers.l1,
        simss_config.layers.l2,
        simss_config.layers.l3,
    ]

    for layer_file in layer_files:
        layer_file.unlink(missing_ok=True)

    sim.delete_folders("tmp", session_path)
