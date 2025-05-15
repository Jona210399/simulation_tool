import shutil
from pathlib import Path

from simulation_tool.constants import NK_FILE_NAME, PLOTTING_ENABLED
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
    if not PLOTTING_ENABLED:
        return
    plots_dir = session_path / "plots"
    plots_dir.mkdir(exist_ok=True)

    for png_file in session_path.glob("*.png"):
        destination = plots_dir / png_file.name
        shutil.move(png_file, destination)


def collect_session_data(session_path: Path):
    data_dir = session_path / "data"
    data_dir.mkdir(exist_ok=True)

    for data_file in session_path.glob("*.parquet"):
        destination = data_dir / data_file.name
        shutil.move(data_file, destination)

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

    layer_files = [
        simss_config.layers.l1,
        simss_config.layers.l2,
        simss_config.layers.l3,
    ]

    for file in layer_files + [session_path / NK_FILE_NAME]:
        file.unlink(missing_ok=True)

    collect_session_plots(session_path)
    collect_session_data(session_path)

    simss_config.setup_file.unlink(missing_ok=True)
