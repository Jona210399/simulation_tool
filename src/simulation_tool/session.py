import os
import shutil
from pathlib import Path

import pySIMsalabim as sim


def set_up_session(
    path_to_simss: Path,
    session_path: Path,
):
    if not path_to_simss.exists():
        raise FileNotFoundError(
            f"Path to simss executable: {path_to_simss} does not exist."
        )

    if not session_path.exists():
        session_path.mkdir(parents=True, exist_ok=True)

    shutil.copy(path_to_simss / "simss", session_path)


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


def clean_session(session_path: Path):
    if not session_path.exists():
        return

    collect_session_plots(session_path)
    collect_session_data(session_path)
    os.remove(session_path / "Var.dat")
    sim.delete_folders("tmp", session_path)
    os.remove(session_path / "simss")
