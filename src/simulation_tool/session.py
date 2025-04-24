import shutil
from pathlib import Path


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
