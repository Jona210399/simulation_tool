import contextlib
from pathlib import Path

import joblib


@contextlib.contextmanager
def joblib_progress_to_file(
    simulation_dir: Path,
    total: int,
    description: str | None = None,
    width: int = 40,
):
    """
    Context manager to track joblib.Parallel progress and write updates
    to a single-line progress.txt file in the given session_path.

    Parameters:
        simulation_dir (Path): Directory where progress.txt will be written.
        total (int): Total number of tasks to track.
        description (str): Optional label for the progress bar.
        width (int): Width of the progress bar.
    """
    simulation_dir.mkdir(parents=True, exist_ok=True)
    progress_file_path = simulation_dir / "progress.txt"

    if description is None:
        description = "Processing"

    original_print_progress = joblib.parallel.Parallel.print_progress

    def update_progress(self: joblib.parallel.Parallel):
        completed = self.n_completed_tasks
        fraction = completed / total if total else 0
        done = int(fraction * width)
        bar = "#" * done + "-" * (width - done)
        percent = int(fraction * 100)

        line = f"{description}: [{bar}] {percent}% ({completed}/{total})"

        with progress_file_path.open("r+") as f:
            f.seek(0)
            f.write(line)
            f.truncate()

        return original_print_progress(self)

    try:
        with progress_file_path.open("w") as f:
            f.write(f"{description}: [{' ' * width}] 0% (0/{total})")

        joblib.parallel.Parallel.print_progress = update_progress
        yield

        # Write final progress line
        final_line = f"{description}: [{'#' * width}] 100% ({total}/{total})"
        with progress_file_path.open("w") as f:
            f.write(final_line + "\n")

    finally:
        joblib.parallel.Parallel.print_progress = original_print_progress
