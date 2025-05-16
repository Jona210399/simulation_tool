from pathlib import Path
from subprocess import TimeoutExpired, run

from simulation_tool.constants import SIMULATION_TIMEOUT
from simulation_tool.exceptions import DeviceParametersIncompleteError, SimulationError
from simulation_tool.templates.simss import SimssOutputFiles, UserInterface


def _construct_command(path_to_executable: Path, cmd_params: dict[str, str]) -> str:
    cmd = str(path_to_executable)
    cmd_params = cmd_params.copy()

    dev_par_file = "dev_par_file"
    if dev_par_file in cmd_params:
        file = cmd_params.pop(dev_par_file)
        cmd += f" {file}"

    for key, value in cmd_params.items():
        cmd += f" -{key} {value}"

    return cmd


def run_simulation(
    session_path: Path,
    path_to_executable: Path,
    cmd_params: dict[str, str],
    simss_ui: UserInterface,
) -> SimulationError | DeviceParametersIncompleteError | None:
    cmd_line = _construct_command(
        path_to_executable=path_to_executable,
        cmd_params=cmd_params,
    )

    stdout = session_path / "sim.out"
    stderr = session_path / "sim.err"
    scPars_file = session_path / simss_ui.scParsFile

    with open(stdout, "w") as stdout_file:
        with open(stderr, "w") as stderr_file:
            try:
                result = run(
                    cmd_line,
                    cwd=session_path,
                    stdout=stdout_file,
                    stderr=stderr_file,
                    check=False,
                    shell=True,
                    timeout=SIMULATION_TIMEOUT,
                )
            except TimeoutExpired:
                return SimulationError(
                    message="Simulation timed out. Check the output files for more details.",
                )

    if result.returncode != 0:
        return SimulationError(
            message=f"Simulation failed with return code {result.returncode}. Check the output files {stdout.name} and {stderr.name} for more details.",
        )

    if not scPars_file.exists():
        return DeviceParametersIncompleteError(
            message=f"Simulation did not produce the expected output files. Missing file: {simss_ui.scParsFile}",
        )

    for file in SimssOutputFiles.get_all() + [simss_ui.JVFile]:
        if not (session_path / file).exists():
            return SimulationError(
                message=f"Simulation did not produce the expected output files. Missing file: {file}",
            )

    return None


def clean_up_simulation_output(
    session_path: Path,
    simss_ui: UserInterface,
) -> None:
    for file in SimssOutputFiles.get_all() + simss_ui.get_files():
        (session_path / file).unlink(missing_ok=True)
