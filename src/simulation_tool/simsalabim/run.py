from pathlib import Path
from subprocess import run

from simulation_tool.exceptions import SimulationError


def construct_command(
    path_to_executable: Path,
    cmd_pars: dict[str, dict[str, str]],
):
    """Construct a single string to use as command to run a SIMsalabim executable."""

    cmd_line = str(path_to_executable)

    # Check whether a device parameters file has been defined
    for i in cmd_pars:
        # When specified, the device parameter file must be placed first, as is required by SIMsalabim
        if i["par"] == "dev_par_file":
            args_single = " " + i["val"] + " "
            cmd_line = cmd_line + args_single
            # After the dev_par_file key had been found once, stop the loop. If more than one dev_par_file is specified, the rest are ignored.
            break

    # Add the parameters
    for i in cmd_pars:
        if i["par"] != "dev_par_file":
            # Add each parameter as " -par_name par_value"
            args_single = " -" + i["par"] + " " + i["val"]
            cmd_line = cmd_line + args_single

    return cmd_line


def run_simulation(
    session_path: Path,
    path_to_executable: Path,
    cmd_pars: dict[str, dict[str, str]],
) -> SimulationError | None:
    cmd_line = construct_command(
        path_to_executable=path_to_executable,
        cmd_pars=cmd_pars,
    )

    stdout = session_path / "stdout.txt"
    stderr = session_path / "stderr.txt"
    scPars_file = session_path / "scPars.dat"

    with open(stdout, "w") as stdout_file:
        with open(stderr, "w") as stderr_file:
            result = run(
                cmd_line,
                cwd=session_path,
                stdout=stdout_file,
                stderr=stderr_file,
                check=False,
                shell=True,
            )

    if result.returncode != 0:
        return SimulationError(
            message=f"Simulation failed with return code {result.returncode}. Check the output files for more details. Error",
        )

    if not scPars_file.exists():
        return SimulationError(
            message="Simulation did not produce the expected output files. Missing file: scPars.dat",
        )

    return None
