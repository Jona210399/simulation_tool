from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import constants

from simulation_tool.constants import (
    EQE_SIM_OUTPUT_FILE_NAME,
    M_TO_NM,
    NANO,
    NUM_PHOTONS,
)
from simulation_tool.EQE.data import EQEData
from simulation_tool.exceptions import DeviceParametersIncompleteError, SimulationError
from simulation_tool.simsalabim.run import clean_up_simulation_output, run_simulation
from simulation_tool.templates.simss import UserInterface
from simulation_tool.utils import save_figure

h = constants.h
c = constants.c
q = constants.e

FatalError = SimulationError | DeviceParametersIncompleteError


@dataclass
class EQEParameters:
    spectrum: Path
    lambda_min: float
    lambda_max: float
    lambda_step: float


def run_EQE_simulation(
    session_path: Path,
    setup_file: Path,
    path_to_executable: Path,
    simss_ui: UserInterface,
    eqe_parameters: EQEParameters,
    plot_output: bool,
    plot_dpi: int,
) -> EQEData | SimulationError | DeviceParametersIncompleteError:
    scPars_file = session_path / simss_ui.scParsFile
    output_file = session_path / EQE_SIM_OUTPUT_FILE_NAME
    lambda_min = eqe_parameters.lambda_min
    lambda_max = eqe_parameters.lambda_max
    lambda_step = eqe_parameters.lambda_step
    spectrum = eqe_parameters.spectrum

    if not scPars_file.with_suffix(".txt").exists():
        # If we already have a scPars file, from a previous jV Simulation, we dont have to run it again to obtain base values for jsc and jsc_err
        result = run_simulation(
            session_path=session_path,
            cmd_params={"dev_par_file", str(setup_file)},
            path_to_executable=path_to_executable,
            simss_ui=simss_ui,
        )
        if isinstance(result, FatalError):
            return result

    jsc_base, jsc_err_base = get_jsc_and_jsc_err(scPars_file.with_suffix(".txt"))

    jscs, jsc_errs = [], []

    wavelengths = np.arange(lambda_min, lambda_max + lambda_step, lambda_step)

    spectrum_data = pl.read_csv(
        spectrum,
        separator=" ",
        truncate_ragged_lines=True,
    )

    for wavelenght in wavelengths:
        spectrum_tmp = (
            session_path
            / f"{spectrum.stem}_{int(wavelenght * M_TO_NM)}nm{spectrum.suffix}"
        )

        spectrum_data_tmp = modify_spectrum(spectrum_data, wavelenght, NUM_PHOTONS)
        spectrum_data_tmp.write_csv(spectrum_tmp, separator=" ")

        result = run_simulation(
            session_path=session_path,
            cmd_params={
                "dev_par_file": str(setup_file),
                "spectrum": str(spectrum_tmp),
            },
            path_to_executable=path_to_executable,
            simss_ui=simss_ui,
        )

        spectrum_tmp.unlink()

        if isinstance(result, FatalError):
            return result

        jsc, jsc_err = get_jsc_and_jsc_err(scPars_file)

        jscs.append(jsc)
        jsc_errs.append(jsc_err)

    eqe_result = calculate_EQE(
        np.array(jscs),
        np.array(jsc_errs),
        jsc_base,
        jsc_err_base,
        wavelengths,
        NUM_PHOTONS,
    )

    clean_up_simulation_output(
        session_path=session_path,
        simss_ui=simss_ui,
    )

    eqe_result.write_parquet(output_file.with_suffix(".parquet"))

    eqe = EQEData.from_dataframe(eqe_result)

    if plot_output:
        save_figure(
            eqe,
            save_path=output_file.with_suffix(".png"),
            dpi=plot_dpi,
        )

    return eqe


def modify_spectrum(
    spectrum_data: pl.DataFrame, wavelenght: float, num_photons: float
) -> pl.DataFrame:
    """Modify the spectrum data by adding a peak at the specified wavelength with a height set by the number of added photons.

    Parameters
    ----------
    spectrum_data : pl.DataFrame
        DataFrame containing the original spectrum data.
    wavelenght : float
        Wavelength at which the peak will be added.
    num_photons : float
        Number of photons that are added to the irradiance.

    Returns
    -------
    pl.DataFrame
        Modified spectrum data with the added peak.
    """

    lambda_diff = (spectrum_data["lambda"] - wavelenght).abs()
    min_idx = lambda_diff.arg_min()

    delta_I = num_photons * ((h * c) / (wavelenght)) / NANO

    spectrum_data = spectrum_data.with_columns(
        pl.when(pl.arange(0, spectrum_data.height) == min_idx)
        .then(spectrum_data["I"] + delta_I)
        .otherwise(spectrum_data["I"])
        .alias("I")
    )

    return spectrum_data


def get_jsc_and_jsc_err(path_to_scPars: Path) -> tuple[float, float]:
    data = pl.read_csv(
        path_to_scPars,
        separator=" ",
    )

    return data["Jsc"].item(), data["ErrJsc"].item()


def calculate_EQE(
    jscs: NDArray,
    jsc_errs: NDArray,
    jsc_base: float,
    jsc_err_base: float,
    wavelengths: NDArray,
    num_photons: int,
):
    """Calculate the EQE values for the monochromatic peaks at each wavelength. Based on the change in the short-circuit current and the number of added photons.

    Parameters
    ----------
    jscs : array
        Short-circuit current density for each monochromatic peak.
    jsc_errs : array
        Error in the short-circuit current density for each monochromatic peak.
    jsc_base : float
        Short-circuit current density for the normal spectrum.
    jsc_err_base : float
        Error in the short-circuit current density for the normal spectrum.
    wavelengths : array
        Array with the wavelengths at which the monochromatic peaks were created.
    num_photons : float
        Number of photons that are added to the irradiance.

    Returns
    -------
    pl.DataFrame
    """

    deltaJ = np.abs(jscs - jsc_base)
    deltaJerr = jsc_errs + jsc_err_base
    E_photon = h * c / wavelengths

    I_diff = num_photons * E_photon
    EQE_val = deltaJ / ((q * I_diff) / E_photon)
    EQE_err = deltaJerr * ((q * I_diff) / E_photon)

    return pl.DataFrame(
        {
            "lambda": wavelengths,
            "Jext": jscs,
            "Jerr": jsc_errs,
            "deltaJ": deltaJ,
            "deltaJerr": deltaJerr,
            "Imonopeak": I_diff,
            "EQE": EQE_val,
            "EQEerr": EQE_err,
        }
    )
