from pathlib import Path

from scipy.constants import nano as NANO

M_TO_NM: float = 1 / NANO
NM_TO_M: float = NANO

ORIGINAL_RANDOMIZATION: bool = False
"""Use the original randomization function. This is used to switch between the original and the kf-adapted randomization function."""

NUM_TODO: int = 5
"""Number of simulations to be run. Note that this is not the number of data points, but the number of simulations. The number of data points is determined by the number of successful simulations. Not all simulations will be successful, so the number of data points will be less than or equal to the number of simulations."""


SET_UP_FILE: str = "simulation_setup"
RUN_DIR_PREFIX: str = "run"

NUM_PHOTONS = 0.03 * 1e21
"""Number of photons that are added to the irradiance during EQE simulation. In this example it is 3% of the number of photons absorbed by a Silicon solar cell in m-2"""

PATH_TO_SIMSS_EXECUTABLE = Path(
    "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/SimSS/simss"
)
"""Path to the simss executable. This is used to run the simulation. The path should be absolute and point to the simss executable."""

PLOTTING_ENABLED = True
"""Enable plotting of the results. This is used to enable or disable plotting of the results."""
DPI = 60
"""DPI for the plots."""
ADJUST_BAND_DIAGRAM_Y_AXIS = True

USE_KK_CONSISTENT_GENERATION = True
"""Use Kramers-Kronig consistent generation of n and k."""
N_OFFSET_RANGE = (1.1, 1.8)
"""Range for the n offset. This is used to generate the refractive index n."""
NK_FILE_NAME = "nk.txt"
"""Name of the generated nk file. This is used to store the refractive index n and extinction coefficient k."""

WAVE_LENGTH_STEP = NANO
"""Step size for the wavelength in nm. This is used for the optical data generation."""
LAMBDA_STEP = 10 * NANO
"""Step size for the wavelength in nm. This is used for the EQE simulation."""

UVVIS_FILE_NAME = "uvvis"
"""Name of the generated UV-Vis file. This is used to store the UV-Vis data. A file ending will be added automatically."""

SIMULATION_TIMEOUT = 30
"""Timeout for the simulation in seconds. This is used to timeout the subprocess call to the simulation executable."""

EQE_SIM_OUTPUT_FILE_NAME = "EQE"
"""Name of the generated EQE file. This is used to store the EQE data. A file ending will be added automatically."""

JV_SIM_OUTPUT_FILE_NAME = "jV"
"""Name of the generated jV file. This is used to store the jV data. A file ending will be added automatically."""

ADDITIONAL_FILES_DIR = Path(__file__).parent / "additional_files"
"""Directory for the additional files. This is used to store the additional files that are needed for the simulation. For example nk.txt and spectrum.txt files"""
