from scipy.constants import nano as NANO

M_TO_NM: float = 1 / NANO
NM_TO_M: float = NANO


SET_UP_FILE: str = "simulation_setup"
RUN_DIR_PREFIX: str = "run"

NUM_PHOTONS = 0.03 * 1e21
"""Number of photons that are added to the irradiance during EQE simulation. In this example it is 3% of the number of photons absorbed by a Silicon solar cell in m-2"""
