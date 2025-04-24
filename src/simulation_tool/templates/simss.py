from simulation_tool.typing_ import PathLike


def get_simss_config(
    session_path: PathLike,
    path_to_simss: PathLike,
) -> dict[str, str | int | float]:
    simss_config = {
        # General
        "T": 295,
        # Layers
        "l1": f"{session_path}/L1_parameters.txt",
        "l2": f"{session_path}/L2_parameters.txt",
        "l3": f"{session_path}/L3_parameters.txt",
        # Contacts
        "W_L": 3.0,
        "W_R": 5.0,
        "S_n_L": -1e-7,
        "S_p_L": -1e-7,
        "S_n_R": -1e-7,
        "S_p_R": -1e-7,
        "R_shunt": -1,
        "R_series": 0,  # np.random.random()*0.1,
        # Optics
        "G_frac": 1.0,
        "genProfile": "calc",
        "L_TCO": 1.1e-7,
        "L_BE": 2e-7,
        "nkSubstrate": f"{path_to_simss}/../Data/nk_SiO2.txt",
        "nkTCO": f"{path_to_simss}/../Data/nk_ITO.txt",
        "nkBE": f"{path_to_simss}/../Data/nk_Al.txt",
        "spectrum": f"{path_to_simss}/../Data/AM15G.txt",
        "lambda_min": 3.5e-7,
        "lambda_max": 8e-7,
        # Numerical Parameters
        "NP": 400,
        "tolPois": 0.00001,
        "maxDelV": 10,
        "maxItPois": 1000,
        "maxItSS": 1000,
        "currDiffInt": 2,
        "tolDens": 1e-8,
        "couplePC": 4,
        "minAcc": 0.5,
        "maxAcc": 1.0,
        "ignoreNegDens": 1,
        "failureMode": 0,
        "grad": 1,
        # Voltage range
        "Vdist": 1,
        "preCond": 0,
        "Vpre": 0,
        "fixIons": 0,
        "Vscan": 1,
        "Vmin": -0.5,
        "Vmax": 1.7,
        "Vstep": 0.025,  # default 0.025
        "Vacc": 0,
        "NJV": 100,
        "untilVoc": 0,
        # User interface
        "timeout": -1,
        "pauseAtEnd": 0,
        "autoTidy": 1,
        "useExpData": 0,
        "expJV": "expJV.csv",
        "fitMode": "lin",
        "fitThreshold": 0.8,
        "JVFile": "JV.dat",
        "varFile": "Var.dat",
        "limitDigits": 1,
        "outputRatio": 1,
        "scParsFile": "scPars.dat",
        "logFile": "log.txt",
    }
    return simss_config
