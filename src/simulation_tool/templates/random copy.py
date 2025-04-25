import numpy as np

from simulation_tool.templates.simss import SimssConfig
from simulation_tool.typing_ import PathLike
from simulation_tool.utils import loguniform


def sample_random_parameters(
    session_path: PathLike,
    path_to_simss: PathLike,
    bandgap: float,
):
    """Bandgap is given in eV"""
    settings = {}

    ### hierarchy: From active layer to l1/l3, and then to the electrodes
    # l2: active layer
    l2_E_c = np.random.randn() * 0.5 + 3.0
    l2_E_v = l2_E_c + bandgap
    # l1: l1_E_c above l2_E_c: don't block holes
    l1_E_c = np.random.randn() * 0.1 + l2_E_c + 0.1
    # l1: l1_E_v above l2_E_v: block electrons
    l1_E_v = np.random.randn() * 0.1 + l2_E_v + 0.1
    # l3: l3_E_c below l2_E_c: block holes
    l3_E_c = np.random.randn() * 0.1 + l2_E_c - 0.1
    # l3: l3_E_v below l2_E_v: don't block electrons
    l3_E_v = np.random.randn() * 0.1 + l2_E_v - 0.1
    # electrodes: W_L equal or above l1_E_c: don't block holes --> hardcoded in simss
    W_L = np.random.randn() * 0.5 + l1_E_c
    if W_L < l1_E_c:
        W_L = l1_E_c
    # electrodes: W_R equal or below l3_E_v: don't block electrons --> hardcoded in simss
    W_R = np.random.randn() * 0.5 + l3_E_v
    if W_R > l3_E_v:
        W_R = l3_E_v

    ### General settings
    # Contacts
    settings["W_L"] = W_L  # default 3
    settings["W_R"] = W_R  # default 5
    settings["S_n_L"] = -1e-7  # default -1e-7
    settings["S_p_L"] = -1e-7  # default -1e-7
    settings["S_n_R"] = -1e-7  # default -1e-7
    settings["S_p_R"] = -1e-7  # default -1e-7
    settings["R_shunt"] = -1  # default -1
    settings["R_series"] = 0  # default 0

    # Optics
    settings["G_frac"] = (
        1  # np.clip(np.random.randn()*0.3+1.0, a_min=0.0, a_max=1.0) # default 1
    )
    settings["genProfile"] = "calc"  #
    settings["L_TCO"] = loguniform(
        low=1e-8, high=3e-7
    )  # default: 1.1e-7, if 0 it is not used
    settings["L_BE"] = loguniform(low=1e-8, high=3e-7)  # default 2e-7, has to be >0

    # Optics
    settings["Vstep"] = 0.025  # default 0.025

    # get the missing default parameters and add them here
    simss_config = SimssConfig.from_session(session_path, path_to_simss).to_dict()
    for key in simss_config.keys():
        if key not in settings.keys():
            settings[key] = simss_config[key]

    ### Layer settings
    # General

    settings["l3.L"] = loguniform(low=3e-9, high=3e-7)  # default 3E-8
    layerfactor_l3 = (np.log(settings["l3.L"]) - np.log(3e-9)) / (
        np.log(3e-7) - np.log(3e-9)
    )
    settings["l3.eps_r"] = np.random.randn() * 1.0 + 4.0  # default 3.5
    settings["l3.E_c"] = l3_E_c  # default 2.5
    settings["l3.E_v"] = l3_E_v  # default 5
    settings["l3.N_c"] = loguniform(low=1e24, high=1e27)  # default 1E26
    settings["l3.N_D"] = loguniform(low=1e16, high=1e21)  # default 0
    settings["l3.N_A"] = loguniform(low=1e16, high=1e21)  # default 0

    # Mobilities

    settings["l3.mu_n"] = 10.0 ** (
        -np.random.randn() * 2.0 - 7.0 + layerfactor_l3 * 2.0
    )  # loguniform(low = 1e-9, high = 1e-2) # default 1e-8
    settings["l3.mu_p"] = 10.0 ** (
        -np.random.randn() * 2.0 - 7.0 + layerfactor_l3 * 2.0
    )  # loguniform(low = 1e-9, high = 1e-2) # default 1e-8
    settings["l3.mobnDep"] = 0  # not randomized
    settings["l3.mobpDep"] = 0  # not randomized
    settings["l3.gamma_n"] = 0  # not randomized
    settings["l3.gamma_p"] = 0  # not randomized

    settings["l2.nu_int_n"] = 1e3  # default 1e3
    settings["l2.nu_int_p"] = 1e3  # default 1e3
    # l2_int_trap_probabilitly = 0.1
    # if l2_int_trap_probabilitly > np.random.random():
    #    l2_N_t_int = loguniform(low = 1e-12, high = 1e-8)
    # else:
    #    l2_N_t_int = 0
    settings["l2.N_t_int"] = 0  # default 0
    # max_E_c_23 = max(randomized_settings["l2.E_c"], randomized_settings["l3.E_c"])
    # min_E_v_23 = min(randomized_settings["l2.E_v"], randomized_settings["l3.E_v"])
    # l2_E_t_int = np.random.uniform(low = max_E_c_23 + 0.01, high = min_E_v_23 - 0.01)
    settings["l2.E_t_int"] = 4  # default 4
    settings["l2.intTrapFile"] = "none"  # not randomized
    settings["l2.intTrapType"] = 1  # not randomized

    settings["l3.nu_int_n"] = 1e3  # default 1e3
    settings["l3.nu_int_p"] = 1e3  # default 1e3
    # l3_int_trap_probabilitly = 0.1
    # if l3_int_trap_probabilitly > np.random.random():
    #    l3_N_t_int = loguniform(low = 1e-12, high = 1e-8)
    # else:
    #    l3_N_t_int = 0
    settings["l3.N_t_int"] = 0  # default 0
    settings["l3.E_t_int"] = 4  # default 4
    settings["l3.intTrapFile"] = "none"  # not randomized
    settings["l3.intTrapType"] = 1  # not randomized
    settings["l3.C_n_int"] = loguniform(low=1e-15, high=1e-11)  # default 1e-13
    settings["l3.C_p_int"] = loguniform(low=1e-15, high=1e-11)  # default 1e-13

    settings["l3.G_ehp"] = 0  # not randomized
    settings["l3.layerGen"] = 0  # not randomized
    settings["l3.nkLayer"] = f"{path_to_simss}/../Data/nk_Ca.txt"  # not randomized here
    settings["l3.fieldDepG"] = 0  # not randomized
    settings["l3.P0"] = 0  # not randomized
    settings["l3.a"] = 1e-9  # not randomized
    settings["l3.thermLengDist"] = 2  # not randomized
    settings["l3.k_f"] = 1e6  # not randomized
    settings["l3.k_direct"] = 1e-18  # not randomized
    settings["l3.preLangevin"] = 1  # not randomized
    settings["l3.useLangevin"] = 1  # not randomized

    # Bulk trapping

    # default 1e-13

    settings["l3.N_t_bulk"] = loguniform(low=1e17, high=1e23)  # default 0
    settings["l3.C_n_bulk"] = loguniform(low=1e-15, high=1e-11)  # default 1e-13
    settings["l3.C_p_bulk"] = loguniform(low=1e-15, high=1e-11)  # default 1e-13
    settings["l3.E_t_bulk"] = np.random.uniform(
        low=settings["l3.E_c"] + 0.01,
        high=settings["l3.E_v"] - 0.01,
    )  # default 4

    return settings
