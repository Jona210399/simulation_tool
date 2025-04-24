from simulation_tool.typing_ import PathLike


def get_layer3_config(path_to_simss: PathLike) -> dict[str, str | int | float]:
    path_to_simss = str(path_to_simss)
    layer3_config = {
        # General
        "l3.L": 3e-8,
        "l3.eps_r": 3.5,
        "l3.E_c": 2.5,
        "l3.E_v": 5,
        "l3.N_c": 1e26,
        "l3.N_D": 0,
        "l3.N_A": 0,
        # Mobilities
        "l3.mu_n": 1e-8,
        "l3.mu_p": 1e-8,
        "l3.mobnDep": 0,
        "l3.mobpDep": 0,
        "l3.gamma_n": 0,
        "l3.gamma_p": 0,
        # Interface-layer-to-right
        "l3.nu_int_n": 1e3,
        "l3.nu_int_p": 1e3,
        "l3.N_t_int": 0,
        "l3.E_t_int": 4,
        "l3.intTrapFile": "none",
        "l3.intTrapType": 1,
        "l3.C_n_int": 1e-13,
        "l3.C_p_int": 1e-13,
        # Ions
        "l3.N_anion": 0,
        "l3.N_cation": 0,
        "l3.mu_anion": 1e-14,
        "l3.mu_cation": 1e-14,
        "l3.ionsMayEnter": 0,
        # Generation and recombination
        "l3.G_ehp": 0,
        "l3.layerGen": 0,
        "l3.nkLayer": f"{path_to_simss}/../Data/nk_Ca.txt",
        "l3.fieldDepG": 0,
        "l3.P0": 0,
        "l3.a": 1e-9,
        "l3.thermLengDist": 2,
        "l3.k_f": 1e6,
        "l3.k_direct": 1e-18,
        "l3.preLangevin": 1,
        "l3.useLangevin": 1,
        # Bulk trapping
        "l3.N_t_bulk": 0,
        "l3.C_n_bulk": 1e-13,
        "l3.C_p_bulk": 1e-13,
        "l3.E_t_bulk": 4,
        "l3.bulkTrapFile": "none",
        "l3.bulkTrapType": -1,
    }
    return layer3_config
