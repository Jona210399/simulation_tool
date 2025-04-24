from simulation_tool.typing_ import PathLike


def get_layer1_config(path_to_simss: PathLike) -> dict[str, str | int | float]:
    layer1_config = {
        # General
        "l1.L": 3e-8,
        "l1.eps_r": 5,
        "l1.E_c": 3.0,
        "l1.E_v": 5.5,
        "l1.N_c": 1e26,
        "l1.N_D": 0,
        "l1.N_A": 0,
        # Mobilities
        "l1.mu_n": 1e-6,
        "l1.mu_p": 1e-6,
        "l1.mobnDep": 0,
        "l1.mobpDep": 0,
        "l1.gamma_n": 0,
        "l1.gamma_p": 0,
        # Interface-layer-to-right
        "l1.nu_int_n": 1e3,
        "l1.nu_int_p": 1e3,
        "l1.N_t_int": 0,
        "l1.E_t_int": 4,
        "l1.intTrapFile": "none",
        "l1.intTrapType": -1,
        "l1.C_n_int": 1e-13,
        "l1.C_p_int": 1e-13,
        # Ions
        "l1.N_anion": 0,
        "l1.N_cation": 0,
        "l1.mu_anion": 1e-14,
        "l1.mu_cation": 1e-14,
        "l1.ionsMayEnter": 0,
        # Generation and recombination
        "l1.G_ehp": 0,
        "l1.layerGen": 0,
        "l1.nkLayer": f"{path_to_simss}/../Data/nk_PEDOT.txt",
        "l1.fieldDepG": 0,
        "l1.P0": 0,
        "l1.a": 1e-9,
        "l1.thermLengDist": 2,
        "l1.k_f": 1e6,
        "l1.k_direct": 1e-18,
        "l1.preLangevin": 1,
        "l1.useLangevin": 1,
        # Bulk trapping
        "l1.N_t_bulk": 0,
        "l1.C_n_bulk": 1e-13,
        "l1.C_p_bulk": 1e-13,
        "l1.E_t_bulk": 4,
        "l1.bulkTrapFile": "none",
        "l1.bulkTrapType": -1,
    }
    return layer1_config
