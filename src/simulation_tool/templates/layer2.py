from simulation_tool.typing_ import PathLike


def get_layer2_config(
    session_path: PathLike,
    bandgap: float,
) -> dict[str, str | int | float]:
    """Bandgap is given in eV."""
    layer2_config = {
        # General
        "l2.L": 5.4e-8,
        "l2.eps_r": 4,
        "l2.E_c": 3.0,
        "l2.E_v": 3 + bandgap,
        "l2.N_c": 1e25,
        "l2.N_D": 0,
        "l2.N_A": 0,
        # Mobilities
        "l2.mu_n": 0.01,
        "l2.mu_p": 0.01,
        "l2.mobnDep": 0,
        "l2.mobpDep": 0,
        "l2.gamma_n": 0.0001,
        "l2.gamma_p": 0,
        # Interface-layer-to-right
        "l2.nu_int_n": 1e3,
        "l2.nu_int_p": 1e3,
        "l2.N_t_int": 0,
        "l2.E_t_int": 4,
        "l2.intTrapFile": "none",
        "l2.intTrapType": 1,
        "l2.C_n_int": 1e-13,
        "l2.C_p_int": 1e-13,
        # Ions
        "l2.N_anion": 0,
        "l2.N_cation": 0,
        "l2.mu_anion": 1e-14,
        "l2.mu_cation": 1e-14,
        "l2.ionsMayEnter": 1,
        # Generation and recombination
        "l2.G_ehp": 7.25e27,
        "l2.layerGen": 1,
        "l2.nkLayer": f"{session_path}/nk.txt",
        "l2.fieldDepG": 0,
        "l2.P0": 0,
        "l2.a": 1e-9,
        "l2.thermLengDist": 2,
        "l2.k_f": 1e6,
        "l2.k_direct": 1e-18,
        "l2.preLangevin": 1,
        "l2.useLangevin": 1,
        # Bulk trapping
        "l2.N_t_bulk": 4.42e21,
        "l2.C_n_bulk": 1e-13,
        "l2.C_p_bulk": 1e-13,
        "l2.E_t_bulk": 4,
        "l2.bulkTrapFile": "none",
        "l2.bulkTrapType": -1,
    }
    return layer2_config
