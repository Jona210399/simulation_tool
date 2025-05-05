from simulation_tool.constants import NANO
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig
from simulation_tool.utils import loguniform, randn, uniform

OPV = False
ENERGY_LEVEL_OFFSET = 0.5  # eV

N_C_RANGE = (1e25, 1e27)
N_D_A_RANGE = (1e16, 1e21)

L_RANGE = (3 * NANO, 300 * NANO)
L2_L_RANGE_OPV = (50 * NANO, 200 * NANO)
L2_L_RANGE_PEROVSKITE = (300 * NANO, 1000 * NANO)
L2_L_RANGE = L2_L_RANGE_OPV if OPV else L2_L_RANGE_PEROVSKITE

MOBILITY_RANGE = (1e-9, 1e-0)
L2_MOBILITY_RANGE_OPV = (5e-9, 1e-6)
L2_MOBILITY_RANGE_PEROVSKITE = (5e-9, 1e-0)
L2_MOBILITY_RANGE = L2_MOBILITY_RANGE_OPV if OPV else L2_MOBILITY_RANGE_PEROVSKITE

INTERFACE_C_N_P_RANGE = (1e-15, 1e-11)

L2_KF_RANGE = (1e4, 1e8)
L2_K_DIRECT_RANGE = (1e-18, 1e-10)

BULK_N_T_RANGE = (1e17, 1e23)
BULK_C_N_RANGE = (1e-15, 1e-11)
BULK_C_P_RANGE = (1e-15, 1e-11)
BULK_E_T_OFFSET = 0.1

EPS_R_OFFSET = 4.0
R_SHUNT_RANGE = (1e-4, 1e-1)
R_SERIES_RANGE = (1e-3, 1e1)

L_TCO_RANGE = (1e-8, 4e-7)
L_BE_RANGE = (1e-8, 1e-7)


def _common_randomization(
    layer: Layer,
    l_range: tuple[float, float],
    mobility_range: tuple[float, float],
) -> Layer:
    layer.L = loguniform(*l_range)
    layer.eps_r = randn() + EPS_R_OFFSET

    layer.N_c = loguniform(*N_C_RANGE)
    layer.N_D = loguniform(*N_D_A_RANGE)
    layer.N_A = loguniform(*N_D_A_RANGE)

    layer.mobilities.mu_n = loguniform(*mobility_range)
    layer.mobilities.mu_p = loguniform(*mobility_range)

    layer.interface.C_n_int = loguniform(*INTERFACE_C_N_P_RANGE)
    layer.interface.C_p_int = loguniform(*INTERFACE_C_N_P_RANGE)

    layer.bulk.N_t_bulk = loguniform(*BULK_N_T_RANGE)
    layer.bulk.C_n_bulk = loguniform(*BULK_C_N_RANGE)
    layer.bulk.C_p_bulk = loguniform(*BULK_C_P_RANGE)
    layer.bulk.E_t_bulk = uniform(
        layer.E_c + BULK_E_T_OFFSET,
        layer.E_v - BULK_E_T_OFFSET,
    )
    return layer


def randomize_layer1(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
) -> Layer:
    layer.E_c = randn() * 0.1 + layer2_E_c + ENERGY_LEVEL_OFFSET
    layer.E_v = randn() * 0.1 + layer2_E_v + ENERGY_LEVEL_OFFSET

    layer = _common_randomization(
        layer=layer,
        l_range=L_RANGE,
        mobility_range=MOBILITY_RANGE,
    )

    return layer


def randomize_layer2(
    layer: Layer,
    bandgap: float,
):
    layer.E_c = randn() * 0.5 + 3.0
    layer.E_v = layer.E_c + bandgap
    layer.mobilities.gamma_n = 0.0
    layer.generation.k_f = loguniform(*L2_KF_RANGE)
    layer.generation.k_direct = loguniform(*L2_K_DIRECT_RANGE)

    layer = _common_randomization(
        layer=layer,
        l_range=L2_L_RANGE,
        mobility_range=L2_MOBILITY_RANGE,
    )

    return layer


def randomize_layer3(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
):
    layer.E_c = randn() * 0.1 + layer2_E_c - ENERGY_LEVEL_OFFSET
    layer.E_v = randn() * 0.1 + layer2_E_v - ENERGY_LEVEL_OFFSET

    layer = _common_randomization(
        layer=layer,
        l_range=L_RANGE,
        mobility_range=MOBILITY_RANGE,
    )

    return layer


def randomize_simss_config(
    simss_config: SimssConfig,
    layer1_E_c: float,
    layer3_E_v: float,
) -> SimssConfig:
    w_l = randn() * 0.5 + layer1_E_c
    w_l = w_l if w_l > layer1_E_c else layer1_E_c
    w_r = randn() * 0.5 + layer3_E_v
    w_r = w_r if w_r < layer3_E_v else layer3_E_v

    simss_config.contacts.W_L = w_l
    simss_config.contacts.W_R = w_r
    simss_config.contacts.R_shunt = loguniform(*R_SHUNT_RANGE)
    simss_config.contacts.R_series = loguniform(*R_SERIES_RANGE)

    simss_config.optics.L_TCO = loguniform(*L_TCO_RANGE)
    simss_config.optics.L_BE = loguniform(*L_BE_RANGE)

    return simss_config
