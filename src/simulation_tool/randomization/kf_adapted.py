import numpy as np

from simulation_tool.randomization.rng_funcs import RangeDict, RNGFunctions
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig

OPV = False
ENERGY_LEVEL_OFFSET = 0.1  # eV
L2_E_C_OFFSET = 3.0  # eV

N_C_RANGE = RangeDict(low=1e25, high=1e27)
N_D_A_RANGE = RangeDict(low=1e16, high=1e21)

L_RANGE = RangeDict(low=3e-9, high=3e-7)
L2_L_RANGE_OPV = RangeDict(low=5e-8, high=2e-7)
L2_L_RANGE_PEROVSKITE = RangeDict(low=3e-7, high=1e-6)
L2_L_RANGE = L2_L_RANGE_OPV if OPV else L2_L_RANGE_PEROVSKITE

MOBILITY_RANGE = RangeDict(low=1e-8, high=1e-3)  # (1e-9, 1e0)
L2_MOBILITY_RANGE_OPV = RangeDict(low=5e-9, high=1e-6)
L2_MOBILITY_RANGE_PEROVSKITE = RangeDict(low=5e-9, high=1e0)
L2_MOBILITY_RANGE = RangeDict(low=1e-8, high=1e-2)
MOBILITY_EXPONENT_OFFSET = 7.0
L2_MOBILITY_EXPONENT_OFFSET = 4.0

"""Interface trapping depends on the product of C_n/p and N_t_int. So only
one of them needs to be randomized. The other one can be set to 1.0. Interface trap energies need to be in between the energy levels of the interfacing layers."""
INTERFACE_C_N_P = 1.0
INTERFACE_N_INT_RANGE = RangeDict(low=1e-5, high=1e2)
INTERFACE_E_T_INT_OFFSET = 0.05

L2_KF_RANGE = RangeDict(low=1e4, high=1e8)
L2_K_DIRECT_RANGE = RangeDict(low=1e-18, high=1e-10)

BULK_N_T_RANGE = RangeDict(low=1e17, high=1e23)
BULK_C_N_RANGE = RangeDict(low=1e-15, high=1e-11)
BULK_C_P_RANGE = RangeDict(low=1e-15, high=1e-11)
BULK_E_T_OFFSET = 0.1

EPS_R_OFFSET = 4.0
R_SHUNT_RANGE = RangeDict(low=1e-4, high=1e-1)
R_SERIES_RANGE = RangeDict(low=1e-3, high=1e1)

L_TCO_RANGE = RangeDict(low=1e-8, high=4e-7)
L_BE_RANGE = RangeDict(low=1e-8, high=1e-7)


def randomize_device(
    layer1: Layer,
    layer2: Layer,
    layer3: Layer,
    simss_config: SimssConfig,
    bandgap: float,
    randomizer: RNGFunctions,
) -> tuple[tuple[Layer, Layer, Layer], SimssConfig]:
    layer2 = _randomize_layer2(
        layer=layer2,
        bandgap=bandgap,
        randomizer=randomizer,
    )
    layer1 = _randomize_layer1(
        layer=layer1,
        layer2_E_c=layer2.E_c,
        layer2_E_v=layer2.E_v,
        randomizer=randomizer,
    )
    layer3 = _randomize_layer3(
        layer=layer3,
        layer2_E_c=layer2.E_c,
        layer2_E_v=layer2.E_v,
        randomizer=randomizer,
    )

    layer2.interface.E_t_int = randomizer.uniform(
        "l2.E_t_int",
        low=layer2.E_c + INTERFACE_E_T_INT_OFFSET,
        high=layer3.E_v - INTERFACE_E_T_INT_OFFSET,
    )

    simss_config = _randomize_simss_config(
        simss_config=simss_config,
        layer1_E_c=layer1.E_c,
        layer3_E_v=layer3.E_v,
        randomizer=randomizer,
    )
    return (layer1, layer2, layer3), simss_config


def _calculate_layer_factor(l_range: RangeDict, layer_L: float) -> float:
    l_low, l_high = l_range["low"], l_range["high"]
    return (np.log(layer_L) - np.log(l_low)) / (np.log(l_high) - np.log(l_low))


def _common_randomization(
    layer: Layer,
    layer_prefix: str,
    l_range: RangeDict,
    mobility_exponent_offset: float,
    randomizer: RNGFunctions,
) -> Layer:
    layer.L = randomizer.loguniform(f"{layer_prefix}.L", **l_range)
    layer_factor = _calculate_layer_factor(l_range, layer.L)
    """Without the layer_factor, the convergence rate slightly dropped."""
    layer.eps_r = randomizer.randn(f"{layer_prefix}.eps_r") + EPS_R_OFFSET

    layer.N_c = randomizer.loguniform(f"{layer_prefix}.N_c", **N_C_RANGE)
    layer.N_D = randomizer.loguniform(f"{layer_prefix}.N_D", **N_D_A_RANGE)
    layer.N_A = randomizer.loguniform(f"{layer_prefix}.N_A", **N_D_A_RANGE)

    layer.mobilities.mu_n = 10.0 ** (
        -randomizer.randn(f"{layer_prefix}.mu_n") * 2.0
        - mobility_exponent_offset
        + layer_factor * 2.0
    )
    layer.mobilities.mu_p = 10.0 ** (
        -randomizer.randn(f"{layer_prefix}.mu_p") * 2.0
        - mobility_exponent_offset
        + layer_factor * 2.0
    )
    layer.interface.N_t_int = randomizer.loguniform(
        f"{layer_prefix}.N_t_int", **INTERFACE_N_INT_RANGE
    )
    layer.interface.C_n_int = INTERFACE_C_N_P
    layer.interface.C_p_int = INTERFACE_C_N_P
    layer.interface.intTrapType = randomizer.randint(
        f"{layer_prefix}.intTrapType",
        low=-1,
        high=2,
    )

    layer.bulk.N_t_bulk = randomizer.loguniform(
        f"{layer_prefix}.N_t_bulk", **BULK_N_T_RANGE
    )
    layer.bulk.C_n_bulk = randomizer.loguniform(
        f"{layer_prefix}.C_n_bulk", **BULK_C_N_RANGE
    )
    layer.bulk.C_p_bulk = randomizer.loguniform(
        f"{layer_prefix}.C_p_bulk", **BULK_C_P_RANGE
    )
    layer.bulk.E_t_bulk = randomizer.uniform(
        f"{layer_prefix}.E_t_bulk",
        low=layer.E_c + BULK_E_T_OFFSET,
        high=layer.E_v - BULK_E_T_OFFSET,
    )
    return layer


def _randomize_layer1(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
    randomizer: RNGFunctions,
) -> Layer:
    layer_prefix = "l1"

    layer.E_c = (
        randomizer.randn(f"{layer_prefix}.E_c") * 0.1 + layer2_E_c + ENERGY_LEVEL_OFFSET
    )
    layer.E_v = (
        randomizer.randn(f"{layer_prefix}.E_v") * 0.1 + layer2_E_v + ENERGY_LEVEL_OFFSET
    )
    layer.interface.E_t_int = randomizer.uniform(
        f"{layer_prefix}.E_t_int",
        low=layer.E_c + INTERFACE_E_T_INT_OFFSET,
        high=layer2_E_v - INTERFACE_E_T_INT_OFFSET,
    )

    layer = _common_randomization(
        layer=layer,
        layer_prefix=layer_prefix,
        l_range=L_RANGE,
        mobility_exponent_offset=MOBILITY_EXPONENT_OFFSET,
        randomizer=randomizer,
    )

    return layer


def _randomize_layer2(
    layer: Layer,
    bandgap: float,
    randomizer: RNGFunctions,
):
    layer_prefix = "l2"
    layer.E_c = randomizer.randn(f"{layer_prefix}.E_c") * 0.5 + L2_E_C_OFFSET
    layer.E_v = layer.E_c + bandgap
    layer.mobilities.gamma_n = 0.0
    layer.generation.k_f = randomizer.loguniform(f"{layer_prefix}.k_f", **L2_KF_RANGE)
    layer.generation.k_direct = randomizer.loguniform(
        f"{layer_prefix}.k_direct", **L2_K_DIRECT_RANGE
    )

    layer = _common_randomization(
        layer=layer,
        layer_prefix=layer_prefix,
        l_range=L2_L_RANGE,
        mobility_exponent_offset=L2_MOBILITY_EXPONENT_OFFSET,
        randomizer=randomizer,
    )

    return layer


def _randomize_layer3(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
    randomizer: RNGFunctions,
):
    layer_prefix = "l3"
    layer.E_c = (
        randomizer.randn(f"{layer_prefix}.E_c") * 0.1 + layer2_E_c - ENERGY_LEVEL_OFFSET
    )
    layer.E_v = (
        randomizer.randn(f"{layer_prefix}.E_v") * 0.1 + layer2_E_v - ENERGY_LEVEL_OFFSET
    )

    layer = _common_randomization(
        layer=layer,
        layer_prefix=layer_prefix,
        l_range=L_RANGE,
        mobility_exponent_offset=MOBILITY_EXPONENT_OFFSET,
        randomizer=randomizer,
    )
    layer.interface.N_t_int = 0.0  # disables trapping at interface to electrode

    return layer


def _randomize_simss_config(
    simss_config: SimssConfig,
    layer1_E_c: float,
    layer3_E_v: float,
    randomizer: RNGFunctions,
) -> SimssConfig:
    w_l = randomizer.randn("W_L") * 0.5 + layer1_E_c
    w_l = w_l if w_l > layer1_E_c else layer1_E_c
    w_r = randomizer.randn("W_R") * 0.5 + layer3_E_v
    w_r = w_r if w_r < layer3_E_v else layer3_E_v

    simss_config.contacts.W_L = w_l
    simss_config.contacts.W_R = w_r
    # simss_config.contacts.R_shunt = randomizer.loguniform(*R_SHUNT_RANGE)
    # simss_config.contacts.R_series = randomizer.loguniform(*R_SERIES_RANGE)
    """Randomization of the Resistances leads to significantly lower simulation convergence rates."""

    simss_config.optics.L_TCO = randomizer.loguniform("L_TCO", **L_TCO_RANGE)
    simss_config.optics.L_BE = randomizer.loguniform("L_BE", **L_BE_RANGE)

    simss_config.numerical.maxItSS *= 2.0

    return simss_config
