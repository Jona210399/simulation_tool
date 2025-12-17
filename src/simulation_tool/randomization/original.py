import numpy as np

from simulation_tool.randomization.rng_funcs import RNGFunctions
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig

ENERGY_LEVEL_OFFSET = 0.1  # eV


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
    simss_config = _randomize_simss_config(
        simss_config=simss_config,
        layer1_E_c=layer1.E_c,
        layer3_E_v=layer3.E_v,
        randomizer=randomizer,
    )
    return (layer1, layer2, layer3), simss_config


def _randomize_layer1(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
    randomizer: RNGFunctions,
) -> Layer:
    layer.L = randomizer.loguniform(low=3e-9, high=3e-7)
    layer_factor = (np.log(layer.L) - np.log(3e-9)) / (np.log(3e-7) - np.log(3e-9))

    layer.eps_r = randomizer.randn() * 1.0 + 4.0

    layer.E_c = randomizer.randn() * 0.1 + layer2_E_c + ENERGY_LEVEL_OFFSET
    layer.E_v = randomizer.randn() * 0.1 + layer2_E_v + ENERGY_LEVEL_OFFSET

    layer.N_c = randomizer.loguniform(low=1e16, high=1e27)
    layer.N_D = randomizer.loguniform(low=1e16, high=1e21)
    layer.N_A = randomizer.loguniform(low=1e16, high=1e21)

    layer.mobilities.mu_n = 10 ** (-randomizer.randn() * 2.0 - 7.0 + layer_factor * 2.0)
    layer.mobilities.mu_p = 10 ** (-randomizer.randn() * 2.0 - 7.0 + layer_factor * 2.0)

    layer.interface.C_n_int = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.interface.C_p_int = randomizer.loguniform(low=1e-15, high=1e-11)

    layer.bulk.N_t_bulk = randomizer.loguniform(low=1e17, high=1e23)
    layer.bulk.C_n_bulk = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.bulk.C_p_bulk = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.bulk.E_t_bulk = randomizer.uniform(low=layer.E_c + 0.1, high=layer.E_v - 0.1)

    return layer


def _randomize_layer2(
    layer: Layer,
    bandgap: float,
    randomizer: RNGFunctions,
):
    layer.L = randomizer.loguniform(low=5e-9, high=7.5e-7)
    layer_factor = (np.log(layer.L) - np.log(5e-9)) / (np.log(7.5e-7) - np.log(5e-9))
    layer.eps_r = randomizer.randn() * 1.0 + 4.0
    layer.E_c = randomizer.randn() * 0.5 + 3.0
    layer.E_v = layer.E_c + bandgap

    layer.N_c = randomizer.loguniform(low=1e24, high=1e27)
    layer.N_D = randomizer.loguniform(low=1e16, high=1e21)
    layer.N_A = randomizer.loguniform(low=1e16, high=1e21)

    layer.mobilities.mu_n = 10.0 ** (-randomizer.randn() * 2 - 4.0 + layer_factor * 2.0)
    layer.mobilities.mu_p = 10.0 ** (
        -randomizer.randn() * 2.0 - 4.0 + layer_factor * 2.0
    )

    layer.mobilities.gamma_n = 0.0

    layer.interface.C_n_int = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.interface.C_p_int = randomizer.loguniform(low=1e-15, high=1e-11)

    layer.generation.a = randomizer.loguniform(low=0.5e-9, high=1e-8)
    layer.generation.P0 = np.clip(randomizer.randn() * 0.2, a_min=0.0, a_max=1.0)
    layer.generation.fieldDepG = 1 if layer.generation.P0 > 0.0 else 0
    layer.generation.k_f = randomizer.loguniform(low=1e4, high=1e8)
    layer.generation.k_direct = randomizer.loguniform(low=1e-19, high=1e-17)

    layer.bulk.N_t_bulk = randomizer.loguniform(low=1e17, high=1e23)
    layer.bulk.C_n_bulk = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.bulk.C_p_bulk = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.bulk.E_t_bulk = randomizer.uniform(low=layer.E_c + 0.1, high=layer.E_v - 0.1)

    return layer


def _randomize_layer3(
    layer: Layer,
    layer2_E_c: float,
    layer2_E_v: float,
    randomizer: RNGFunctions,
):
    layer.L = randomizer.loguniform(low=3e-9, high=3e-7)
    layer_factor = (np.log(layer.L) - np.log(3e-9)) / (np.log(3e-7) - np.log(3e-9))

    layer.eps_r = randomizer.randn() * 1.0 + 4.0
    layer.E_c = randomizer.randn() * 0.1 + layer2_E_c - ENERGY_LEVEL_OFFSET
    layer.E_v = randomizer.randn() * 0.1 + layer2_E_v - ENERGY_LEVEL_OFFSET

    layer.N_c = randomizer.loguniform(low=1e24, high=1e27)
    layer.N_D = randomizer.loguniform(low=1e16, high=1e21)
    layer.N_A = randomizer.loguniform(low=1e16, high=1e21)

    layer.mobilities.mu_n = 10.0 ** (
        -randomizer.randn() * 2.0 - 7.0 + layer_factor * 2.0
    )
    layer.mobilities.mu_p = 10.0 ** (
        -randomizer.randn() * 2.0 - 7.0 + layer_factor * 2.0
    )

    layer.interface.C_n_int = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.interface.C_p_int = randomizer.loguniform(low=1e-15, high=1e-11)

    layer.bulk.N_t_bulk = randomizer.loguniform(low=1e17, high=1e23)
    layer.bulk.C_n_bulk = randomizer.loguniform(low=1e-15, high=1e-11)
    layer.bulk.C_p_bulk = randomizer.loguniform(low=1e-15, high=1e-11)

    layer.bulk.E_t_bulk = randomizer.uniform(
        low=layer.E_c + 0.01, high=layer.E_v - 0.01
    )

    return layer


def _randomize_simss_config(
    simss_config: SimssConfig,
    layer1_E_c: float,
    layer3_E_v: float,
    randomizer: RNGFunctions,
) -> SimssConfig:
    W_L = randomizer.randn() * 0.5 + layer1_E_c
    W_L = W_L if W_L > layer1_E_c else layer1_E_c
    W_R = randomizer.randn() * 0.5 + layer3_E_v
    W_R = W_R if W_R < layer3_E_v else layer3_E_v

    simss_config.contacts.W_L = W_L
    simss_config.contacts.W_R = W_R
    simss_config.optics.L_TCO = randomizer.loguniform(low=1e-8, high=3e-7)
    simss_config.optics.L_BE = randomizer.loguniform(low=1e-8, high=3e-7)

    return simss_config
