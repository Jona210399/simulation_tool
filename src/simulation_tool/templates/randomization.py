import numpy as np

from simulation_tool.templates.layer import Layer
from simulation_tool.typing_ import PathLike
from simulation_tool.utils import loguniform


def get_random_layer1(
    path_to_simss: PathLike,
    layer2_E_c: float,
    layer2_E_v: float,
) -> Layer:
    layer = Layer.get_default_layer1(path_to_simss)

    layer.L = loguniform(low=3e-9, high=3e-7)
    layer_factor = (np.log(layer.L) - np.log(3e-9)) / (np.log(3e-7) - np.log(3e-9))

    layer.eps_r = np.random.randn() * 1.0 + 4.0

    layer.E_c = np.random.randn() * 0.1 + layer2_E_c + 0.1
    layer.E_v = np.random.randn() * 0.1 + layer2_E_v + 0.1

    layer.N_c = loguniform(low=1e16, high=1e27)
    layer.N_D = loguniform(low=1e16, high=1e21)
    layer.N_A = loguniform(low=1e16, high=1e21)

    layer.mobilities.mu_n = 10 ** (-np.random.randn() * 2.0 - 7.0 + layer_factor * 2.0)
    layer.mobilities.mu_p = 10 ** (-np.random.randn() * 2.0 - 7.0 + layer_factor * 2.0)

    layer.interface.C_n_int = loguniform(low=1e-15, high=1e-11)
    layer.interface.C_p_int = loguniform(low=1e-15, high=1e-11)

    layer.bulk.N_t_bulk = loguniform(low=1e17, high=1e23)
    layer.bulk.C_n_bulk = loguniform(low=1e-15, high=1e-11)
    layer.bulk.C_p_bulk = loguniform(low=1e-15, high=1e-11)
    layer.bulk.E_t_bulk = np.random.uniform(layer.E_c + 0.1, layer.E_v - 0.1)

    return layer


def get_random_layer2(
    session_path: PathLike,
    bandgap: float,
):
    layer = Layer.get_default_layer2(session_path, bandgap)

    layer.L = loguniform(low=5e-9, high=7.5e-7)
    layer_factor = (np.log(layer.L) - np.log(5e-9)) / (np.log(7.5e-7) - np.log(5e-9))
    layer.eps_r = np.random.randn() * 1.0 + 4.0
    layer.E_c = np.random.randn() * 0.5 + 3.0
    layer.E_v = layer.E_c + bandgap

    layer.N_c = loguniform(low=1e24, high=1e27)
    layer.N_D = loguniform(low=1e16, high=1e21)
    layer.N_A = loguniform(low=1e16, high=1e21)

    layer.mobilities.mu_n = 10.0 ** (-np.random.randn() * 2 - 4.0 + layer_factor * 2.0)
    layer.mobilities.mu_p = 10.0 ** (
        -np.random.randn() * 2.0 - 4.0 + layer_factor * 2.0
    )

    layer.mobilities.gamma_n = 0.0

    layer.interface.C_n_int = loguniform(low=1e-15, high=1e-11)
    layer.interface.C_p_int = loguniform(low=1e-15, high=1e-11)

    layer.generation.a = loguniform(low=0.5e-9, high=1e-8)
    layer.generation.P0 = np.clip(np.random.randn() * 0.2, a_min=0.0, a_max=1.0)
    layer.generation.fieldDepG = 1 if layer.generation.P0 > 0.0 else 0
    layer.generation.k_f = loguniform(low=1e4, high=1e8)
    layer.generation.k_direct = loguniform(low=1e-19, high=1e-17)

    layer.bulk.N_t_bulk = loguniform(low=1e17, high=1e23)
    layer.bulk.C_n_bulk = loguniform(low=1e-15, high=1e-11)
    layer.bulk.C_p_bulk = loguniform(low=1e-15, high=1e-11)
    layer.bulk.E_t_bulk = np.random.uniform(layer.E_c + 0.1, layer.E_v - 0.1)

    return layer
