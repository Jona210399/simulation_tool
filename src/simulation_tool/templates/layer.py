from dataclasses import dataclass

from simulation_tool.templates.serializeable import JSONSerializable
from simulation_tool.typing_ import PathLike


@dataclass
class Mobilities:
    mu_n: float
    mu_p: float
    mobnDep: int
    mobpDep: int
    gamma_n: float
    gamma_p: float

    @classmethod
    def get_default(cls) -> "Mobilities":
        return cls(
            mu_n=None,
            mu_p=None,
            mobnDep=0,
            mobpDep=0,
            gamma_n=None,
            gamma_p=0.0,
        )


@dataclass
class Interface:
    nu_int_n: float
    nu_int_p: float
    N_t_int: float
    E_t_int: float
    intTrapFile: str
    intTrapType: int
    C_n_int: float
    C_p_int: float

    @classmethod
    def get_default(cls) -> "Interface":
        return cls(
            nu_int_n=1000.0,
            nu_int_p=1000.0,
            N_t_int=0.0,
            E_t_int=4.0,
            intTrapFile="none",
            intTrapType=None,
            C_n_int=1e-13,
            C_p_int=1e-13,
        )


@dataclass
class Ions:
    N_anion: float
    N_cation: float
    mu_anion: float
    mu_cation: float
    ionsMayEnter: int

    @classmethod
    def get_default(cls) -> "Ions":
        return cls(
            N_anion=0.0,
            N_cation=0.0,
            mu_anion=1e-14,
            mu_cation=1e-14,
            ionsMayEnter=None,
        )


@dataclass
class Generation:
    G_ehp: float
    layerGen: int
    nkLayer: str
    fieldDepG: int
    P0: float
    a: float
    thermLengDist: int
    k_f: float
    k_direct: float
    preLangevin: float
    useLangevin: int

    @classmethod
    def get_default(cls) -> "Generation":
        return cls(
            G_ehp=None,
            layerGen=None,
            nkLayer=None,
            fieldDepG=0,
            P0=0.0,
            a=1e-9,
            thermLengDist=2,
            k_f=1e6,
            k_direct=1e-18,
            preLangevin=1.0,
            useLangevin=1,
        )


@dataclass
class BulkTrapping:
    N_t_bulk: float
    C_n_bulk: float
    C_p_bulk: float
    E_t_bulk: float
    bulkTrapFile: str
    bulkTrapType: int

    @classmethod
    def get_default(cls) -> "BulkTrapping":
        return cls(
            N_t_bulk=None,
            C_n_bulk=1e-13,
            C_p_bulk=1e-13,
            E_t_bulk=4.0,
            bulkTrapFile="none",
            bulkTrapType=-1,
        )


@dataclass
class Layer(JSONSerializable):
    L: float
    eps_r: float
    E_c: float
    E_v: float
    N_c: float
    N_D: float
    N_A: float
    mobilities: Mobilities
    interface: Interface
    ions: Ions
    generation: Generation
    bulk: BulkTrapping

    @classmethod
    def get_default(cls) -> "Layer":
        return cls(
            L=None,
            eps_r=None,
            E_c=None,
            E_v=None,
            N_c=None,
            N_D=0.0,
            N_A=0.0,
            mobilities=Mobilities.get_default(),
            interface=Interface.get_default(),
            ions=Ions.get_default(),
            generation=Generation.get_default(),
            bulk=BulkTrapping.get_default(),
        )

    @classmethod
    def get_default_layer1(cls, path_to_simss: PathLike) -> "Layer":
        layer = cls.get_default()

        layer.L = 3e-8
        layer.eps_r = 5.0
        layer.E_c = 3.0
        layer.E_v = 5.5
        layer.N_c = 1e26

        layer.mobilities.mu_n = 1e-6
        layer.mobilities.mu_p = 1e-6
        layer.mobilities.gamma_n = 0.0

        layer.interface.intTrapType = -1

        layer.ions.ionsMayEnter = 0

        layer.generation.G_ehp = 0.0
        layer.generation.layerGen = 0
        layer.generation.nkLayer = f"{path_to_simss}/../Data/nk_PEDOT.txt"
        layer.bulk.N_t_bulk = 0.0

        return layer

    @classmethod
    def get_default_layer2(cls, session_path: PathLike, bandgap: float) -> "Layer":
        layer = cls.get_default()

        layer.L = 5.4e-8
        layer.eps_r = 4.0
        layer.E_c = 3.0
        layer.E_v = 3.0 + bandgap
        layer.N_c = 1e25

        layer.mobilities.mu_n = 0.01
        layer.mobilities.mu_p = 0.01
        layer.mobilities.gamma_n = 0.0001

        layer.interface.intTrapType = 1

        layer.ions.ionsMayEnter = 1

        layer.generation.G_ehp = 7.25e27
        layer.generation.layerGen = 1
        layer.generation.nkLayer = f"{session_path}/nk.txt"

        layer.bulk.N_t_bulk = 4.42e21

        return layer

    @classmethod
    def get_default_layer3(cls, path_to_simss: PathLike) -> "Layer":
        layer = cls.get_default()

        layer.L = 3e-8
        layer.eps_r = 3.5
        layer.E_c = 2.5
        layer.E_v = 5.0
        layer.N_c = 1e26

        layer.mobilities.mu_n = 1e-8
        layer.mobilities.mu_p = 1e-8
        layer.mobilities.gamma_n = 0.0

        layer.interface.intTrapType = 1

        layer.ions.ionsMayEnter = 0

        layer.generation.G_ehp = 0.0
        layer.generation.layerGen = 0
        layer.generation.nkLayer = f"{path_to_simss}/../Data/nk_Ca.txt"

        layer.bulk.N_t_bulk = 0.0

        return layer
