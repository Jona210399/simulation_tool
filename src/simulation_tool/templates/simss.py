from dataclasses import asdict, dataclass

from simulation_tool.typing_ import PathLike, SettingsDict


@dataclass
class Layers:
    l1: str
    l2: str
    l3: str

    @classmethod
    def from_session_path(cls, session_path: PathLike) -> "Layers":
        return cls(
            f"{session_path}/L1_parameters.txt",
            f"{session_path}/L2_parameters.txt",
            f"{session_path}/L3_parameters.txt",
        )


@dataclass
class Contacts:
    W_L: float
    W_R: float
    S_n_L: float
    S_p_L: float
    S_n_R: float
    S_p_R: float
    R_shunt: float
    R_series: float

    @classmethod
    def get_default(cls) -> "Contacts":
        return cls(
            W_L=3.0,
            W_R=5.0,
            S_n_L=-1e-7,
            S_p_L=-1e-7,
            S_n_R=-1e-7,
            S_p_R=-1e-7,
            R_shunt=-1.0,
            R_series=0.0,
        )


@dataclass
class Optics:
    G_frac: float
    genProfile: str
    L_TCO: float
    L_BE: float
    nkSubstrate: str
    nkTCO: str
    nkBE: str
    spectrum: str
    lambda_min: float
    lambda_max: float

    @classmethod
    def get_default(cls, path_to_simss: PathLike) -> "Optics":
        return cls(
            G_frac=1.0,
            genProfile="calc",
            L_TCO=1.1e-7,
            L_BE=2e-7,
            nkSubstrate=f"{path_to_simss}/../Data/nk_SiO2.txt",
            nkTCO=f"{path_to_simss}/../Data/nk_ITO.txt",
            nkBE=f"{path_to_simss}/../Data/nk_Al.txt",
            spectrum=f"{path_to_simss}/../Data/AM15G.txt",
            lambda_min=3.5e-7,
            lambda_max=8e-7,
        )


@dataclass
class Numerical:
    NP: int
    tolPois: float
    maxDelV: float
    maxItPois: int
    maxItSS: int
    currDiffInt: int
    tolDens: float
    couplePC: int
    minAcc: float
    maxAcc: float
    ignoreNegDens: int
    failureMode: int
    grad: int

    @classmethod
    def get_default(cls) -> "Numerical":
        return cls(
            NP=400,
            tolPois=0.00001,
            maxDelV=10.0,
            maxItPois=1000,
            maxItSS=1000,
            currDiffInt=2,
            tolDens=1e-8,
            couplePC=4,
            minAcc=0.5,
            maxAcc=1.0,
            ignoreNegDens=1,
            failureMode=0,
            grad=1,
        )


@dataclass
class VoltageRange:
    Vdist: int
    preCond: int
    Vpre: float
    fixIons: int
    Vscan: int
    Vmin: float
    Vmax: float
    Vstep: float
    Vacc: int
    NJV: int
    untilVoc: int

    @classmethod
    def get_default(cls) -> "VoltageRange":
        return cls(
            Vdist=1,
            preCond=0,
            Vpre=0.0,
            fixIons=0,
            Vscan=1,
            Vmin=-0.5,
            Vmax=1.7,
            Vstep=0.025,
            Vacc=0,
            NJV=100,
            untilVoc=0,
        )


@dataclass
class UserInterface:
    timeout: int
    pauseAtEnd: int
    autoTidy: int
    useExpData: int
    expJV: str
    fitMode: str
    fitThreshold: float
    JVFile: str
    varFile: str
    limitDigits: int
    outputRatio: int
    scParsFile: str
    logFile: str

    @classmethod
    def get_default(cls) -> "UserInterface":
        return cls(
            timeout=-1,
            pauseAtEnd=0,
            autoTidy=1,
            useExpData=0,
            expJV="expJV.csv",
            fitMode="lin",
            fitThreshold=0.8,
            JVFile="JV.dat",
            varFile="Var.dat",
            limitDigits=1,
            outputRatio=1,
            scParsFile="scPars.dat",
            logFile="log.txt",
        )


@dataclass
class SimssConfig:
    T: float
    layers: Layers
    contacts: Contacts
    optics: Optics
    numerical: Numerical
    voltage: VoltageRange
    ui: UserInterface

    @classmethod
    def from_session(
        cls,
        session_path: PathLike,
        path_to_simss: PathLike,
    ) -> "SimssConfig":
        return cls(
            T=295.0,
            layers=Layers.from_session_path(session_path),
            contacts=Contacts.get_default(),
            optics=Optics.get_default(path_to_simss),
            numerical=Numerical.get_default(),
            voltage=VoltageRange.get_default(),
            ui=UserInterface.get_default(),
        )

    def to_dict(self) -> SettingsDict:
        dict_ = {
            "T": self.T,
            **asdict(self.layers),
            **asdict(self.contacts),
            **asdict(self.optics),
            **asdict(self.numerical),
            **asdict(self.voltage),
            **asdict(self.ui),
        }
        return dict_
