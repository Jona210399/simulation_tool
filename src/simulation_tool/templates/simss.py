from dataclasses import asdict, dataclass
from pathlib import Path

from simulation_tool.constants import SET_UP_FILE
from simulation_tool.templates.serializeable import JSONSerializable

ADDITIONAL_FILES_DIR = Path(__file__).parent.parent / "additional_files"


@dataclass
class LayerFiles:
    l1: Path
    l2: Path
    l3: Path

    @classmethod
    def from_session_path(cls, session_path: Path) -> "LayerFiles":
        return cls(
            session_path / "L1_parameters",
            session_path / "L2_parameters",
            session_path / "L3_parameters",
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

    def numeric_parameters(self) -> dict[str, float]:
        return asdict(self)


@dataclass
class Optics:
    G_frac: float
    genProfile: str
    L_TCO: float
    L_BE: float
    nkSubstrate: str
    nkTCO: str
    nkBE: str
    spectrum: Path
    lambda_min: float
    lambda_max: float

    @classmethod
    def get_default(cls) -> "Optics":
        return cls(
            G_frac=1.0,
            genProfile="calc",
            L_TCO=1.1e-7,
            L_BE=2e-7,
            nkSubstrate=str(ADDITIONAL_FILES_DIR / "nk_SiO2.txt"),
            nkTCO=str(ADDITIONAL_FILES_DIR / "nk_ITO.txt"),
            nkBE=str(ADDITIONAL_FILES_DIR / "nk_Al.txt"),
            spectrum=ADDITIONAL_FILES_DIR / "AM15G.txt",
            lambda_min=3.5e-7,
            lambda_max=8e-7,
        )

    def numeric_parameters(self) -> dict[str, float]:
        return {
            "G_frac": self.G_frac,
            "L_TCO": self.L_TCO,
            "L_BE": self.L_BE,
        }


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
        """The filenames are hardcoded in the simss code and should therefore not be changed."""
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
class SimssConfig(JSONSerializable):
    T: float
    setup_file: Path
    layers: LayerFiles
    contacts: Contacts
    optics: Optics
    numerical: Numerical
    voltage: VoltageRange
    ui: UserInterface

    @classmethod
    def from_session(
        cls,
        session_path: Path,
    ) -> "SimssConfig":
        return cls(
            T=295.0,
            setup_file=session_path / SET_UP_FILE,
            layers=LayerFiles.from_session_path(session_path),
            contacts=Contacts.get_default(),
            optics=Optics.get_default(),
            numerical=Numerical.get_default(),
            voltage=VoltageRange.get_default(),
            ui=UserInterface.get_default(),
        )

    def numeric_parameters(self) -> dict[str, float]:
        return {
            **self.contacts.numeric_parameters(),
            **self.optics.numeric_parameters(),
        }
