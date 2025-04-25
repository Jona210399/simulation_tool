from pathlib import Path

import numpy as np
from numpy.typing import NDArray

PathLike = str | Path
Array1D = NDArray[np.float64]
SettingsDict = dict[str, str | int | float]
