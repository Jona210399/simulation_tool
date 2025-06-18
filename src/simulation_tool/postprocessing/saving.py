from pathlib import Path

import numpy as np
import polars as pl

from simulation_tool.constants import (
    EQE_SIM_OUTPUT_FILE_NAME,
    JV_SIM_OUTPUT_FILE_NAME,
    NK_FILE_NAME,
    PARAMS_OUTPUT_FILE_NAME,
    UVVIS_FILE_NAME,
)
from simulation_tool.data import OpticalData, UVVisData
from simulation_tool.EQE.data import EQEData
from simulation_tool.jV.data import JVData
from simulation_tool.typing_ import Array1D, Array2D


def save_jV_data(
    jVs: list[JVData], save_dir: Path
) -> tuple[
    pl.DataFrame,
    tuple[
        Array2D,
        Array2D,
        Array1D,
        Array1D,
        Array1D,
        Array1D,
    ],
]:
    v, j, ff, voc, jsc, pce = combined = JVData.combine(jVs)

    df = pl.concat(
        [
            pl.from_numpy(v, schema=["Vext"]),
            pl.from_numpy(j, schema=["Jext"]),
            pl.from_numpy(ff, schema=["FF"]),
            pl.from_numpy(voc, schema=["Voc"]),
            pl.from_numpy(jsc, schema=["Jsc"]),
            pl.from_numpy(pce, schema=["PCE"]),
        ],
        how="horizontal",
    )

    df.write_parquet(file=(save_dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))

    return df, combined


def save_uvvis_data(
    uvvis: list[UVVisData], save_dir: Path
) -> tuple[
    pl.DataFrame,
    tuple[
        Array2D,
        Array2D,
        Array2D,
        Array2D,
    ],
]:
    a, r, _, wavelenghts = combined = UVVisData.combine(uvvis)

    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(a, schema=["absorption"]),
            pl.from_numpy(r, schema=["reflection"]),
        ],
        how="horizontal",
    )

    df.write_parquet(file=(save_dir / UVVIS_FILE_NAME).with_suffix(".parquet"))

    return df, combined


def save_optical_data(
    optical: list[OpticalData], save_dir: Path
) -> tuple[pl.DataFrame, tuple[Array2D, Array2D, Array2D]]:
    n, k, wavelenghts = combined = OpticalData.combine(optical)
    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(n, schema=["n"]),
            pl.from_numpy(k, schema=["k"]),
        ],
        how="horizontal",
    )

    df.write_parquet(file=(save_dir / NK_FILE_NAME).with_suffix(".parquet"))

    return df, combined


def save_eqe_data(
    eqes: list[EQEData], save_dir: Path
) -> tuple[pl.DataFrame, tuple[Array2D, Array2D]]:
    eqe, wavelenghts = combined = EQEData.combine(eqes)

    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(eqe, schema=["EQE"]),
        ],
        how="horizontal",
    )
    df.write_parquet(file=(save_dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))

    return df, combined


def save_parameters(parameters: list[dict], save_dir: Path) -> pl.DataFrame:
    data = np.stack([np.array(list(params.values())) for params in parameters])
    columns = list(parameters[0].keys())
    schema = {col: pl.Float32 for col in columns}
    df = pl.DataFrame(
        data=data,
        schema=schema,
    )
    df.write_parquet(file=(save_dir / PARAMS_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    return df
