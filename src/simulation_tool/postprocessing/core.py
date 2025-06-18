from pathlib import Path

import polars as pl

from simulation_tool.constants import (
    EQE_SIM_OUTPUT_FILE_NAME,
    JV_SIM_OUTPUT_FILE_NAME,
    NK_FILE_NAME,
    PARAMS_OUTPUT_FILE_NAME,
    SET_UP_FILE,
    UVVIS_FILE_NAME,
)
from simulation_tool.data.optical import OpticalData
from simulation_tool.data.uvvis import UVVisData
from simulation_tool.EQE.data import EQEData
from simulation_tool.jV.data import JVData
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig, UserInterface
from simulation_tool.utils import add_prefix_to_keys


def _read_session_configuration(session_path: Path) -> tuple[SimssConfig, list[Layer]]:
    simss_config: SimssConfig = SimssConfig.from_json_file(
        (session_path / SET_UP_FILE).with_suffix(".json")
    )  # type: ignore

    layer_files = [
        simss_config.layers.l1,
        simss_config.layers.l2,
        simss_config.layers.l3,
    ]

    layers: list[Layer] = [
        Layer.from_json_file(layer_file.with_suffix(".json"))
        for layer_file in layer_files
    ]  # type: ignore

    return simss_config, layers


def _read_simulation_data(
    session_path: Path,
    simss_ui: UserInterface,
) -> tuple[JVData, UVVisData, EQEData, OpticalData]:
    data_dir = session_path / "data"

    jv_data = JVData.from_files(
        device_characteristics_file=(data_dir / simss_ui.scParsFile).with_suffix(
            ".txt"
        ),
        jv_file=(data_dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"),
    )

    uvvis_data = UVVisData.from_file(
        (data_dir / UVVIS_FILE_NAME).with_suffix(".parquet")
    )
    eqe_data = EQEData.from_file(
        (data_dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet")
    )

    optical_data = OpticalData.from_file(
        (data_dir / NK_FILE_NAME).with_suffix(".parquet")
    )

    return jv_data, uvvis_data, eqe_data, optical_data


def _get_relevant_parameters(
    simss_config: SimssConfig, layers: list[Layer]
) -> dict[str, float]:
    parameters = SimssConfig.numeric_parameters(simss_config)

    for i, layer in enumerate(layers, start=1):
        parameters.update(add_prefix_to_keys(layer.numeric_parameters(), f"l{i}."))

    return parameters


def post_process_session(session_path: Path):
    simss_config, layers = _read_session_configuration(session_path)
    jv_data, uvvis_data, eqe_data, optical_data = _read_simulation_data(
        session_path, simss_config.ui
    )

    parameters = _get_relevant_parameters(simss_config, layers)

    return jv_data, uvvis_data, eqe_data, optical_data, parameters


def read_post_processed_data(
    dir: Path,
) -> tuple[
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
    pl.DataFrame,
]:
    jVs = pl.read_parquet((dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    uvvis = pl.read_parquet((dir / UVVIS_FILE_NAME).with_suffix(".parquet"))
    eqes = pl.read_parquet((dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    optical = pl.read_parquet((dir / NK_FILE_NAME).with_suffix(".parquet"))
    parameters = pl.read_parquet(
        (dir / PARAMS_OUTPUT_FILE_NAME).with_suffix(".parquet")
    )

    return jVs, uvvis, eqes, optical, parameters
