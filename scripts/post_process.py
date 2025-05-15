import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib import pyplot as plt
from numpy.typing import NDArray

from simulation_tool.constants import (
    EQE_SIM_OUTPUT_FILE_NAME,
    JV_SIM_OUTPUT_FILE_NAME,
    NK_FILE_NAME,
    RUN_DIR_PREFIX,
    SET_UP_FILE,
    UVVIS_FILE_NAME,
)
from simulation_tool.data import OpticalData, UVVisData
from simulation_tool.EQE import EQEData
from simulation_tool.jV import JVData
from simulation_tool.templates.layer import Layer
from simulation_tool.templates.simss import SimssConfig, UserInterface
from simulation_tool.utils import add_prefix_to_keys, get_sorted_run_dirs

LINE_WIDTH: float = 1.0
LINE_ALPHA: float = 0.5


def read_session_configuration(session_path: Path) -> tuple[SimssConfig, list[Layer]]:
    simss_config = SimssConfig.from_json_file(
        (session_path / SET_UP_FILE).with_suffix(".json")
    )

    layer_files = [
        simss_config.layers.l1,
        simss_config.layers.l2,
        simss_config.layers.l3,
    ]

    layers = [
        Layer.from_json_file(layer_file.with_suffix(".json"))
        for layer_file in layer_files
    ]

    return simss_config, layers


def read_simulation_data(
    session_path: Path,
    simss_ui: UserInterface,
) -> tuple[JVData, UVVisData, EQEData, OpticalData]:
    data_dir = session_path / "data"

    jv_data = JVData.from_files(
        device_characteristics_file=data_dir / simss_ui.scParsFile,
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


def get_relevant_parameters(
    simss_config: SimssConfig, layers: list[Layer]
) -> dict[str, float]:
    parameters = SimssConfig.numeric_parameters(simss_config)

    for i, layer in enumerate(layers, start=1):
        parameters.update(add_prefix_to_keys(layer.numeric_parameters(), f"l{i}."))

    return parameters


def post_process_session(session_path: Path):
    simss_config, layers = read_session_configuration(session_path)
    jv_data, uvvis_data, eqe_data, optical_data = read_simulation_data(
        session_path, simss_config.ui
    )

    parameters = get_relevant_parameters(simss_config, layers)

    return jv_data, uvvis_data, eqe_data, optical_data, parameters


def plot_device_characteristics_distributions(
    pces: NDArray,
    jscs: NDArray,
    vocs: NDArray,
    ffs: NDArray,
    save_path: Path,
    dpi: int,
):
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)

    plt.hist(pces, bins=20, color="C0", edgecolor="black")
    plt.xlabel("PCE [%]")
    plt.ylabel("Counts")
    plt.title("PCE distribution")

    plt.subplot(2, 2, 2)
    plt.hist(jscs, bins=20, color="C1", edgecolor="black")
    plt.xlabel("Jsc [A/m²]")
    plt.ylabel("Counts")
    plt.title("Jsc distribution")

    plt.subplot(2, 2, 3)
    plt.hist(vocs, bins=20, color="C2", edgecolor="black")
    plt.xlabel("Voc [V]")
    plt.ylabel("Counts")
    plt.title("Voc distribution")

    plt.subplot(2, 2, 4)
    plt.hist(ffs * 100, bins=20, color="C3", edgecolor="black")
    plt.xlabel("Fill Factor [%]")
    plt.ylabel("Counts")
    plt.title("FF distribution")

    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_jVs(
    jVs: list[JVData],
    save_path: Path,
    dpi: int,
):
    _, ax = plt.subplots(figsize=(8, 6))

    for jV in jVs:
        jV.plot(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA, label=False)

    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_EQEs(eqes: list[EQEData], save_path: Path, dpi: int = 300):
    _, ax = plt.subplots(figsize=(8, 6))

    for eqe in eqes:
        eqe.plot(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_path, dpi=dpi)
    plt.close()


def plot_uvvis(
    uvvis: list[UVVisData],
    save_dir: Path,
    dpi: int,
):
    _, ax = plt.subplots(figsize=(8, 6))

    for uvv in uvvis:
        uvv.plot_absorption(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_dir / "absorption.png", dpi=dpi)
    plt.close()

    _, ax = plt.subplots(figsize=(8, 6))
    for uvv in uvvis:
        uvv.plot_reflection(ax=ax, linewidth=LINE_WIDTH, alpha=LINE_ALPHA)

    plt.savefig(save_dir / "reflection.png", dpi=dpi)
    plt.close()


def save_jV_data(
    jVs: list[JVData], save_dir: Path
) -> tuple[
    NDArray,
    NDArray,
    NDArray,
    NDArray,
    NDArray,
    NDArray,
]:
    v = np.stack([jv.Vext for jv in jVs])
    j = np.stack([jv.Jext for jv in jVs])
    ff = np.array([jv.ff for jv in jVs])
    voc = np.array([jv.voc for jv in jVs])
    jsc = np.array([jv.jsc for jv in jVs])
    pce = np.array([jv.calculate_pce() for jv in jVs])

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

    return v, j, ff, voc, jsc, pce


def save_uvvis_data(
    uvvis: list[UVVisData], save_dir: Path
) -> tuple[
    NDArray,
    NDArray,
    NDArray,
]:
    wavelenghts = np.stack((uvv.wavelengths for uvv in uvvis))
    a = np.stack((uvv.absorption for uvv in uvvis))
    r = np.stack((uvv.reflection for uvv in uvvis))

    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(a, schema=["absorption"]),
            pl.from_numpy(r, schema=["reflection"]),
        ],
        how="horizontal",
    )

    df.write_parquet(file=(save_dir / UVVIS_FILE_NAME).with_suffix(".parquet"))

    return a, r, wavelenghts


def save_optical_data(
    optical: list[OpticalData], save_dir: Path
) -> tuple[
    NDArray,
    NDArray,
    NDArray,
]:
    n = np.stack([opt.n for opt in optical])
    k = np.stack([opt.k for opt in optical])
    wavelenghts = np.stack([opt.wavelenghts for opt in optical])

    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(n, schema=["n"]),
            pl.from_numpy(k, schema=["k"]),
        ],
        how="horizontal",
    )

    df.write_parquet(file=(save_dir / NK_FILE_NAME).with_suffix(".parquet"))

    return n, k, wavelenghts


def save_eqe_data(eqes: list[EQEData], save_dir: Path) -> tuple[NDArray, NDArray]:
    eqe = np.stack([eqe.EQE for eqe in eqes])
    wavelenghts = np.stack([eqe.wavelenghts for eqe in eqes])
    df = pl.concat(
        [
            pl.from_numpy(wavelenghts, schema=["wavelenghts"]),
            pl.from_numpy(eqe, schema=["EQE"]),
        ],
        how="horizontal",
    )
    df.write_parquet(file=(save_dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))

    return eqe, wavelenghts


def save_parameters(parameters: list[dict], save_dir: Path) -> pl.DataFrame:
    data = np.stack([np.array(list(params.values())) for params in parameters])
    columns = list(parameters[0].keys())
    df = pl.DataFrame(
        data=data,
        schema=columns,
    )
    df.write_parquet(file=save_dir / "params.parquet")
    return df


def find_devices_with_best_pce(
    pces: NDArray,
    run_dirs: list[Path],
    save_dir: Path,
    top_k: int,
    jvs: list[JVData],
) -> list[int]:
    run_ids: list[int] = np.argsort(pces)[-top_k:][::-1].tolist()

    save_dir = save_dir / "best_pces"
    save_dir.mkdir(exist_ok=True)

    d = defaultdict(list)

    for run_id in run_ids:
        run = run_dirs[run_id]
        jv = jvs[run_id]
        pce = pces[run_id]
        d["run"].append(run.name)
        d["ff"].append(jv.ff)
        d["voc"].append(jv.voc)
        d["jsc"].append(jv.jsc)
        d["pce"].append(pce)
        shutil.copytree(run, save_dir / run.name)

    pl.DataFrame(data=d).write_csv(file=save_dir / "best_pces.csv")

    return run_ids


def post_process_simulations(simulation_dir: Path):
    MAX_CURVES_TO_PLOT: int | None = 50
    TOP_K: int = 5
    DPI: int = 300
    SAVE_DIR = simulation_dir / "post_processed"
    data_dir = SAVE_DIR / "data"
    plot_dir = SAVE_DIR / "plots"
    data_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    jVs: list[JVData] = []
    uvvis: list[UVVisData] = []
    eqes: list[EQEData] = []
    optical: list[OpticalData] = []
    parameters = []

    run_dirs = get_sorted_run_dirs(simulation_dir, run_dir_prefix=RUN_DIR_PREFIX)

    for session_path in run_dirs:
        jv_data, uvvis_data, eqe_data, optical_data, params = post_process_session(
            session_path
        )
        jVs.append(jv_data)
        uvvis.append(uvvis_data)
        eqes.append(eqe_data)
        optical.append(optical_data)
        parameters.append(params)

    v, j, ff, voc, jsc, pces = save_jV_data(jVs, data_dir)
    absorption, reflection, wavelenghts_abs = save_uvvis_data(uvvis, data_dir)
    n, k, wavelenghts_opt = save_optical_data(optical, data_dir)
    eqe, wavelenghts_eqe = save_eqe_data(eqes, data_dir)
    parameters = save_parameters(parameters, data_dir)

    plot_device_characteristics_distributions(
        pces=pces,
        jscs=jsc,
        vocs=voc,
        ffs=ff,
        save_path=plot_dir / "device_characteristics_distributions.png",
        dpi=DPI,
    )

    num_to_plot = (
        len(jVs) if MAX_CURVES_TO_PLOT is None else min(MAX_CURVES_TO_PLOT, len(jVs))
    )

    plot_jVs(
        jVs=jVs[:num_to_plot],
        save_path=plot_dir / "jV.png",
        dpi=DPI,
    )

    plot_EQEs(
        eqes=eqes[:num_to_plot],
        save_path=plot_dir / "EQE.png",
        dpi=DPI,
    )

    plot_uvvis(
        uvvis=uvvis[:num_to_plot],
        save_dir=plot_dir,
        dpi=DPI,
    )

    find_devices_with_best_pce(
        pces=pces,
        run_dirs=run_dirs,
        save_dir=SAVE_DIR,
        top_k=TOP_K,
        jvs=jVs,
    )


def main():
    simulation_dir = Path.cwd() / "simulation_runs"
    post_process_simulations(simulation_dir)


if __name__ == "__main__":
    main()
