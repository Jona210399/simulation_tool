import os
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import polars as pl
from numpy.typing import NDArray
from tqdm import tqdm

from simulation_tool.constants import (
    EQE_SIM_OUTPUT_FILE_NAME,
    JV_SIM_OUTPUT_FILE_NAME,
    NK_FILE_NAME,
    PARAMS_OUTPUT_FILE_NAME,
    RUN_DIR_PREFIX,
    UVVIS_FILE_NAME,
)
from simulation_tool.data import OpticalData, UVVisData
from simulation_tool.EQE import EQEData
from simulation_tool.jV import JVData
from simulation_tool.postprocessing.core import post_process_session
from simulation_tool.postprocessing.plotting import (
    plot_device_characteristics_distributions,
    plot_EQEs,
    plot_jVs,
    plot_uvvis,
)
from simulation_tool.postprocessing.saving import (
    save_eqe_data,
    save_jV_data,
    save_optical_data,
    save_parameters,
    save_uvvis_data,
)
from simulation_tool.utils import get_sorted_run_dirs


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


def post_process_simulations(
    simulation_dir: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    MAX_CURVES_TO_PLOT: int | None = 50
    TOP_K: int = 5
    DPI: int = 300
    SAVE_DIR = simulation_dir / "post_processed"
    if SAVE_DIR.exists():
        shutil.rmtree(SAVE_DIR)
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

    for session_path in tqdm(run_dirs):
        jv_data, uvvis_data, eqe_data, optical_data, params = post_process_session(
            session_path
        )
        jVs.append(jv_data)
        uvvis.append(uvvis_data)
        eqes.append(eqe_data)
        optical.append(optical_data)
        parameters.append(params)

    jV_df, (v, j, ff, voc, jsc, pces) = save_jV_data(jVs, data_dir)
    uvvis_df, (absorption, reflection, transmission, wavelenghts_abs) = save_uvvis_data(
        uvvis, data_dir
    )
    eqes_df, (eqe, wavelenghts_eqe) = save_eqe_data(eqes, data_dir)
    optical_df, (n, k, wavelenghts_opt) = save_optical_data(optical, data_dir)
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

    return jV_df, uvvis_df, eqes_df, optical_df, parameters


def read_post_processed_data(
    dir: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    jVs = pl.read_parquet((dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    uvvis = pl.read_parquet((dir / UVVIS_FILE_NAME).with_suffix(".parquet"))
    eqes = pl.read_parquet((dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    optical = pl.read_parquet((dir / NK_FILE_NAME).with_suffix(".parquet"))
    parameters = pl.read_parquet(
        (dir / PARAMS_OUTPUT_FILE_NAME).with_suffix(".parquet")
    )

    return jVs, uvvis, eqes, optical, parameters


def main():
    simulation_dir = Path.cwd() / "simulation_runs"
    max_workers = int(os.getenv("SLURM_CPUS_PER_TASK", 1))

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        names = []
        for dir_ in list(simulation_dir.iterdir()):
            if not dir_.is_dir():
                continue

            if (already_processed := dir_ / "post_processed").exists():
                print(f"Skipping {dir_} as it has already been post-processed.")
                futures.append(
                    executor.submit(read_post_processed_data, already_processed)
                )

            print(f"Found run directory: {dir_}")
            names.append(dir_.name)
            futures.append(executor.submit(post_process_simulations, dir_))

        jv_dfs = []
        uvvis_dfs = []
        eqe_dfs = []
        optical_dfs = []
        parameters_dfs = []

        with tqdm(total=len(futures), desc="Processing directories") as pbar:
            for future in as_completed(futures):
                pbar.update(1)
                jVs, uvvis, eqes, optical, parameters = future.result()
                jv_dfs.append(jVs)
                uvvis_dfs.append(uvvis)
                eqe_dfs.append(eqes)
                optical_dfs.append(optical)
                parameters_dfs.append(parameters)

    jv_df: pl.DataFrame = pl.concat(jv_dfs)
    uvvis_df: pl.DataFrame = pl.concat(uvvis_dfs)
    eqe_df: pl.DataFrame = pl.concat(eqe_dfs)
    optical_df: pl.DataFrame = pl.concat(optical_dfs)
    parameters_df: pl.DataFrame = pl.concat(parameters_dfs)

    output_dir = simulation_dir / "post_processed"
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "jobs.txt", "w") as f:
        f.write("\n".join(names))

    jv_df.write_parquet((output_dir / JV_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet"))
    uvvis_df.write_parquet((output_dir / UVVIS_FILE_NAME).with_suffix(".parquet"))
    eqe_df.write_parquet(
        (output_dir / EQE_SIM_OUTPUT_FILE_NAME).with_suffix(".parquet")
    )
    optical_df.write_parquet((output_dir / NK_FILE_NAME).with_suffix(".parquet"))
    parameters_df.write_parquet(
        (output_dir / PARAMS_OUTPUT_FILE_NAME).with_suffix(".parquet")
    )
    print(f"Post-processing completed. Data saved to {output_dir}")


if __name__ == "__main__":
    main()
