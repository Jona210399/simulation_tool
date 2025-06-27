from pathlib import Path

from simulation_tool.postprocessing.core import post_process_session
from simulation_tool.postprocessing.plotting import (
    plot_device_characteristics_distributions,
    plot_EQEs,
    plot_jVs,
    plot_optical_data,
    plot_uvvis,
)
from simulation_tool.utils import get_sorted_run_dirs, save_figure


def main():
    simualtion_dir = Path.cwd() / "simulation_runs" / "697856" / "post_processed"

    best_pces_dir = simualtion_dir / "best_pces"

    jvs = []
    uvvis = []
    eqes = []
    opticals = []

    for run_dir in get_sorted_run_dirs(best_pces_dir):
        plot_dir = run_dir / "plots"
        plot_dir.mkdir(parents=True, exist_ok=True)
        jv, uv, eqe, optical, _ = post_process_session(run_dir)

        save_figure(jv, save_path=plot_dir / "jV.png", dpi=300)
        save_figure(uv, save_path=plot_dir / "uvvis.png", dpi=300)
        save_figure(eqe, save_path=plot_dir / "eqe.png", dpi=300)
        save_figure(optical, save_path=plot_dir / "nk.png", dpi=300)

        jvs.append(jv)
        uvvis.append(uv)
        eqes.append(eqe)
        opticals.append(optical)

    plot_dir = best_pces_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plot_jVs(
        jVs=jvs,
        save_path=plot_dir / "jVs.png",
        dpi=300,
        label=True,
        ylim=(-250, 200),
    )

    plot_uvvis(
        uvvis=uvvis,
        save_dir=plot_dir,
        dpi=300,
    )

    plot_EQEs(
        eqes=eqes,
        save_path=plot_dir / "eqes.png",
        dpi=300,
    )

    plot_optical_data(
        opticals=opticals,
        save_path=plot_dir / "nks.png",
        dpi=300,
    )


def plot_post_processed_statistics():
    post_processed_dir = Path.cwd() / "simulation_runs" / "post_processed"

    from simulation_tool.postprocessing.core import read_post_processed_data

    jv_df, *rest = read_post_processed_data(post_processed_dir)

    plot_device_characteristics_distributions(
        pces=jv_df["PCE"].to_numpy(),
        jscs=jv_df["Jsc"].to_numpy(),
        vocs=jv_df["Voc"].to_numpy(),
        ffs=jv_df["FF"].to_numpy(),
        save_path=post_processed_dir / "device_characteristics_distributions.png",
        dpi=300,
    )


if __name__ == "__main__":
    main()
