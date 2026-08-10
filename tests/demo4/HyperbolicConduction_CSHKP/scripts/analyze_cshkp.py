#!/usr/bin/env python3
"""Validate and plot the fixed-resistivity CSHKP transport comparison."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yt


GAMMA_MINUS_ONE = 2.0 / 3.0
UNIT_TIME_S = 77.82914611581296
UNIT_LENGTH_MM = 10.0
SHEET_Y = 6.0

CASES = {
    "parallel": {
        "stem": "flare00_parallel_",
        "label": "matched parallel",
        "color": "#3F3F3F",
        "linestyle": "-",
    },
    "electron": {
        "stem": "flare01_electron_magnetization_",
        "label": "electron magnetization",
        "color": "#4C78A8",
        "linestyle": "--",
    },
    "effective_1e-3": {
        "stem": "flare03_effective_1e-3_",
        "label": r"effective $f_\perp=10^{-3}$",
        "color": "#E6A32F",
        "linestyle": "-.",
    },
    "effective_1e-2": {
        "stem": "flare04_effective_1e-2_",
        "label": r"effective $f_\perp=10^{-2}$",
        "color": "#D55E5E",
        "linestyle": "-",
    },
}

LOG_COLUMNS = [
    "time",
    "E_excess",
    "E_core",
    "E_halo",
    "halo_fraction",
    "width_rms",
    "Tmax",
    "Mdot_up",
    "Erec_max",
    "vdrift_max",
]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def read_logs(data_dir: Path) -> dict[str, np.ndarray]:
    logs: dict[str, np.ndarray] = {}
    reference_time = None
    for key, meta in CASES.items():
        path = data_dir / f"{meta['stem']}.log"
        values = np.loadtxt(path)
        if values.ndim != 2 or values.shape[1] != len(LOG_COLUMNS):
            raise ValueError(f"Unexpected log shape in {path}: {values.shape}")
        if reference_time is None:
            reference_time = values[:, 0]
        elif not np.allclose(values[:, 0], reference_time, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"Time samples do not match the parallel baseline: {path}")
        logs[key] = values
    return logs


def relative_to_parallel(logs: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    baseline = logs["parallel"]
    denominator = np.where(np.abs(baseline[:, 1:]) > 1.0e-14, np.abs(baseline[:, 1:]), np.nan)
    return {
        key: (values[:, 1:] - baseline[:, 1:]) / denominator
        for key, values in logs.items()
    }


def load_uniform_snapshot(path: Path) -> dict[str, np.ndarray | float]:
    logging.getLogger("yt").setLevel(logging.ERROR)
    ds = yt.load(path)
    level = int(ds.max_level)
    dims = np.asarray(ds.domain_dimensions, dtype=int) * (2**level)
    grid = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    def field(name: str) -> np.ndarray:
        return np.asarray(grid[("amrvac", name)])[:, :, 0]

    rho = field("rho")
    internal_energy = field("e")
    temperature = GAMMA_MINUS_ONE * internal_energy / rho
    nx, ny = temperature.shape
    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    x = np.linspace(left[0], right[0], nx, endpoint=False) + (right[0] - left[0]) / (2 * nx)
    y = np.linspace(left[1], right[1], ny, endpoint=False) + (right[1] - left[1]) / (2 * ny)
    b1 = field("b1")
    b2 = field("b2")
    by0 = -40.0 * np.tanh((20.0 / 3.0) * x)[:, None]
    return {
        "time": float(ds.current_time),
        "x": x,
        "y": y,
        "temperature": temperature,
        "bx": b1,
        "by": b2 + by0,
    }


def write_source_data(
    out_dir: Path,
    logs: dict[str, np.ndarray],
    rel: dict[str, np.ndarray],
    snapshots: dict[str, dict[str, np.ndarray | float]],
) -> None:
    source_dir = out_dir / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    time = logs["parallel"][:, 0]
    with (source_dir / "time_series_relative.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_code",
                "time_s",
                "case",
                "halo_fraction",
                "width_rms_code",
                "Tmax_MK",
                "Erec_max_code",
                "delta_halo_fraction_percent",
                "delta_width_rms_percent",
                "delta_Tmax_percent",
                "delta_Erec_max_percent",
            ]
        )
        for key in CASES:
            for i, t in enumerate(time):
                writer.writerow(
                    [
                        f"{t:.10e}",
                        f"{t * UNIT_TIME_S:.10e}",
                        key,
                        f"{logs[key][i, 4]:.10e}",
                        f"{logs[key][i, 5]:.10e}",
                        f"{logs[key][i, 6]:.10e}",
                        f"{logs[key][i, 8]:.10e}",
                        f"{100.0 * rel[key][i, 3]:.10e}",
                        f"{100.0 * rel[key][i, 4]:.10e}",
                        f"{100.0 * rel[key][i, 5]:.10e}",
                        f"{100.0 * rel[key][i, 7]:.10e}",
                    ]
                )

    x = np.asarray(snapshots["parallel"]["x"])
    y = np.asarray(snapshots["parallel"]["y"])
    iy = int(np.argmin(np.abs(y - SHEET_Y)))
    with (source_dir / "cross_sheet_temperature_t0p5.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x_code", "x_Mm", "y_code"] + [f"T_MK_{key}" for key in CASES])
        for ix, x_value in enumerate(x):
            writer.writerow(
                [
                    f"{x_value:.10e}",
                    f"{x_value * UNIT_LENGTH_MM:.10e}",
                    f"{y[iy]:.10e}",
                ]
                + [f"{np.asarray(snapshots[key]['temperature'])[ix, iy]:.10e}" for key in CASES]
            )

    crop_x = np.abs(x) <= 1.0
    crop_y = (y >= 5.0) & (y <= 7.0)
    np.savez_compressed(
        source_dir / "temperature_maps_t0p5.npz",
        x_code=x[crop_x],
        y_code=y[crop_y],
        **{
            f"T_MK_{key}": np.asarray(snapshots[key]["temperature"])[np.ix_(crop_x, crop_y)]
            for key in CASES
        },
    )


def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def make_figure(
    out_dir: Path,
    logs: dict[str, np.ndarray],
    rel: dict[str, np.ndarray],
    snapshots: dict[str, dict[str, np.ndarray | float]],
) -> Path:
    configure_matplotlib()
    fig = plt.figure(figsize=(7.20, 4.35), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1.03, 0.86], width_ratios=[1.0, 1.0, 1.18])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1], sharex=ax_a, sharey=ax_a)
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[1, :2])
    ax_e = fig.add_subplot(grid[1, 2])

    base = snapshots["parallel"]
    eff = snapshots["effective_1e-2"]
    x = np.asarray(base["x"])
    y = np.asarray(base["y"])
    temperature = np.asarray(base["temperature"])
    delta_temperature = np.asarray(eff["temperature"]) - temperature
    bx = np.asarray(base["bx"])
    by = np.asarray(base["by"])
    crop_x = np.abs(x) <= 1.0
    crop_y = (y >= 5.0) & (y <= 7.0)
    x_crop = x[crop_x]
    y_crop = y[crop_y]
    temp_crop = temperature[np.ix_(crop_x, crop_y)]
    delta_crop = delta_temperature[np.ix_(crop_x, crop_y)]

    image_a = ax_a.pcolormesh(
        x_crop * UNIT_LENGTH_MM,
        y_crop * UNIT_LENGTH_MM,
        temp_crop.T,
        cmap="magma",
        shading="auto",
        vmin=1.0,
        vmax=2.0,
        rasterized=True,
    )
    stride = 4
    ax_a.streamplot(
        x_crop[::stride] * UNIT_LENGTH_MM,
        y_crop[::stride] * UNIT_LENGTH_MM,
        bx[np.ix_(crop_x, crop_y)].T[::stride, ::stride],
        by[np.ix_(crop_x, crop_y)].T[::stride, ::stride],
        color="white",
        linewidth=0.25,
        density=0.75,
        arrowsize=0.35,
    )
    colorbar_a = fig.colorbar(image_a, ax=ax_a, pad=0.02, fraction=0.048)
    colorbar_a.set_label(r"$T$ (MK)")
    ax_a.set_title(r"matched parallel, $t=0.5$")
    ax_a.set_xlabel(r"$x$ (Mm)")
    ax_a.set_ylabel(r"$y$ (Mm)")
    panel_label(ax_a, "a")

    delta_limit = max(0.02, float(np.nanmax(np.abs(delta_crop))))
    image_b = ax_b.pcolormesh(
        x_crop * UNIT_LENGTH_MM,
        y_crop * UNIT_LENGTH_MM,
        delta_crop.T,
        cmap="RdBu_r",
        shading="auto",
        vmin=-delta_limit,
        vmax=delta_limit,
        rasterized=True,
    )
    ax_b.contour(
        x_crop * UNIT_LENGTH_MM,
        y_crop * UNIT_LENGTH_MM,
        temp_crop.T,
        levels=[1.1, 1.5],
        colors=["#303030", "#303030"],
        linewidths=[0.35, 0.55],
    )
    colorbar_b = fig.colorbar(image_b, ax=ax_b, pad=0.02, fraction=0.048)
    colorbar_b.set_label(r"$\Delta T$ (MK)")
    ax_b.set_title(r"$f_\perp=10^{-2}$ minus parallel")
    ax_b.set_xlabel(r"$x$ (Mm)")
    ax_b.tick_params(labelleft=False)
    panel_label(ax_b, "b")

    iy = int(np.argmin(np.abs(y - SHEET_Y)))
    profile_mask = np.abs(x) <= 0.55
    baseline_profile = np.asarray(snapshots["parallel"]["temperature"])[:, iy]
    for key in ["electron", "effective_1e-3", "effective_1e-2"]:
        meta = CASES[key]
        values = np.asarray(snapshots[key]["temperature"])
        ax_c.plot(
            x[profile_mask] * UNIT_LENGTH_MM,
            1.0e3 * (values[profile_mask, iy] - baseline_profile[profile_mask]),
            color=meta["color"],
            linestyle=meta["linestyle"],
            linewidth=1.15 if key == "effective_1e-2" else 0.9,
            label=meta["label"],
        )
    ax_c.axhline(0.0, color="#A0A0A0", linewidth=0.55)
    ax_c.set_xlabel(r"$x$ at $y\simeq60$ Mm (Mm)")
    ax_c.set_ylabel(r"$\Delta T$ (kK)")
    ax_c.set_title("cross-sheet thermal response")
    ax_c.legend(loc="lower left", handlelength=2.2)
    panel_label(ax_c, "c")

    time_s = logs["parallel"][:, 0] * UNIT_TIME_S
    valid = logs["parallel"][:, 0] <= 0.5000001
    for key in ["electron", "effective_1e-3", "effective_1e-2"]:
        meta = CASES[key]
        ax_d.plot(
            time_s[valid],
            100.0 * rel[key][valid, 3],
            color=meta["color"],
            linestyle=meta["linestyle"],
            linewidth=1.2 if key == "effective_1e-2" else 0.95,
            label=meta["label"],
        )
    ax_d.axhline(0.0, color="#A0A0A0", linewidth=0.55)
    ax_d.set_xlabel("time (s)")
    ax_d.set_ylabel(r"$\Delta$ halo fraction (%)")
    ax_d.set_title("cross-field thermal redistribution")
    ax_d.legend(loc="upper right", ncol=1)
    panel_label(ax_d, "d")

    for key in ["electron", "effective_1e-3", "effective_1e-2"]:
        meta = CASES[key]
        ax_e.plot(
            time_s[valid],
            100.0 * rel[key][valid, 7],
            color=meta["color"],
            linestyle=meta["linestyle"],
            linewidth=1.2 if key == "effective_1e-2" else 0.95,
        )
    ax_e.axhline(0.0, color="#A0A0A0", linewidth=0.55)
    ax_e.set_xlabel("time (s)")
    ax_e.set_ylabel(r"$\Delta E_{\rm rec,max}$ (%)")
    ax_e.set_title("matched reconnection control")
    panel_label(ax_e, "e")

    for ax in [ax_a, ax_b]:
        ax.set_aspect("equal", adjustable="box")
    for ax in [ax_c, ax_d, ax_e]:
        ax.grid(color="#D8D8D8", linewidth=0.35, alpha=0.7)

    figure_base = out_dir / "cshkp_transport_candidate"
    fig.savefig(figure_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(figure_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(figure_base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    preview = figure_base.with_suffix(".png")
    fig.savefig(preview, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return preview


def write_validation_summary(
    out_dir: Path,
    logs: dict[str, np.ndarray],
    rel: dict[str, np.ndarray],
) -> None:
    baseline = logs["parallel"]
    valid = baseline[:, 0] <= 0.5000001
    rows = []
    for key in CASES:
        i02 = int(np.argmin(np.abs(baseline[:, 0] - 0.2)))
        i05 = int(np.argmin(np.abs(baseline[:, 0] - 0.5)))
        rows.append(
            {
                "case": key,
                "halo_change_t0p2_percent": 100.0 * rel[key][i02, 3],
                "halo_change_t0p5_percent": 100.0 * rel[key][i05, 3],
                "max_abs_halo_change_pre0p5_percent": 100.0 * np.nanmax(np.abs(rel[key][valid, 3])),
                "max_abs_width_change_pre0p5_percent": 100.0 * np.nanmax(np.abs(rel[key][valid, 4])),
                "max_abs_Tmax_change_pre0p5_percent": 100.0 * np.nanmax(np.abs(rel[key][valid, 5])),
                "max_abs_Erec_change_pre0p5_percent": 100.0 * np.nanmax(np.abs(rel[key][valid, 7])),
            }
        )
    path = out_dir / "validation_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("production/datamr"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--snapshot", type=int, default=5, help="Snapshot index used for spatial panels")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    logs = read_logs(args.data_dir)
    rel = relative_to_parallel(logs)
    snapshots = {
        key: load_uniform_snapshot(args.data_dir / f"{meta['stem']}{args.snapshot:04d}.dat")
        for key, meta in CASES.items()
    }
    write_source_data(args.output_dir, logs, rel, snapshots)
    write_validation_summary(args.output_dir, logs, rel)
    preview = make_figure(args.output_dir, logs, rel, snapshots)
    print(preview)


if __name__ == "__main__":
    main()
