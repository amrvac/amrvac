#!/usr/bin/env python3
"""Generate the fixed, finite Fourier benchmark matrix and parameter overlays."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARAMS = ROOT / "params" / "production"
MATRIX = ROOT / "matrix.csv"

ANGLES = {
    "parallel": (2, 1, 0.0),
    "mixed": (1, 2, math.degrees(math.acos(4.0 / 5.0))),
    "perpendicular": (-1, 2, 90.0),
}
RESOLUTIONS = (64, 128, 256)
TAUS = (2.5e-4, 1.25e-4, 6.25e-5)
KAPPA_PARALLEL = 150.0


def tau_tag(tau: float) -> str:
    return f"{tau:.3e}".replace(".", "p").replace("-", "m").replace("+", "")


def write_overlay(row: dict[str, object]) -> None:
    path = ROOT / str(row["overlay"])
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "&filelist",
        f"  base_filename='runs/{row['case']}/output'",
        "/",
        "&meshlist",
        f"  domain_nx1={row['N']}",
        f"  domain_nx2={row['N']}",
        "/",
    ]
    if row["method"] == "hyperbolic":
        lines.extend(
            [
                "&paramlist",
                f"  courantpar={float(row['courantpar']):.16e}",
                f"  dtpar={float(row['dtpar']):.16e}",
                "/",
                "&mhd_list",
                f"  mhd_hyperbolic_tc_kappa_perp_factor={float(row['ratio']):.16e}",
                "/",
            ]
        )
    else:
        lines.extend(
            [
                "&tc_list",
                f"  tc_k_para={KAPPA_PARALLEL:.16e}",
                f"  tc_k_perp={KAPPA_PARALLEL * float(row['ratio']):.16e}",
                "/",
            ]
        )
    lines.extend(
        [
            "&usr_list",
            f"  kx_mode={row['kx']}",
            f"  ky_mode={row['ky']}",
            f"  kappa_perp_ratio_usr={float(row['ratio']):.16e}",
            f"  tau_target={float(row['tau_target'] or 0.0):.16e}",
            "/",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def make_row(
    suite: str,
    method: str,
    angle: str,
    n: int,
    ratio: float,
    tau: float | None = None,
) -> dict[str, object]:
    kx, ky, theta = ANGLES[angle]
    if method == "parabolic":
        suffix = "p"
        dtpar = ""
        courant = ""
        parfile = "parabolic.par"
    else:
        assert tau is not None
        suffix = f"h_tau{tau_tag(tau)}"
        dtpar = tau / 8.0
        courant = math.sqrt(tau / 60.0)
        parfile = "hyperbolic.par"
    case = f"{suite}_{angle}_n{n:03d}_{suffix}"
    return {
        "case": case,
        "suite": suite,
        "method": method,
        "angle": angle,
        "kx": kx,
        "ky": ky,
        "theta_deg": theta,
        "ratio": ratio,
        "N": n,
        "tau_target": "" if tau is None else tau,
        "dtpar": dtpar,
        "courantpar": courant,
        "parfile": parfile,
        "overlay": f"params/production/{case}.par",
        "run_dir": f"runs/{case}",
    }


def build_matrix() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in RESOLUTIONS:
        for angle in ANGLES:
            rows.append(make_row("aniso", "parabolic", angle, n, 0.1))
            for tau in TAUS:
                rows.append(make_row("aniso", "hyperbolic", angle, n, 0.1, tau))
    for angle in ANGLES:
        rows.append(make_row("isotropic", "parabolic", angle, 128, 1.0))
        rows.append(make_row("isotropic", "hyperbolic", angle, 128, 1.0, TAUS[-1]))
    return rows


def main() -> None:
    PARAMS.mkdir(parents=True, exist_ok=True)
    rows = build_matrix()
    for row in rows:
        write_overlay(row)
    with MATRIX.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} cases to {MATRIX}")


if __name__ == "__main__":
    main()
