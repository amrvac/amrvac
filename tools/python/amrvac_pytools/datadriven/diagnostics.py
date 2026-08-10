"""Magnetofrictional relaxation diagnostics helpers.

This module keeps the notebook-facing diagnostics path lightweight.  Pandas is
not required; CSV files are read with the standard library and converted to
NumPy arrays for plotting and snapshot selection.
"""

from __future__ import print_function

import csv
import math
import re
from pathlib import Path


_ALIASES = {
    "iteration": ("iteration", "iter", "it", "itmf", "step", "nstep", "istep"),
    "time": ("time", "t", "physical_time", "simulation_time"),
    "divb": (
        "divb",
        "div_b",
        "divergence",
        "divergence_b",
        "max_divb",
        "mean_divb",
        "f_i",
        "fi",
    ),
    "cw_sin_theta": (
        "cw_sin_theta",
        "cw sintheta",
        "cw sin theta",
        "cwsintheta",
        "cw_sin",
        "sigmaj",
        "sigma_j",
    ),
    "current": ("current", "mean_current", "average_current", "j"),
    "lorentz_force": (
        "lorentz_force",
        "lorenz_force",
        "lorentz",
        "lorenz",
        "force",
        "jxb",
        "j_cross_b",
    ),
    # Backward-compatible canonical name used by the first V1 implementation.
    "jxb": (
        "jxb",
        "j_cross_b",
        "force",
        "lorentz",
        "lorenz",
        "lorentz_force",
        "lorenz_force",
    ),
    "magnetic_energy": (
        "magnetic_energy",
        "energy",
        "eb",
        "e_mag",
        "emag",
        "mag_energy",
    ),
    "dt": ("dt", "delta_t", "timestep", "time_step"),
}


_SNAPSHOT_RE = re.compile(r"^(?P<stem>.*?)(?P<number>\d+)\.dat$")


def read_relaxation_diagnostics(csv_path):
    """Read a magnetofrictional relaxation diagnostics CSV.

    Returns a notebook-friendly dictionary with original columns preserved in
    ``data`` and recognized aliases listed in ``canonical_columns``.
    """

    csv_path = Path(csv_path).expanduser().resolve()
    if not csv_path.exists():
        raise FileNotFoundError("cannot find relaxation diagnostics CSV: {}".format(csv_path))

    try:
        import numpy as np
    except ImportError as error:  # pragma: no cover
        raise ImportError("numpy is required to read relaxation diagnostics") from error

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = _sniff_csv_dialect(sample)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError("diagnostics CSV has no header: {}".format(csv_path))
        columns = [_clean_column_name(name) for name in reader.fieldnames]
        values = {column: [] for column in columns}
        for row in reader:
            for raw_name, column in zip(reader.fieldnames, columns):
                values[column].append(_to_float(row.get(raw_name)))

    data = {name: np.asarray(vals, dtype=float) for name, vals in values.items()}
    canonical = _canonical_columns(columns)
    row_count = len(next(iter(data.values()))) if data else 0
    return {
        "csv_path": str(csv_path),
        "columns": columns,
        "canonical_columns": canonical,
        "data": data,
        "row_count": row_count,
    }


def plot_relaxation_diagnostics(
    diagnostics,
    output_path=None,
    columns=None,
    x_column=None,
    case_dir=None,
    output_dir=None,
    base_filename="output/data_driven_mfr",
    mf_ditsave=None,
    restart_markers=None,
    figsize=(10, 5.5),
    dpi=150,
):
    """Plot key relaxation diagnostics and mark available restart snapshots."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as error:  # pragma: no cover
        raise ImportError("matplotlib is required to plot relaxation diagnostics") from error

    diag = _as_diagnostics(diagnostics)
    data = diag["data"]
    canonical = diag["canonical_columns"]
    warnings = []

    x_name = _resolve_column(diag, x_column) if x_column else None
    if x_name is None:
        x_name = canonical.get("iteration") or canonical.get("time")
    if x_name is None:
        try:
            import numpy as np
        except ImportError as error:  # pragma: no cover
            raise ImportError("numpy is required to plot relaxation diagnostics") from error

        x_values = np.arange(diag["row_count"], dtype=float)
        x_label = "row"
    else:
        x_values = data[x_name]
        x_label = _diagnostic_label(diag, x_name)

    y_columns = _default_plot_columns(diag) if columns is None else [
        _resolve_column(diag, column) for column in columns
    ]
    y_columns = [column for column in y_columns if column is not None and column != x_name]
    if not y_columns:
        raise ValueError("no diagnostics columns available to plot")

    if restart_markers is None:
        marker_summary = find_relaxation_restart_snapshots(
            case_dir=case_dir,
            output_dir=output_dir,
            base_filename=base_filename,
            mf_ditsave=mf_ditsave,
            diagnostics=diag,
        )
        restart_markers = marker_summary["markers"]
        warnings.extend(marker_summary["warnings"])
    else:
        restart_markers = list(restart_markers)

    marker_x_values = _marker_x_values(restart_markers, x_name)
    selected_iteration = None
    for marker in restart_markers:
        if marker.get("selected"):
            selected_iteration = marker.get("iteration")
            break

    if len(y_columns) <= 2:
        fig, primary_axis = plt.subplots(figsize=figsize)
        axes = [primary_axis]
        if len(y_columns) == 2:
            axes.append(primary_axis.twinx())

        colors = ("tab:blue", "tab:orange")
        legend_handles = []
        legend_labels = []
        for axis, column, color in zip(axes, y_columns, colors):
            label = _diagnostic_label(diag, column)
            line, = axis.plot(x_values, data[column], color=color, linewidth=1.3)
            axis.set_ylabel(label, color=color)
            axis.tick_params(axis="y", labelcolor=color)
            legend_handles.append(line)
            legend_labels.append(label)

        primary_axis.grid(True, alpha=0.25)
        saved_handle = None
        for marker_x in marker_x_values:
            marker_line = primary_axis.axvline(
                marker_x,
                color="tab:red",
                linewidth=0.8,
                alpha=0.35,
            )
            if saved_handle is None:
                saved_handle = marker_line
        if saved_handle is not None:
            legend_handles.append(saved_handle)
            legend_labels.append("saved restart")
        if selected_iteration is not None:
            selected_handle = primary_axis.axvline(
                selected_iteration,
                color="tab:green",
                linewidth=1.4,
                alpha=0.8,
            )
            legend_handles.append(selected_handle)
            legend_labels.append("selected restart")

        primary_axis.set_xlabel(x_label)
        primary_axis.set_title("Magnetofrictional relaxation diagnostics")
        primary_axis.legend(legend_handles, legend_labels, loc="best")
    else:
        fig, subplot_axes = plt.subplots(
            len(y_columns),
            1,
            sharex=True,
            figsize=(figsize[0], max(figsize[1], 2.6 * len(y_columns))),
        )
        axes = list(subplot_axes)
        for axis, column in zip(axes, y_columns):
            axis.plot(x_values, data[column], linewidth=1.2)
            axis.set_ylabel(_diagnostic_label(diag, column))
            axis.grid(True, alpha=0.25)
            for marker_x in marker_x_values:
                axis.axvline(marker_x, color="tab:red", linewidth=0.8, alpha=0.35)
            if selected_iteration is not None:
                axis.axvline(selected_iteration, color="tab:green", linewidth=1.4, alpha=0.8)
        axes[-1].set_xlabel(x_label)
        axes[0].set_title("Magnetofrictional relaxation diagnostics")
    fig.tight_layout()

    saved_path = None
    if output_path is not None:
        output_path = Path(output_path).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(output_path), dpi=dpi)
        saved_path = str(output_path)

    return {
        "figure": fig,
        "axes": axes,
        "output_path": saved_path,
        "x_column": x_name,
        "plotted_columns": y_columns,
        "restart_markers": restart_markers,
        "warnings": warnings,
    }


def select_relaxation_restart(
    diagnostics=None,
    case_dir=None,
    output_dir=None,
    base_filename="output/data_driven_mfr",
    mf_ditsave=None,
    metric=None,
    target_iteration=None,
    snapshot_number=None,
):
    """Select a restart snapshot from saved MFR checkpoints.

    If ``target_iteration`` or ``snapshot_number`` is provided, the nearest or
    exact saved snapshot is selected.  Otherwise Cw sin(theta) and Lorentz force
    are combined with equal rank weight at the saved checkpoints.  Divergence is
    used as a stability guard, while current is reported but is deliberately not
    minimized because a force-free field may carry electric current.
    """

    diag = _as_diagnostics(diagnostics) if diagnostics is not None else None
    snapshot_summary = find_relaxation_restart_snapshots(
        case_dir=case_dir,
        output_dir=output_dir,
        base_filename=base_filename,
        mf_ditsave=mf_ditsave,
        diagnostics=diag,
    )
    markers = snapshot_summary["markers"]
    warnings = list(snapshot_summary["warnings"])
    if not markers:
        return {
            "restart_file": None,
            "selected_marker": None,
            "metric": None,
            "value": None,
            "markers": markers,
            "warnings": warnings + ["no restart snapshot files were found"],
        }

    selected = None
    selected_metric = None
    selected_value = None
    selection_reason = None
    if snapshot_number is not None:
        wanted = int(snapshot_number)
        for marker in markers:
            if marker["snapshot_number"] == wanted:
                selected = marker
                break
        if selected is None:
            raise ValueError("snapshot_number {} was not found".format(wanted))
        selection_reason = "restart snapshot was selected explicitly"
    elif target_iteration is not None:
        target = float(target_iteration)
        with_iterations = [marker for marker in markers if marker.get("iteration") is not None]
        if with_iterations:
            selected = min(with_iterations, key=lambda marker: abs(marker["iteration"] - target))
            selection_reason = "nearest saved restart to the requested iteration"
        else:
            warnings.append("cannot select nearest iteration because snapshot iterations were not inferred")
    else:
        if metric is None:
            trend_selection = _select_from_relaxation_trends(diag, markers)
            selected = trend_selection["selected"]
            selected_metric = trend_selection["metric"]
            selected_value = trend_selection["value"]
            selection_reason = trend_selection["reason"]
            warnings.extend(trend_selection["warnings"])
        else:
            selected_metric = _resolve_selection_metric(diag, metric)
            candidates = []
            for marker in markers:
                value = _marker_metric_value(marker, selected_metric)
                if value is not None and math.isfinite(value):
                    candidates.append((value, marker))
            if candidates:
                selected_value, selected = min(candidates, key=lambda item: item[0])
                selection_reason = "minimum saved-snapshot value of {}".format(selected_metric)
            else:
                warnings.append("no finite saved-snapshot values found for {}".format(selected_metric))

    if selected is None:
        selected = _latest_marker(markers)
        warnings.append("selected the latest available restart snapshot as a fallback")
        selection_reason = "latest available restart fallback"
    if selected_metric is not None and selected_value is None:
        selected_value = _marker_metric_value(selected, selected_metric)

    return {
        "restart_file": selected["path"],
        "selected_marker": selected,
        "metric": selected_metric,
        "value": selected_value,
        "reason": selection_reason,
        "diagnostic_values": _selected_diagnostic_values(diag, selected),
        "markers": markers,
        "warnings": warnings,
    }


def find_relaxation_restart_snapshots(
    case_dir=None,
    output_dir=None,
    base_filename="output/data_driven_mfr",
    mf_ditsave=None,
    diagnostics=None,
):
    """Find saved MFR restart snapshots and infer their CSV locations."""

    diag = _as_diagnostics(diagnostics) if diagnostics is not None else None
    search_dir, stem = _snapshot_search_location(case_dir, output_dir, base_filename)
    warnings = []
    if search_dir is None:
        warnings.append("case_dir or output_dir is required to find restart snapshots")
        return {"markers": [], "warnings": warnings}
    if not search_dir.exists():
        warnings.append("restart output directory does not exist: {}".format(search_dir))
        return {"markers": [], "warnings": warnings}

    files = []
    for path in search_dir.glob("{}*.dat".format(stem)):
        number = _snapshot_number(path.name, stem)
        if number is not None:
            files.append((number, path.resolve()))
    files.sort(key=lambda item: item[0])
    if not files:
        warnings.append("no restart snapshots matching {}*.dat in {}".format(stem, search_dir))
        return {"markers": [], "warnings": warnings}

    iteration_map, map_warnings = _infer_snapshot_iterations(files, mf_ditsave, diag)
    warnings.extend(map_warnings)
    markers = []
    for number, path in files:
        marker = {
            "path": str(path),
            "name": path.name,
            "snapshot_number": number,
            "iteration": iteration_map.get(number),
            "row_index": None,
            "values": {},
        }
        if diag is not None and marker["iteration"] is not None:
            row_index = _nearest_iteration_row(diag, marker["iteration"])
            if row_index is not None:
                marker["row_index"] = row_index
                for column, values in diag["data"].items():
                    marker["values"][column] = float(values[row_index])
        markers.append(marker)
    return {"markers": markers, "warnings": warnings}


def _sniff_csv_dialect(sample):
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t ")
    except csv.Error:
        return csv.excel


def _clean_column_name(name):
    return str(name).strip()


def _normalize_column_name(name):
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def _canonical_columns(columns):
    normalized = {_normalize_column_name(column): column for column in columns}
    canonical = {}
    for key, aliases in _ALIASES.items():
        for alias in aliases:
            column = normalized.get(_normalize_column_name(alias))
            if column is not None:
                canonical[key] = column
                break
    return canonical


def _to_float(value):
    if value is None:
        return float("nan")
    value = str(value).strip()
    if value == "":
        return float("nan")
    try:
        return float(value.replace("D", "e").replace("d", "e"))
    except ValueError:
        return float("nan")


def _as_diagnostics(value):
    if isinstance(value, dict) and "data" in value and "canonical_columns" in value:
        return value
    return read_relaxation_diagnostics(value)


def _resolve_column(diag, name):
    if name is None:
        return None
    if name in diag["data"]:
        return name
    canonical = diag["canonical_columns"].get(str(name))
    if canonical is not None:
        return canonical
    normalized = _normalize_column_name(name)
    for column in diag["columns"]:
        if _normalize_column_name(column) == normalized:
            return column
    return None


def _default_plot_columns(diag):
    sin_theta = diag["canonical_columns"].get("cw_sin_theta")
    if sin_theta is not None:
        return [sin_theta]
    columns = []
    lorentz_force = diag["canonical_columns"].get("lorentz_force")
    if lorentz_force is not None:
        columns.append(lorentz_force)
    if columns:
        return columns
    for key in ("jxb", "divb", "magnetic_energy"):
        column = diag["canonical_columns"].get(key)
        if column is not None and column not in columns:
            columns.append(column)
    if columns:
        return columns
    skip = set(
        column
        for column in (
            diag["canonical_columns"].get("iteration"),
            diag["canonical_columns"].get("time"),
        )
        if column is not None
    )
    return [column for column in diag["columns"] if column not in skip][:4]


def _snapshot_search_location(case_dir, output_dir, base_filename):
    if output_dir is not None:
        output_dir = Path(output_dir).expanduser().resolve()
        return output_dir, Path(base_filename).name
    if case_dir is None:
        return None, None
    base_path = Path(base_filename)
    case_dir = Path(case_dir).expanduser().resolve()
    search_dir = case_dir / base_path.parent
    return search_dir, base_path.name


def _snapshot_number(filename, stem):
    match = _SNAPSHOT_RE.match(filename)
    if match is None:
        return None
    if match.group("stem") != stem:
        return None
    try:
        return int(match.group("number"))
    except ValueError:
        return None


def _infer_snapshot_iterations(files, mf_ditsave, diag):
    if mf_ditsave is None:
        direct = _direct_iteration_map(files, diag)
        if direct:
            return direct, []
        return {}, ["mf_ditsave was not provided; snapshot iterations were not inferred"]

    mf_ditsave = int(mf_ditsave)
    if mf_ditsave <= 0:
        return {}, ["mf_ditsave must be positive; snapshot iterations were not inferred"]

    zero_based = {number: number * mf_ditsave for number, _ in files}
    one_based = {number: (number + 1) * mf_ditsave for number, _ in files}
    if diag is None or diag["canonical_columns"].get("iteration") is None:
        return zero_based, []

    zero_score = _iteration_hit_score(zero_based.values(), diag)
    one_score = _iteration_hit_score(one_based.values(), diag)
    if one_score > zero_score:
        return one_based, []
    if one_score == zero_score and one_score > 0:
        return zero_based, ["zero-based and one-based snapshot iteration mappings both match CSV rows; using zero-based mapping"]
    return zero_based, []


def _direct_iteration_map(files, diag):
    if diag is None:
        return {}
    iteration_column = diag["canonical_columns"].get("iteration")
    if iteration_column is None:
        return {}
    iterations = set(int(round(value)) for value in diag["data"][iteration_column] if math.isfinite(value))
    direct = {}
    for number, _ in files:
        if number in iterations:
            direct[number] = number
    return direct


def _iteration_hit_score(iterations, diag):
    iteration_column = diag["canonical_columns"].get("iteration")
    if iteration_column is None:
        return 0
    csv_iterations = set(
        int(round(value))
        for value in diag["data"][iteration_column]
        if math.isfinite(value)
    )
    return sum(1 for iteration in iterations if int(round(iteration)) in csv_iterations)


def _nearest_iteration_row(diag, iteration):
    iteration_column = diag["canonical_columns"].get("iteration")
    if iteration_column is None:
        return None
    try:
        import numpy as np
    except ImportError as error:  # pragma: no cover
        raise ImportError("numpy is required to match restart snapshots to diagnostics rows") from error

    values = diag["data"][iteration_column]
    finite = np.isfinite(values)
    if not np.any(finite):
        return None
    finite_indices = np.where(finite)[0]
    local = np.argmin(np.abs(values[finite] - float(iteration)))
    return int(finite_indices[local])


def _marker_x_values(markers, x_name):
    values = []
    for marker in markers:
        if x_name is not None and marker.get("values", {}).get(x_name) is not None:
            value = marker["values"][x_name]
        elif x_name is None:
            value = marker.get("row_index")
        else:
            value = marker.get("iteration")
        if value is not None and math.isfinite(float(value)):
            values.append(float(value))
    return values


def _resolve_selection_metric(diag, metric):
    if diag is None:
        return None
    if metric is not None:
        column = _resolve_column(diag, metric)
        if column is None:
            raise ValueError("metric column was not found: {}".format(metric))
        return column
    for key in ("lorentz_force", "jxb", "cw_sin_theta", "divb"):
        column = diag["canonical_columns"].get(key)
        if column is not None:
            return column
    return None


def _marker_metric_value(marker, metric):
    value = marker.get("values", {}).get(metric)
    if value is None:
        return None
    return float(value)


def _latest_marker(markers):
    def key(marker):
        iteration = marker.get("iteration")
        if iteration is None:
            iteration = -1
        return (iteration, marker.get("snapshot_number", -1))

    return max(markers, key=key)


def _diagnostic_label(diag, column):
    labels = {
        "iteration": "MFR iteration",
        "cw_sin_theta": r"$\sin(\theta)$",
        "lorentz_force": "Lorentz force",
        "divb": "f_i (dimensionless divergence)",
        "current": "Current",
    }
    for key, canonical_column in diag["canonical_columns"].items():
        if canonical_column == column and key in labels:
            return labels[key]
    return column


def _select_from_relaxation_trends(diag, markers):
    """Recommend a saved checkpoint from the MFR convergence trends."""

    if diag is None:
        return {
            "selected": None,
            "metric": None,
            "value": None,
            "reason": None,
            "warnings": ["diagnostics are unavailable for automatic trend selection"],
        }

    canonical = diag["canonical_columns"]
    primary = [
        canonical.get("cw_sin_theta"),
        canonical.get("lorentz_force") or canonical.get("jxb"),
    ]
    primary = list(dict.fromkeys(column for column in primary if column is not None))
    if not primary:
        fallback = _resolve_selection_metric(diag, None)
        if fallback is None:
            return {
                "selected": None,
                "metric": None,
                "value": None,
                "reason": None,
                "warnings": ["no sin(theta), Lorentz-force, or divergence column was recognized"],
            }
        primary = [fallback]

    candidates = [marker for marker in markers if _marker_has_metrics(marker, primary)]
    # Iteration zero normally represents the imported potential field before
    # magnetofriction has had a chance to impose and relax the vector boundary.
    noninitial = [
        marker for marker in candidates
        if marker.get("iteration") is None or float(marker["iteration"]) > 0.0
    ]
    warnings = []
    if noninitial:
        candidates = noninitial
    if not candidates:
        return {
            "selected": None,
            "metric": "relaxation_trend_score",
            "value": None,
            "reason": None,
            "warnings": ["no saved checkpoint has finite values for the recognized convergence metrics"],
        }

    divergence = canonical.get("divb")
    stable_candidates = _exclude_divergence_blowup(candidates, divergence)
    if stable_candidates:
        if len(stable_candidates) < len(candidates):
            warnings.append("excluded restart checkpoints whose divergence grew by more than 10x")
        candidates = stable_candidates
    elif divergence is not None:
        warnings.append("all saved checkpoints exceed the divergence stability guard; using unfiltered trends")

    scored = []
    for marker in candidates:
        ranks = [_normalized_rank(marker, candidates, column) for column in primary]
        finite_ranks = [rank for rank in ranks if rank is not None]
        if finite_ranks:
            scored.append((sum(finite_ranks) / len(finite_ranks), marker))
    if not scored:
        return {
            "selected": None,
            "metric": "relaxation_trend_score",
            "value": None,
            "reason": None,
            "warnings": warnings + ["could not compute a finite relaxation trend score"],
        }

    # In an exact tie, prefer the later checkpoint: if all monitored quantities
    # are still falling, this naturally recommends the latest available state.
    selected_value, selected = min(
        scored,
        key=lambda item: (item[0], -_marker_order_value(item[1])),
    )
    if len(primary) >= 2:
        minima = [_minimum_marker_index(candidates, column) for column in primary]
        at_right_edge = all(index == len(candidates) - 1 for index in minima)
        has_interior_minimum = any(0 < index < len(candidates) - 1 for index in minima)
        if at_right_edge:
            reason = "sin(theta) and Lorentz force are still decreasing; latest saved restart recommended"
            warnings.append("no relaxation minimum has appeared yet; consider continuing the relaxation")
        elif has_interior_minimum:
            reason = "best equal-rank compromise near the sin(theta)/Lorentz-force low point"
        else:
            reason = "best equal-rank compromise; sin(theta) and Lorentz force do not yet share a clear minimum"
            warnings.append("sin(theta) and Lorentz force do not show a common low point; inspect the plot before accepting the automatic choice")
    else:
        reason = "minimum saved-checkpoint rank of {}".format(primary[0])
        selected_value = _marker_metric_value(selected, primary[0])

    return {
        "selected": selected,
        "metric": "relaxation_trend_score" if len(primary) >= 2 else primary[0],
        "value": selected_value,
        "reason": reason,
        "warnings": warnings,
    }


def _marker_has_metrics(marker, columns):
    return all(
        column in marker.get("values", {})
        and math.isfinite(float(marker["values"][column]))
        for column in columns
    )


def _exclude_divergence_blowup(markers, column):
    if column is None:
        return list(markers)
    values = [abs(float(marker["values"].get(column, float("nan")))) for marker in markers]
    finite_positive = [value for value in values if math.isfinite(value) and value > 0.0]
    if not finite_positive:
        return list(markers)
    baseline = finite_positive[0]
    limit = 10.0 * baseline
    return [
        marker for marker, value in zip(markers, values)
        if math.isfinite(value) and value <= limit
    ]


def _normalized_rank(marker, markers, column):
    values = [
        abs(float(item["values"][column]))
        for item in markers
        if column in item.get("values", {})
        and math.isfinite(float(item["values"][column]))
    ]
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    target = abs(float(marker["values"][column]))
    lower = sum(1 for value in values if value < target)
    equal = sum(1 for value in values if value == target)
    average_rank = lower + 0.5 * (equal - 1)
    return float(average_rank) / float(len(values) - 1)


def _minimum_marker_index(markers, column):
    return min(
        range(len(markers)),
        key=lambda index: abs(float(markers[index]["values"][column])),
    )


def _marker_order_value(marker):
    iteration = marker.get("iteration")
    if iteration is not None:
        return float(iteration)
    return float(marker.get("snapshot_number", -1))


def _selected_diagnostic_values(diag, marker):
    if diag is None or marker is None:
        return {}
    values = marker.get("values", {})
    selected = {}
    for key in ("divb", "cw_sin_theta", "current", "lorentz_force"):
        column = diag["canonical_columns"].get(key)
        if column is not None and column in values:
            selected[key] = float(values[column])
    return selected
