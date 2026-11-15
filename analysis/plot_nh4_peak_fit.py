#!/usr/bin/env python3
"""Plot NH4 components at the exact charged x0 values in the configured peak list.

From reproduce/:
    python ../analysis/plot_nh4_peak_fit_fixed_centers.py ../configs/soar.yaml

HDF5: run_root/run_id/data_dir/paths.analysis.data_spectra_h5.
Peak list: paths.core.in_peak_list_dir / paths.core.in_peak_list_filename,
using the same working-directory convention as data_preparation.cpp.
Figures: run_root/run_id/paths.analysis.plot_dir.
Numerical data: plots_dir/paths.analysis.plot_demonstration_data_dir.

Centers are copied from the charged NH4 peak-list entries, never optimized.
The HDF5 supplies measured spectra, calibration, profile and existing labels.
No electron correction is applied to the measured axis or to already charged
peak-list masses. The original HDF5 reference masses are retained for audit.
This constrained demonstration does not independently identify formulas or
recover centroids. Amplitudes, widths and a linear baseline are fitted.
Save and display by default; --no-show saves without opening a window.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import re
import sys

import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import scipy
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
import yaml


SCRIPT_VERSION = "charged-peaklist-centers-v1"
PREFIX = "fig04"
EXTRACT_RANGE = (173.96, 174.30)
FIT_RANGE = (173.99, 174.26)
SEED = 20261112
N_STARTS = 20
ORIGINAL_INPUT_SHA256 = (
    "a0be97e41c326930440e7be485036b6d9b99f1bdfcc52cc7fe43e1689362a154"
)
EXPECTED_FORMULAS = (
    "(NH4)C6H4O5+",
    "(NH4)C7H8O4+",
    "(NH4)C8H12O3+",
    "(NH4)C9H16O2+",
    "(NH4)C10H20O+",
)
# Same monoisotopic masses used to prepare the corrected reference list.
ATOMIC_MASSES = {"C": 12.0, "H": 1.00782503223,
                 "N": 14.00307400443, "O": 15.99491461957}
ELECTRON_MASS = 0.0005485799090441


def require(condition, message):
    """Validation remains active even when Python is invoked with -O."""
    if not condition:
        raise ValueError(message)


def decode(value):
    return value.decode("utf-8").strip() if isinstance(value, bytes) else str(value).strip()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def paths_from_config(config_path):
    """Match analysis.py: relative run_root is relative to the process CWD."""
    with config_path.open("r", encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    require(isinstance(config, dict), "The YAML config must contain a mapping.")
    keys = {
        "run_root": ("paths", "core", "run_root"),
        "in_peak_list_dir": ("paths", "core", "in_peak_list_dir"),
        "in_peak_list_filename": ("paths", "core", "in_peak_list_filename"),
        "run_id": ("run", "run_id"),
        "data_dir": ("paths", "core", "data_dir"),
        "data_spectra_h5": ("paths", "analysis", "data_spectra_h5"),
        "plot_dir": ("paths", "analysis", "plot_dir"),
        "plot_demonstration_data_dir": ("paths", "analysis", "plot_demonstration_data_dir"),
    }
    values = {}
    for name, key_path in keys.items():
        value = config
        for key in key_path:
            require(isinstance(value, dict) and key in value,
                    f"Missing YAML key: {'.'.join(key_path)}")
            value = value[key]
        require(isinstance(value, str) and bool(value.strip()),
                f"YAML key {'.'.join(key_path)} must be a nonempty string.")
        values[name] = value
    run_dir = Path(values["run_root"]) / values["run_id"]
    source = run_dir / values["data_dir"] / values["data_spectra_h5"]
    plots_dir = run_dir / values["plot_dir"]
    data_output_dir = plots_dir / values["plot_demonstration_data_dir"]
    provenance = {
        "config_filename": config_path.name, "config_sha256": sha256(config_path),
        "relative_path_base": "process working directory", "values": values,
    }
    peak_list = Path(values["in_peak_list_dir"]) / values["in_peak_list_filename"]
    return source.resolve(), peak_list.resolve(), plots_dir.resolve(), data_output_dir.resolve(), provenance


def reference_mass(formula):
    """Calculate +1 ion mass for the five explicitly supported NH4 formulas.

    NH4 atoms are already part of the ion composition: subtract one electron,
    without adding another proton. Round to the reference list's six decimals.
    This correction is NOT applied to the measured spectral mass axis or the
    already charged centers copied from the corrected peak list.
    """
    match = re.fullmatch(r"\(NH4\)C(\d+)H(\d+)O(\d*)\+", formula)
    require(match is not None, f"Unsupported reference formula: {formula}")
    carbon, hydrogen, oxygen = match.groups()
    mass = (int(carbon) * ATOMIC_MASSES["C"]
            + (int(hydrogen) + 4) * ATOMIC_MASSES["H"]
            + ATOMIC_MASSES["N"]
            + int(oxygen or "1") * ATOMIC_MASSES["O"] - ELECTRON_MASS)
    return round(mass, 6)


def read_charged_peak_list(path):
    """Read the master list's original x0 values, without rounding or conversion.

    Duplicate NH4 entries are accepted only if their x0 values agree exactly.
    The narrow formula-mass check validates this reference-list demonstration;
    it is not an accuracy threshold for measured experimental centroids.
    """
    matches = {formula: [] for formula in EXPECTED_FORMULAS}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for line, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if len(fields) != 29 or fields[28].strip() != "NH4":
                continue
            formula = fields[2].strip()
            if formula not in matches:
                continue
            mz = float(fields[3])
            reference = reference_mass(formula)
            require(np.isfinite(mz) and abs((mz - reference) / reference) * 1e6 < .1,
                    f"Peak-list line {line}: {formula} x0={mz!r} is not consistent "
                    "with the charged reference mass. Use the charge-corrected master list.")
            matches[formula].append((line, mz))
    entries = {}
    for formula, rows in matches.items():
        require(bool(rows), f"Missing charged NH4 assignment {formula} in {path}. "
                "The plot requires the corrected master peak list, with explicit '+' labels.")
        values = {mz for _, mz in rows}
        require(len(values) == 1, f"Conflicting x0 values for {formula}: {rows}")
        entries[formula] = {"peak_list_mz": rows[0][1],
                            "peak_list_lines_1based": [line for line, _ in rows]}
    return entries


def read_and_align(path):
    """Align every acquisition record using its own stored calibration."""
    with h5py.File(path, "r") as source:
        required = ("FullSpectra/MassAxis", "FullSpectra/SumSpectrum",
                    "FullSpectra/TofData", "FullSpectra/MassCalibration",
                    "Peaks/PeakShape", "PeakData/PeakTable")
        for name in required:
            require(name in source, f"HDF5 dataset is missing: {name}")
        group = source["FullSpectra"]
        equation = decode(group.attrs.get("MassCalibration Function", ""))
        require(equation.replace(" ", "") == "p1*sqrt(m)+p2",
                f"Unsupported calibration function: {equation!r}")
        ion_mode = decode(source.attrs.get("IonMode", ""))
        require(ion_mode.lower() == "positive", "Expected positive ion-mode data.")
        axis = group["MassAxis"][:]
        raw_data = group["TofData"]
        calibration = group["MassCalibration"][:]
        require(axis.ndim == 1 and np.all(np.isfinite(axis))
                and np.all(axis > 0) and np.all(np.diff(axis) > 0),
                "MassAxis must be finite, positive, and strictly increasing.")
        require(raw_data.ndim == 4 and raw_data.shape[1:3] == (1, 1)
                and raw_data.shape[-1] == axis.size,
                f"Expected TofData shape (records, 1, 1, bins), got {raw_data.shape}.")
        n_records = raw_data.shape[0]
        require(n_records > 0 and calibration.shape == (n_records, 2)
                and np.all(np.isfinite(calibration))
                and np.all(calibration[:, 0] > 0), "Invalid per-record calibration.")
        reference_p1 = float(group.attrs["MassCalibration p1"])
        reference_p2 = float(group.attrs["MassCalibration p2"])
        reconstructed_axis = ((np.arange(axis.size) - reference_p2) / reference_p1) ** 2
        require(np.allclose(axis, reconstructed_axis, rtol=1e-12, atol=1e-9),
                "MassAxis disagrees with the reference p1/p2 calibration.")
        region = (axis >= EXTRACT_RANGE[0]) & (axis <= EXTRACT_RANGE[1])
        indices = np.flatnonzero(region)
        require(indices.size > 20, "The requested m/z region is absent or too sparse.")
        x = axis[region]
        aligned = np.empty((n_records, x.size), dtype=np.float64)
        unaligned = np.empty_like(aligned)
        for k, (p1, p2) in enumerate(calibration):
            query = p1 * np.sqrt(x) + p2
            lo = int(np.floor(min(query.min(), indices.min()))) - 2
            hi = int(np.ceil(max(query.max(), indices.max()))) + 3
            require(0 <= lo < hi <= axis.size,
                    f"Calibration requests out-of-range samples in record {k}.")
            values = raw_data[k, 0, 0, lo:hi]
            require(np.all(np.isfinite(values)), f"Nonfinite intensity in record {k}.")
            aligned[k] = np.interp(query, np.arange(lo, hi), values)
            unaligned[k] = values[indices - lo]
        stored_sum = group["SumSpectrum"][indices[0]:indices[-1] + 1]
        unaligned_sum = unaligned.sum(axis=0)
        require(np.allclose(unaligned_sum, stored_sum, rtol=1e-7, atol=.01),
                "The unaligned records do not reproduce the stored SumSpectrum.")

        table = source["PeakData/PeakTable"][:]
        names = ("label", "mass", "lower integration limit", "upper integration limit")
        require(table.dtype.names is not None and all(n in table.dtype.names for n in names),
                "Unrecognized PeakData/PeakTable columns.")
        rows = np.flatnonzero((table["mass"] >= FIT_RANGE[0])
                             & (table["mass"] <= FIT_RANGE[1]))
        rows = rows[np.argsort(table["mass"][rows])]
        assignments = []
        for row in rows:
            entry = table[row]
            native_label = decode(entry["label"])
            label = native_label if native_label.endswith("+") else native_label + "+"
            assignments.append({
                "formula": label, "hdf5_label": native_label,
                "hdf5_peak_table_row_1based": int(row) + 1,
                # float32 -> Python float is exact: no rounding or correction.
                "hdf5_stored_mass": float(entry["mass"]),
                "hdf5_integration_limits": [float(entry[n]) for n in names[2:]],
            })
        require(tuple(a["formula"] for a in assignments) == EXPECTED_FORMULAS,
                "This script targets the original five NH4 assignments near m/z 174; "
                "the input peak table has different entries in that interval.")

        shape = source["Peaks/PeakShape"][:]
        require(shape.ndim == 2 and shape.shape[1] == 2
                and np.all(np.isfinite(shape)) and np.all(np.diff(shape[:, 0]) > 0)
                and np.all(shape[:, 1] >= 0), "Invalid empirical peak-shape profile.")
        sx, sy = shape.T
        peak = int(np.argmax(sy))
        require(0 < peak < len(sy) - 1 and np.isclose(sx[peak], 0, atol=1e-12)
                and np.isclose(sy[peak], 1, atol=1e-12),
                "Expected a unit-height empirical profile centered at zero.")
        require(np.all(np.diff(sy[:peak + 1]) >= 0)
                and np.all(np.diff(sy[peak:]) <= 0)
                and sy[0] < .5 and sy[-1] < .5,
                "Expected a single-peaked profile with both half-height crossings.")
        left = np.interp(.5, sy[:peak + 1], sx[:peak + 1])
        right = np.interp(.5, sy[peak:][::-1], sx[peak:][::-1])
        native_fwhm = float(right - left)
        require(native_fwhm > 0, "Invalid instrument-profile FWHM.")
        shape_x, shape_y = sx / native_fwhm, sy
        warning_path = "FullSpectra/SaturationWarning"
        saturation_count = (int(np.count_nonzero(source[warning_path][:]))
                            if warning_path in source else None)
    return {
        "mz": x, "aligned_sum": aligned.sum(axis=0), "shape_x": shape_x,
        "shape_y": shape_y, "assignments": assignments,
        "provenance": {
            "input_filename": path.name, "input_sha256": sha256(path),
            "original_input_sha256": ORIGINAL_INPUT_SHA256,
            "ion_mode": ion_mode, "records_summed": n_records,
            "calibration_equation": equation,
            "reference_p1": reference_p1, "reference_p2": reference_p2,
            "alignment": "linear interpolation at p1_record*sqrt(reference_mz)+p2_record",
            "smoothing": "none", "record_selection": "all records",
            "unaligned_sum_max_abs_difference": float(np.max(abs(unaligned_sum - stored_sum))),
            "saturation_warning_nonzero_entries": saturation_count,
            "empirical_profile_dataset": "Peaks/PeakShape",
            "empirical_profile_native_fwhm": native_fwhm,
        },
    }


def evaluate(parameters, positions, centers, profile, fit_endpoints):
    """Evaluate profiles at immutable, absolute charged peak-list coordinates.

    Only [amplitude, FWHM] pairs and two baseline endpoints are optimized.
    Keeping centers outside the parameter vector prevents any optimizer from
    changing them. Absolute coordinates also avoid offset/reconstruction steps.
    """
    pars = parameters[:-2].reshape(len(centers), 2)
    components = np.array([
        amplitude * np.nan_to_num(profile((positions - center) / width), nan=0.)
        for (amplitude, width), center in zip(pars, centers)
    ])
    fraction = (positions - fit_endpoints[0]) / (fit_endpoints[1] - fit_endpoints[0])
    baseline = parameters[-2] * (1 - fraction) + parameters[-1] * fraction
    return components.sum(axis=0) + baseline, components, baseline


def fit_spectrum(data):
    mask = (data["mz"] >= FIT_RANGE[0]) & (data["mz"] <= FIT_RANGE[1])
    x, y = data["mz"][mask], data["aligned_sum"][mask]
    scale = float(y.max())
    require(scale > 0, "No positive signal in the fitting interval.")
    v = y / scale
    profile = PchipInterpolator(data["shape_x"], data["shape_y"], extrapolate=False)
    centers = np.array([a["peak_list_mz"] for a in data["assignments"]], dtype=np.float64)
    require(len(centers) == 5 and np.all(np.isfinite(centers))
            and np.all(np.diff(centers) > 0)
            and np.all((centers >= x[0]) & (centers <= x[-1])),
            "Expected five ordered charged peak-list centers inside the fitting interval.")
    centers.setflags(write=False)
    lo = np.r_[np.column_stack([np.zeros(5), np.full(5, .008)]).ravel(), 0., 0.]
    hi = np.r_[np.column_stack([np.full(5, 2.), np.full(5, .045)]).ravel(), .08, .08]
    initial = np.r_[np.column_stack([np.interp(centers, x, v), np.full(5, .021)]).ravel(), .006, .006]
    endpoints = (x[0], x[-1])
    rng = np.random.default_rng(SEED)
    trials = []
    for k in range(N_STARTS):
        start = initial.copy()
        if k:
            p = start[:-2].reshape(5, 2)
            p[:, 0] *= rng.uniform(.75, 1.25, 5)
            p[:, 1] *= rng.uniform(.75, 1.25, 5)
        start = np.clip(start, lo + 1e-9, hi - 1e-9)
        trials.append(least_squares(
            lambda p: evaluate(p, x, centers, profile, endpoints)[0] - v,
            start, bounds=(lo, hi), x_scale="jac", ftol=1e-12, xtol=1e-12,
            gtol=1e-12, max_nfev=2000,
        ))
    best = min(trials, key=lambda trial: np.sum(trial.fun ** 2))
    require(best.success, f"The best optimization did not converge: {best.message}")
    pred, components, baseline = evaluate(best.x, x, centers, profile, endpoints)
    sse = float(np.sum((pred - v) ** 2))
    stable = [trial for trial in trials if np.sum(trial.fun ** 2) <= sse * 1.000001]
    pars = best.x[:-2].reshape(5, 2)
    peaks = []
    for j, assignment in enumerate(data["assignments"]):
        mu = float(centers[j])
        lower, upper = assignment["hdf5_integration_limits"]
        require(lower <= mu <= upper,
                f"Stored peak {j + 1} is outside its native assignment window.")
        require(mu == assignment["peak_list_mz"],
                f"Peak {j + 1} no longer exactly equals its charged peak-list x0.")
        reference = reference_mass(assignment["formula"])
        peaks.append({**assignment, "peak": j + 1, "fixed_center_mz": mu,
                      "center_fitted": False,
                      "center_source": "configured corrected peak list, x0 column",
                      "center_difference_from_peak_list_mz": mu - assignment["peak_list_mz"],
                      "center_difference_from_hdf5_mz": mu - assignment["hdf5_stored_mass"],
                      "center_matches_peak_list_exactly": True,
                      "charge_corrected_reference_mz": reference,
                      "fixed_center_minus_ionic_reference_ppm": 1e6 * (mu - reference) / reference,
                      "amplitude": float(pars[j, 0] * scale),
                      "fwhm": float(pars[j, 1])})
    return {
        "mz": x, "measured": y, "normalized_measured": v, "scale": scale,
        "normalized_fit": pred, "normalized_components": components,
        "normalized_baseline": baseline, "parameters": best.x, "profile": profile,
        "centers": centers, "peaks": peaks,
        "details": {
            "model": "five empirical profiles at fixed charged reference centers plus a linear baseline",
            "objective": "unweighted nonlinear least squares",
            "center_constraint": "exactly fixed to the corrected peak-list x0; not optimized",
            "center_mass_convention": "charged peak-list x0 copied verbatim; no second electron correction",
            "all_centers_match_peak_list_exactly": all(p["center_matches_peak_list_exactly"] for p in peaks),
            "fixed_centers_mz": centers.tolist(),
            "number_of_free_parameters": len(best.x),
            "requested_fit_interval": list(FIT_RANGE),
            "actual_fit_interval": [float(x[0]), float(x[-1])], "bins": len(x),
            "normalization": "divide by maximum measured intensity in fit interval",
            "intensity_scale": scale,
            "parameter_order": "[amplitude, FWHM] for each peak; baseline_left, baseline_right",
            "parameters_normalized": best.x.tolist(),
            "initial_parameters": initial.tolist(),
            "parameter_lower_bounds": lo.tolist(), "parameter_upper_bounds": hi.tolist(),
            "starts": N_STARTS, "random_seed": SEED,
            "optimizer": "scipy.optimize.least_squares; method=trf; x_scale=jac",
            "ftol": 1e-12, "xtol": 1e-12, "gtol": 1e-12, "max_nfev": 2000,
            "optimizer_success": bool(best.success), "optimizer_message": str(best.message),
            "successful_starts": sum(bool(trial.success) for trial in trials),
            "starts_within_1ppm_of_best_sse": len(stable),
            "bound_hits": np.flatnonzero((best.x - lo < 1e-7) | (hi - best.x < 1e-7)).tolist(),
            "sse_normalized": sse,
            "r_squared": 1 - sse / float(np.sum((v - v.mean()) ** 2)),
            "normalized_rmse": float(np.sqrt(np.mean((pred - v) ** 2))),
            "normalized_max_abs_residual": float(np.max(abs(pred - v))),
            "baseline_endpoints": (best.x[-2:] * scale).tolist(),
        },
    }


def formula_mathtext(formula):
    body = re.sub(r"([A-Za-z]+)(\d*)",
                  lambda m: r"\mathrm{" + m[1] + "}" + ("_{" + m[2] + "}" if m[2] else ""),
                  formula.removesuffix("+"))
    return "$" + body + "^+$"


def draw_figure(fit, output):
    x, y = fit["mz"], fit["normalized_measured"]
    # Include each exact center in the display grid, so even the sampled
    # component maxima are located at the charged peak-list x0 values.
    fine = np.unique(np.r_[np.linspace(x[0], x[-1], 2400), fit["centers"]])
    fitted, components, baseline = evaluate(
        fit["parameters"], fine, fit["centers"], fit["profile"], (x[0], x[-1]))
    for j, mu in enumerate(fit["centers"]):
        require(fine[np.argmax(components[j])] == mu,
                f"Displayed component {j + 1} is not centered at its charged peak-list x0.")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelsize": 10, "axes.titlesize": 11, "axes.linewidth": .7,
                         "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    colors = ["#3475B7", "#158C7D", "#B58A16", "#8C61A4", "#CE6B37"]
    raw_color, fit_color, ink = "#BB8A26", "#BE263B", "#222A32"
    fig = plt.figure(figsize=(8.2, 8.2), facecolor="white")
    grid = fig.add_gridspec(2, 2, left=.10, right=.975, top=.945, bottom=.19,
                           height_ratios=[1.0, 1.85], hspace=.42, wspace=.27)
    a = fig.add_subplot(grid[0, 0])
    b = fig.add_subplot(grid[0, 1], sharex=a, sharey=a)
    lower = grid[1, :].subgridspec(2, 1, height_ratios=[4, .8], hspace=.07)
    c = fig.add_subplot(lower[0], sharex=a)
    r = fig.add_subplot(lower[1], sharex=a)
    for ax in (a, b, c, r):
        ax.set_xlim(*FIT_RANGE)
        ax.set_xticks(np.arange(174., 174.251, .05))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(direction="out", length=3, width=.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#E1E5E8", lw=.45, alpha=.8)
        ax.set_axisbelow(True)
    for ax, label, title in ((a, "(a)", "Measured spectrum"),
                              (b, "(b)", "Spectrum and fit"),
                              (c, "(c)", "Resolved peaks with fixed charged reference centers")):
        ax.text(-.052 if ax is c else -.12, 1.08, label, transform=ax.transAxes,
                fontweight="bold", fontsize=12, ha="left", va="bottom", color=ink)
        ax.set_title(title, pad=8, color=ink)
    a.plot(x, y, color=raw_color, lw=1.1, marker=".", markersize=2.)
    a.set_ylim(0, 1.20)
    a.set_yticks([0, .25, .5, .75, 1.])
    a.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    a.set_ylabel("Normalized intensity")
    a.set_xlabel("m/z", labelpad=3)
    b.plot(x, y, color=raw_color, lw=.9, marker=".", markersize=2., label="Measured")
    b.plot(fine, fitted, color=fit_color, lw=1.2, label="Fitted sum")
    b.set_xlabel("m/z", labelpad=3)
    b.legend(loc="upper left", ncol=2, frameon=False, fontsize=8.1, handlelength=1.4,
             handletextpad=.5, columnspacing=1.1, borderaxespad=.2)
    for j, component in enumerate(components):
        c.fill_between(fine, 0, component, color=colors[j], alpha=.09, lw=0)
        c.plot(fine, component, color=colors[j], lw=1.1, zorder=2)
    c.plot(fine, baseline, color="#7B838A", lw=.85, ls="--", zorder=2)
    c.plot(fine, fitted, color=fit_color, lw=1.25, zorder=3)
    c.plot(x, y, "o", mfc="white", mec=ink, mew=.65, ms=2.7, zorder=4)
    offsets = [(-.005, .17), (-.004, .16), (.004, .16), (.004, .18), (.008, .16)]
    annotations = []
    for j, peak in enumerate(fit["peaks"]):
        mu = peak["fixed_center_mz"]
        height = float(np.interp(mu, fine, fitted))
        dx, dy = offsets[j]
        c.vlines(mu, 0, height, color=colors[j], ls=":", lw=.8, zorder=1)
        c.plot(mu, height, "o", color=colors[j], ms=3.7, zorder=5)
        annotations.append(c.annotate(
            f"Peak {j + 1}\n{mu:.6f}", xy=(mu, height), xytext=(mu + dx, height + dy),
            ha="center", va="bottom", fontsize=8.5, color=ink, linespacing=1.4,
            arrowprops={"arrowstyle": "-", "color": colors[j], "lw": .75,
                        "shrinkA": 1.5, "shrinkB": 3},
        ))
    c.set_ylim(0, 1.30)
    c.set_yticks([0, .25, .5, .75, 1., 1.25])
    c.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    c.set_ylabel("Normalized intensity")
    c.tick_params(labelbottom=False)
    c.legend(handles=[
        Line2D([], [], marker="o", ls="None", mfc="white", mec=ink, ms=3, label="Measured"),
        Line2D([], [], color=fit_color, lw=1.25, label="Fitted sum"),
        Line2D([], [], color="#7B838A", ls="--", lw=.85, label="Baseline"),
    ], loc="upper left", ncol=3, frameon=False, fontsize=8.2, handlelength=1.5,
        handletextpad=.55, columnspacing=1.3, borderaxespad=.25)
    r.axhline(0, color="#737F89", lw=.6)
    residual = fit["measured"] - fit["normalized_fit"] * fit["scale"]
    r.plot(x, 100 * residual / fit["scale"], color="#34495D", marker=".", ms=3, lw=.65)
    residual_limit = max(5.2, float(np.ceil(np.max(np.abs(100 * residual / fit["scale"]))) + .2))
    r.set_ylim(-residual_limit, residual_limit)
    residual_tick = float(np.floor(residual_limit - .2))
    r.set_yticks([-residual_tick, 0, residual_tick])
    r.set_ylabel("Residual\n(%)", fontsize=8.7, labelpad=7)
    r.set_xlabel("m/z", labelpad=3)
    handles = [Line2D([], [], color=colors[j], lw=1.8,
                      label=f"{j + 1}: {formula_mathtext(peak['formula'])}")
               for j, peak in enumerate(fit["peaks"])]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.53, .035), ncol=3,
               frameon=False, fontsize=8.7, handlelength=1.4, handletextpad=.65,
               columnspacing=1.75, labelspacing=1.2)
    fig.canvas.draw()
    boxes = [text.get_window_extent(fig.canvas.get_renderer()) for text in annotations]
    require(all(not boxes[i].overlaps(boxes[j]) for i in range(5) for j in range(i)),
            "Peak annotations overlap; review the figure layout.")
    require(all(fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1)
                for box in boxes), "Peak annotation extends beyond the figure.")
    fig.savefig(output / f"{PREFIX}.png", dpi=300, facecolor="white")
    fig.savefig(output / f"{PREFIX}.pdf", facecolor="white", metadata={
        "Title": "Overlapping NH4-adduct peaks near m/z 174", "Author": "",
        "Subject": "Measured spectrum and empirical-profile fit with centers fixed to charged peak-list x0 values",
        "CreationDate": None, "ModDate": None,
    })
    return fig


def write_supporting_data(data, fit, output):
    scale = fit["scale"]
    predicted = fit["normalized_fit"] * scale
    residual = fit["measured"] - predicted
    matrix = np.column_stack([
        fit["mz"], fit["measured"], predicted, fit["normalized_baseline"] * scale,
        (fit["normalized_components"] * scale).T, residual,
        fit["normalized_measured"], fit["normalized_fit"],
    ])
    columns = ["mz", "measured_aligned_sum", "fitted_sum", "baseline"]
    columns += [f"component_{j}" for j in range(1, 6)]
    columns += ["residual_measured_minus_fit", "measured_normalized", "fit_normalized"]
    np.savetxt(output / f"{PREFIX}_fit_data.txt", matrix, delimiter="\t",
               header="\t".join(columns), comments="", fmt="%.17g")
    np.savetxt(output / f"{PREFIX}_instrument_profile.txt",
               np.column_stack([data["shape_x"], data["shape_y"]]), delimiter="\t",
               header="unit_fwhm_coordinate\trelative_intensity", comments="", fmt="%.17g")
    peak_columns = [
        "peak", "formula", "hdf5_label", "fixed_center_mz", "hdf5_stored_mass",
        "center_difference_from_hdf5_mz", "peak_list_mz", "peak_list_lines_1based",
        "center_difference_from_peak_list_mz", "center_matches_peak_list_exactly", "center_fitted",
        "charge_corrected_reference_mz", "fixed_center_minus_ionic_reference_ppm",
        "hdf5_peak_table_row_1based", "fwhm", "amplitude",
    ]
    with (output / f"{PREFIX}_fitted_peaks.txt").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=peak_columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        # Python float strings preserve exact round-trip values. In particular,
        # do not format the fixed center to a limited number of decimal places.
        writer.writerows(fit["peaks"])
    details = {
        "source": data["provenance"], "fit": fit["details"], "peaks": fit["peaks"],
        "script_version": SCRIPT_VERSION,
        "script_filename": Path(__file__).name,
        "script_sha256": sha256(Path(__file__).resolve()),
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                     "h5py": h5py.__version__, "PyYAML": yaml.__version__},
        "reference_mass_calculation": {"charge": 1, "atomic_masses_u": ATOMIC_MASSES,
                                       "electron_mass_u": ELECTRON_MASS,
                                       "reference_rounding_decimals": 6},
        "interpretation": [
            "Formula labels are existing HDF5 assignments; this is not independent chemical identification.",
            "Centers are fixed inputs copied exactly from the corrected peak-list x0 values, not independently estimated centroids.",
            "Original HDF5 reference masses are retained only for audit; the plotted centers use the charged peak-list convention.",
            "The input peak-list masses already include the charge correction. No additional electron correction is applied to them or to the measured spectrum.",
            "Amplitudes, widths and a linear baseline are fitted; this is not a reproduction of saved Igor/Tofware intensities.",
            "Center standard errors are not estimated because centers are fixed; exact copying does not imply zero measurement uncertainty.",
            "Displayed centers are rounded to six decimals; JSON and the peak table retain full round-trip precision.",
            "R-squared measures this spectral fit, not machine-learning performance.",
            "This constrained empirical-profile illustration is not a test of independent centroid detection or formula assignment.",
        ],
    }
    (output / f"{PREFIX}_fit_details.json").write_text(
        json.dumps(details, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify_exported_centers(source, fit, output):
    """Re-read source and both exports, requiring exact equality, not tolerance.

    This also catches accidental decimal truncation in a future output-format
    change. The comparison uses the corrected peak-list numeric values as authority.
    """
    entries = read_charged_peak_list(source)
    expected = np.array([entries[p["formula"]]["peak_list_mz"] for p in fit["peaks"]])
    with (output / f"{PREFIX}_fitted_peaks.txt").open(encoding="utf-8", newline="") as stream:
        text_peaks = list(csv.DictReader(stream, delimiter="\t"))
    details = json.loads((output / f"{PREFIX}_fit_details.json").read_text(encoding="utf-8"))
    values_to_check = {
        "fit": fit["centers"],
        "fit peak records": [p["fixed_center_mz"] for p in fit["peaks"]],
        "text export": [float(p["fixed_center_mz"]) for p in text_peaks],
        "JSON peak records": [p["fixed_center_mz"] for p in details["peaks"]],
        "JSON fit details": details["fit"]["fixed_centers_mz"],
    }
    for name, values in values_to_check.items():
        require(np.array_equal(np.asarray(values, dtype=np.float64), expected),
                f"{name} centers do not exactly match the corrected peak-list x0 values.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    parser.add_argument("config", type=Path, help="Path to the pipeline YAML config")
    parser.add_argument("--no-show", action="store_true",
                        help="Save all outputs without opening the figure (headless runs)")
    args = parser.parse_args()
    try:
        source, peak_list, plots_dir, data_output_dir, path_config = paths_from_config(args.config.expanduser())
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.error(str(error))
    if not source.is_file():
        parser.error(f"Input file does not exist: {source}")
    if not peak_list.is_file():
        parser.error(f"Configured peak list does not exist: {peak_list}")
    if args.no_show:
        matplotlib.use("Agg")
    print(f"Reading and aligning {source}", flush=True)
    entries = read_charged_peak_list(peak_list)
    data = read_and_align(source)
    for assignment in data["assignments"]:
        assignment.update(entries[assignment["formula"]])
    data["provenance"]["peak_list"] = {"filename": peak_list.name, "sha256": sha256(peak_list),
        "formula_column_1based": 3, "x0_column_1based": 4, "reagent_column_1based": 29}
    data["provenance"]["path_configuration"] = path_config
    if data["provenance"]["input_sha256"] != ORIGINAL_INPUT_SHA256:
        print("Note: input bytes differ from the original demonstration HDF5.", file=sys.stderr)
    if data["provenance"]["saturation_warning_nonzero_entries"]:
        print("Note: the input contains saturation-warning entries; inspect the provenance report.", file=sys.stderr)
    print(f"Fitting amplitudes, widths and baseline with {N_STARTS} starts; "
          "centers fixed to exact charged peak-list x0 values...", flush=True)
    fit = fit_spectrum(data)
    plots_dir.mkdir(parents=True, exist_ok=True)
    data_output_dir.mkdir(parents=True, exist_ok=True)
    fig = draw_figure(fit, plots_dir)
    write_supporting_data(data, fit, data_output_dir)
    verify_exported_centers(peak_list, fit, data_output_dir)
    for peak in fit["peaks"]:
        print(f"{peak['formula']:18s}  {peak['fixed_center_mz']!r}  exact charged peak-list match")
    print("Verified: all five centers match the corrected peak list exactly in memory, text and JSON.")
    print(f"R² = {fit['details']['r_squared']:.12f}")
    print(f"Saved PDF: {plots_dir / (PREFIX + '.pdf')}", flush=True)
    print(f"Saved PNG: {plots_dir / (PREFIX + '.png')}", flush=True)
    print(f"Saved numerical outputs to {data_output_dir}", flush=True)
    if not args.no_show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
