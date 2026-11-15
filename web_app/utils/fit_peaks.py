#!/usr/bin/env python3
"""Detect and fit Gaussian peaks in a calibrated mass spectrum.

Drop-in command-line interface (also used by the PyInstaller executable):
    fit_peaks mz_base.csv mz_av.csv output.txt

Inputs are single-column CSVs headed 'mz_base' and 'mz_av', respectively.
Output is ascending fitted m/z centroids, one number per line, without a header.
There is no charge correction: these are measured coordinates, not formula masses.

Detection defaults: raw peak height >= 1500 and prominence >= 1000, both in
the input intensity units. The previous script declared, but did not enforce,
the 1000 threshold. Override with --height / --prominence when appropriate.
No smoothing or intensity resampling is applied. Detection requires an interior
local maximum: an unresolved shoulder or a maximum outside the acquired axis
does not create an extra fitted component.

Overlapping fitting windows are merged. Each group is fitted to the original
samples with positive-amplitude Gaussians plus a free linear baseline. Centers
stay in non-overlapping cells around detected maxima; widths are positive and
bounded. Three deterministic starting widths are tried. Failed/undersampled
fits raise an error instead of silently omitting a group or substituting an
unfitted maximum. A successful no-peaks result writes an empty output file.

Output is replaced atomically only after all fits succeed. On error the command
returns nonzero and leaves any existing output untouched; callers must check
the return code. Normal operation is headless; --plot is optional.

Runtime dependencies: NumPy and SciPy; Matplotlib only for --plot.
Compatible with the existing SoarAI web_app/utils/requirements.txt.
Scientific assumptions: Gaussian component shapes, locally linear background,
and enough sampled, distinguishable maxima. This is not a formula assignment
algorithm or an instrument-specific validation of resolving power.

SciPy references:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks as _detect_peaks


DEFAULT_HEIGHT = 1500.0
DEFAULT_PROMINENCE = 1000.0
_FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
_WINDOW_SIGMAS = 5.0


@dataclass(frozen=True)
class _Candidate:
    index: int
    sigma: float
    prominence: float
    first: int
    stop: int  # exclusive


@dataclass(frozen=True)
class _GroupFit:
    first: int
    stop: int
    centers: np.ndarray
    fitted: np.ndarray


def _read_column(path: str | os.PathLike, header: str) -> np.ndarray:
    """Read strictly one numeric column; accept BOMs and CRLF line endings."""
    with open(path, "r", encoding="utf-8-sig", newline="") as stream:
        rows = csv.reader(stream)
        names = next(rows, None)
        if names is None or len(names) != 1 or names[0].strip() != header:
            raise ValueError(f"{path}: expected one column headed {header!r}")

        def values():
            for row in rows:
                if not row:  # pandas.read_csv also ignores completely blank lines
                    continue
                if len(row) != 1 or not row[0].strip():
                    raise ValueError(f"{path}:{rows.line_num}: expected one numeric value")
                try:
                    value = float(row[0])
                except ValueError as exc:
                    raise ValueError(f"{path}:{rows.line_num}: invalid numeric value") from exc
                if not math.isfinite(value):
                    raise ValueError(f"{path}:{rows.line_num}: NaN and infinity are not allowed")
                yield value

        return np.fromiter(values(), dtype=np.float64)


def _validate_spectrum(mz, intensity) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(mz, dtype=np.float64)
    y = np.asarray(intensity, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
        raise ValueError("m/z and intensity must be one-dimensional arrays of equal length")
    if x.size < 7:
        raise ValueError("at least seven spectrum samples are required")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("spectrum contains NaN or infinity")
    if np.any(x <= 0) or np.any(np.diff(x) <= 0):
        raise ValueError("m/z values must be positive and strictly increasing")
    return x, y


def _candidate_groups(x, y, height, prominence) -> list[list[_Candidate]]:
    # Request widths at half prominence, but do not filter on a guessed width.
    indices, properties = _detect_peaks(
        y, height=height, prominence=prominence, width=(None, None)
    )
    if indices.size == 0:
        return []
    samples = np.arange(x.size)
    left = np.interp(properties["left_ips"], samples, x)
    right = np.interp(properties["right_ips"], samples, x)
    candidates = []
    for j, index in enumerate(indices):
        spacing = min(x[index] - x[index - 1], x[index + 1] - x[index])
        sigma = max((right[j] - left[j]) * _FWHM_TO_SIGMA, spacing)
        radius = _WINDOW_SIGMAS * sigma
        first = max(0, min(int(np.searchsorted(x, x[index] - radius)), int(index) - 4))
        stop = min(x.size, max(int(np.searchsorted(x, x[index] + radius, side="right")), int(index) + 5))
        candidates.append(_Candidate(int(index), sigma, float(properties["prominences"][j]), first, stop))
    # A broad later peak can reach back across several earlier windows. Merge
    # intervals by their left edge, then restore mass order inside each group.
    groups: list[list[_Candidate]] = []
    group_stop = -1
    for candidate in sorted(candidates, key=lambda p: p.first):
        first, stop = candidate.first, candidate.stop
        if not groups or first >= group_stop:
            groups.append([candidate])
            group_stop = stop
        else:
            groups[-1].append(candidate)
            group_stop = max(group_stop, stop)
    return [sorted(group, key=lambda p: p.index) for group in groups]


def _model(parameters: np.ndarray, t: np.ndarray, count: int) -> np.ndarray:
    amplitudes = parameters[:count]
    centers = parameters[count:2 * count]
    sigmas = parameters[2 * count:3 * count]
    z = (t[:, None] - centers) / sigmas
    components = np.exp(-0.5 * z * z) * amplitudes
    return components.sum(axis=1) + parameters[-2] * (1 - t) + parameters[-1] * t


def _jacobian(parameters: np.ndarray, t: np.ndarray, count: int) -> np.ndarray:
    amplitudes = parameters[:count]
    centers = parameters[count:2 * count]
    sigmas = parameters[2 * count:3 * count]
    z = (t[:, None] - centers) / sigmas
    exponentials = np.exp(-0.5 * z * z)
    components = amplitudes * exponentials
    return np.column_stack((
        exponentials,
        components * z / sigmas,
        components * z * z / sigmas,
        1 - t,
        t,
    ))


def _fit_group(x, y, candidates: list[_Candidate]) -> _GroupFit:
    first = min(p.first for p in candidates)
    stop = max(p.stop for p in candidates)
    xx, yy = x[first:stop], y[first:stop]
    count = len(candidates)
    location = f"m/z {xx[0]:.9g}–{xx[-1]:.9g}"
    if xx.size <= 3 * count + 2:
        raise RuntimeError(f"{location}: too few samples to fit {count} peaks and a baseline")

    span = float(xx[-1] - xx[0])
    scale = float(max(np.ptp(yy), max(p.prominence for p in candidates)))
    if not math.isfinite(span) or not math.isfinite(scale) or scale <= 0:
        raise ValueError(f"{location}: invalid mass span or intensity range")
    t = (xx - xx[0]) / span
    edge_count = min(7, max(2, xx.size // 20))
    baseline_left = float(np.median(yy[:edge_count]))
    baseline_right = float(np.median(yy[-edge_count:]))
    offset = 0.5 * baseline_left + 0.5 * baseline_right
    target = (yy - offset) / scale
    peaks = np.array([p.index for p in candidates], dtype=int)
    centers = (x[peaks] - xx[0]) / span
    sigmas = np.array([p.sigma for p in candidates]) / span
    baseline = (baseline_left - offset) / scale, (baseline_right - offset) / scale
    background_at_peaks = baseline_left * (1 - centers) + baseline_right * centers
    amplitudes = np.maximum(y[peaks] - background_at_peaks,
                            [p.prominence for p in candidates]) / scale

    # Disjoint center cells prevent component swaps and multiple fits of one peak.
    midpoints = 0.5 * centers[:-1] + 0.5 * centers[1:]
    center_lower = np.r_[0., midpoints]
    center_upper = np.r_[midpoints, 1.]
    minimum_sigma = np.array([
        0.25 * min(x[p.index] - x[p.index - 1], x[p.index + 1] - x[p.index]) / span
        for p in candidates
    ])
    minimum_sigma = np.maximum(minimum_sigma, 16 * np.finfo(float).eps)
    maximum_sigma = np.full(count, 0.5)
    if np.any(minimum_sigma >= maximum_sigma):
        raise RuntimeError(f"{location}: insufficient mass-axis resolution")
    lower = np.r_[np.zeros(count), center_lower, minimum_sigma, -np.inf, -np.inf]
    upper = np.r_[np.full(count, np.inf), center_upper, maximum_sigma, np.inf, np.inf]

    best = None
    for width_factor in (1.0, 0.7, 1.6):
        initial_widths = np.clip(sigmas * width_factor, minimum_sigma * 1.01, maximum_sigma * 0.99)
        initial = np.r_[amplitudes, centers, initial_widths, baseline]
        result = least_squares(
            lambda p: _model(p, t, count) - target,
            initial,
            jac=lambda p: _jacobian(p, t, count),
            bounds=(lower, upper),
            method="trf",
            loss="linear",
            x_scale="jac",
            ftol=1e-10, xtol=1e-10, gtol=1e-10,
            max_nfev=3000,
        )
        if result.success and np.all(np.isfinite(result.x)) and np.isfinite(result.cost):
            if best is None or result.cost < best.cost:
                best = result
            if result.cost < 1e-24 * xx.size:
                break  # exact noiseless fit; further starts provide no benefit
    if best is None:
        raise RuntimeError(f"{location}: Gaussian fit did not converge; no output was written")

    fitted_centers = xx[0] + best.x[count:2 * count] * span
    if np.any(best.x[:count] <= np.finfo(float).eps):
        raise RuntimeError(f"{location}: a fitted component has negligible amplitude")
    if count > 1 and np.any(np.diff(fitted_centers) <= 16 * np.finfo(float).eps * np.max(xx)):
        raise RuntimeError(f"{location}: fitted centers are indistinguishable")
    return _GroupFit(first, stop, fitted_centers,
                     offset + scale * _model(best.x, t, count))


def _fit_spectrum(mz, intensity, *, height=DEFAULT_HEIGHT, prominence=DEFAULT_PROMINENCE):
    x, y = _validate_spectrum(mz, intensity)
    if not math.isfinite(height) or not math.isfinite(prominence) or prominence < 0:
        raise ValueError("height must be finite; prominence must be finite and nonnegative")
    groups = _candidate_groups(x, y, height, prominence)
    fits = [_fit_group(x, y, group) for group in groups]
    centers = np.sort(np.concatenate([fit.centers for fit in fits])) if fits else np.empty(0)
    return x, y, centers, fits


def fit_centroids(mz, intensity, *, height=DEFAULT_HEIGHT, prominence=DEFAULT_PROMINENCE):
    """Return sorted fitted centroids from arrays; do not read or write files."""
    return _fit_spectrum(mz, intensity, height=height, prominence=prominence)[2]


def _write_centroids(output_path, centers):
    destination = Path(output_path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="ascii", newline="\n", delete=False,
            dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp",
        ) as stream:
            temporary = Path(stream.name)
            for center in centers:
                stream.write(f"{float(center):.17g}\n")
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def find_peaks(mz_b, mz_av, output_path, plot=False, *,
               height=DEFAULT_HEIGHT, prominence=DEFAULT_PROMINENCE):
    """Preserve the original file-based entry point and three-path CLI contract."""
    destination = Path(output_path).resolve()
    if destination in (Path(mz_b).resolve(), Path(mz_av).resolve()):
        raise ValueError("output path must differ from both input paths")
    x = _read_column(mz_b, "mz_base")
    y = _read_column(mz_av, "mz_av")
    x, y, centers, fits = _fit_spectrum(x, y, height=height, prominence=prominence)
    _write_centroids(output_path, centers)
    if plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(x, y, color="0.35", linewidth=0.8, label="Measured")
        for i, fit in enumerate(fits):
            ax.plot(x[fit.first:fit.stop], fit.fitted, color="tab:red",
                    linewidth=1.2, label="Gaussian fit + baseline" if i == 0 else None)
        for center in centers:
            ax.axvline(center, color="tab:blue", alpha=0.35, linewidth=0.7)
        ax.set(xlabel="m/z", ylabel="Intensity", title="Mass spectrum peak fits")
        ax.legend()
        fig.tight_layout()
        plt.show()
        plt.close(fig)
    return centers


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mz_base", help="single-column CSV headed mz_base")
    parser.add_argument("mz_av", help="single-column CSV headed mz_av")
    parser.add_argument("output", help="output text file: one fitted centroid per line")
    parser.add_argument("--height", type=float, default=DEFAULT_HEIGHT,
                        help="minimum raw peak height (default: %(default)g)")
    parser.add_argument("--prominence", type=float, default=DEFAULT_PROMINENCE,
                        help="minimum peak prominence (default: %(default)g)")
    parser.add_argument("--plot", action="store_true", help="display measured and fitted spectra")
    args = parser.parse_args(argv)
    try:
        find_peaks(args.mz_base, args.mz_av, args.output, args.plot,
                   height=args.height, prominence=args.prominence)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        print(f"fit_peaks: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
