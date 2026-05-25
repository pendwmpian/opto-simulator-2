from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class OpticalParams:
    surface_irradiance_mw_mm2: float = 1.0
    mu_eff_mm_inv: float = 2.12
    lateral_sigma0_um: float = 8.0
    lateral_spread_per_depth: float = 0.10
    background_mw_mm2: float = 0.0


def make_pattern(kind: str, n: int = 10) -> np.ndarray:
    pattern = np.zeros((n, n), dtype=float)
    if kind == "single":
        pattern[n // 2, n // 2] = 1.0
    elif kind == "block":
        pattern[4:7, 4:7] = 1.0
    elif kind == "separated":
        for ix, iz in [(2, 2), (7, 2), (2, 7), (7, 7), (5, 5)]:
            pattern[iz, ix] = 1.0
    else:
        raise ValueError(f"Unknown pattern kind: {kind}")
    return pattern


def pixel_centers(n: int, pitch_um: float) -> np.ndarray:
    extent_um = n * pitch_um
    return np.linspace(-extent_um / 2 + pitch_um / 2, extent_um / 2 - pitch_um / 2, n)


def irradiance_xz_at_depth(
    pattern: np.ndarray,
    x_um: np.ndarray,
    z_um: np.ndarray,
    depth_um: float,
    pitch_um: float,
    params: OpticalParams,
) -> np.ndarray:
    """Return irradiance on an x-z plane at a given cortical depth.

    The Gaussian kernel is normalized to preserve total source power, so
    larger depth spreads light laterally while depth attenuation reduces total
    power.
    """
    n_z, n_x = pattern.shape
    centers_x = pixel_centers(n_x, pitch_um)
    centers_z = pixel_centers(n_z, pitch_um)
    xx, zz = np.meshgrid(x_um, z_um, indexing="xy")
    depth_mm = depth_um / 1000.0
    attenuation = math.exp(-params.mu_eff_mm_inv * depth_mm)
    sigma_um = max(1.0, params.lateral_sigma0_um + params.lateral_spread_per_depth * depth_um)
    kernel_norm = pitch_um * pitch_um / (2.0 * math.pi * sigma_um * sigma_um)
    field = np.full_like(xx, params.background_mw_mm2, dtype=float)
    for iz, cz in enumerate(centers_z):
        for ix, cx in enumerate(centers_x):
            amp = pattern[iz, ix]
            if amp == 0:
                continue
            r2 = (xx - cx) ** 2 + (zz - cz) ** 2
            field += (
                amp
                * params.surface_irradiance_mw_mm2
                * attenuation
                * kernel_norm
                * np.exp(-0.5 * r2 / (sigma_um * sigma_um))
            )
    return field


def irradiance_xy_slice(
    pattern: np.ndarray,
    x_um: np.ndarray,
    depth_um: np.ndarray,
    z_slice_um: float,
    pitch_um: float,
    params: OpticalParams,
) -> np.ndarray:
    n_z, n_x = pattern.shape
    centers_x = pixel_centers(n_x, pitch_um)
    centers_z = pixel_centers(n_z, pitch_um)
    xx, dd = np.meshgrid(x_um, depth_um, indexing="xy")
    field = np.full_like(xx, params.background_mw_mm2, dtype=float)
    for depth in depth_um:
        pass
    for iz, cz in enumerate(centers_z):
        for ix, cx in enumerate(centers_x):
            amp = pattern[iz, ix]
            if amp == 0:
                continue
            sigma = np.maximum(1.0, params.lateral_sigma0_um + params.lateral_spread_per_depth * dd)
            attenuation = np.exp(-params.mu_eff_mm_inv * dd / 1000.0)
            r2 = (xx - cx) ** 2 + (z_slice_um - cz) ** 2
            kernel_norm = pitch_um * pitch_um / (2.0 * np.pi * sigma * sigma)
            field += amp * params.surface_irradiance_mw_mm2 * attenuation * kernel_norm * np.exp(-0.5 * r2 / (sigma * sigma))
    return field


def summarize_field(pattern: np.ndarray, pitch_um: float, params: OpticalParams, out_path: Path) -> None:
    x = np.linspace(-160, 160, 161)
    z = np.linspace(-160, 160, 161)
    depths = [0, 100, 300, 600, 900]
    summary = []
    for depth in depths:
        field = irradiance_xz_at_depth(pattern, x, z, depth, pitch_um, params)
        summary.append(
            {
                "depth_um": depth,
                "peak_mw_mm2": float(field.max()),
                "mean_mw_mm2": float(field.mean()),
                "total_relative_power": float(field.sum()),
            }
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(summary, f, indent=2)


def plot_pattern_and_slice(
    pattern: np.ndarray,
    field_xy: np.ndarray,
    x_um: np.ndarray,
    depth_um: np.ndarray,
    pitch_um: float,
    z_slice_um: float,
    params: OpticalParams,
    title: str,
    out_path: Path,
) -> None:
    n_z, n_x = pattern.shape
    extent_x = [-n_x * pitch_um / 2, n_x * pitch_um / 2]
    extent_z = [-n_z * pitch_um / 2, n_z * pitch_um / 2]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    ax = axes[0]
    ax.imshow(
        pattern,
        origin="lower",
        extent=[extent_x[0], extent_x[1], extent_z[0], extent_z[1]],
        cmap="gray_r",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    ax.axhline(z_slice_um, color="#d62728", lw=2, label=f"x-y slice at z={z_slice_um:g} um")
    for edge in np.linspace(extent_x[0], extent_x[1], n_x + 1):
        ax.axvline(edge, color="0.8", lw=0.4)
    for edge in np.linspace(extent_z[0], extent_z[1], n_z + 1):
        ax.axhline(edge, color="0.8", lw=0.4)
    ax.set_xlabel("surface x (um)")
    ax.set_ylabel("surface z (um)")
    ax.set_title("input surface pattern")
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[1]
    positive = field_xy[field_xy > 0]
    norm = matplotlib.colors.LogNorm(vmin=max(float(positive.min()), 1e-5), vmax=max(float(field_xy.max()), 1e-4))
    im = ax.imshow(
        field_xy,
        origin="upper",
        extent=[x_um[0], x_um[-1], depth_um[-1], depth_um[0]],
        aspect="auto",
        cmap="viridis",
        norm=norm,
    )
    ax.set_xlabel("x (um)")
    ax.set_ylabel("depth from cortical surface (um)")
    ax.set_title("irradiance in x-y slice")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("irradiance (mW/mm2)")

    fig.suptitle(
        f"{title}: I0={params.surface_irradiance_mw_mm2:g} mW/mm2, "
        f"mu_eff={params.mu_eff_mm_inv:g} mm^-1, "
        f"sigma(z)= {params.lateral_sigma0_um:g} + {params.lateral_spread_per_depth:g}*depth_um"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/phase3_irradiance"))
    parser.add_argument("--pitch-um", type=float, default=20.0)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--surface-irradiance-mw-mm2", type=float, default=1.0)
    parser.add_argument("--mu-eff-mm-inv", type=float, default=2.12)
    parser.add_argument("--lateral-sigma0-um", type=float, default=8.0)
    parser.add_argument("--lateral-spread-per-depth", type=float, default=0.10)
    parser.add_argument("--z-slice-um", type=float, default=0.0)
    args = parser.parse_args()

    params = OpticalParams(
        surface_irradiance_mw_mm2=args.surface_irradiance_mw_mm2,
        mu_eff_mm_inv=args.mu_eff_mm_inv,
        lateral_sigma0_um=args.lateral_sigma0_um,
        lateral_spread_per_depth=args.lateral_spread_per_depth,
    )
    x = np.linspace(-180, 180, 241)
    depth = np.linspace(0, 1000, 241)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for kind in ["single", "block", "separated"]:
        pattern = make_pattern(kind, args.n)
        field_xy = irradiance_xy_slice(pattern, x, depth, args.z_slice_um, args.pitch_um, params)
        plot_pattern_and_slice(
            pattern=pattern,
            field_xy=field_xy,
            x_um=x,
            depth_um=depth,
            pitch_um=args.pitch_um,
            z_slice_um=args.z_slice_um,
            params=params,
            title=kind,
            out_path=args.out_dir / f"{kind}_pattern_xy_slice.png",
        )
        summarize_field(pattern, args.pitch_um, params, args.out_dir / f"{kind}_depth_summary.json")
        print(f"{kind}: peak={field_xy.max():.4f} mW/mm2 at x-y slice, output={args.out_dir / f'{kind}_pattern_xy_slice.png'}")

    with (args.out_dir / "optical_params.json").open("w") as f:
        json.dump(asdict(params) | {"pitch_um": args.pitch_um, "n": args.n, "z_slice_um": args.z_slice_um}, f, indent=2)


if __name__ == "__main__":
    main()
