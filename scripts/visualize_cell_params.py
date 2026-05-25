from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection


def load_pickle(path: Path):
    with path.open("rb") as f:
        try:
            return pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)
            return pickle.load(f, encoding="latin1")


def sec_color(name: str) -> str:
    if name.startswith("soma"):
        return "#111111"
    if name.startswith("axon"):
        return "#2ca02c"
    if name.startswith("apic") or name.startswith("Adend"):
        return "#d62728"
    if name.startswith("dend") or name.startswith("Bdend"):
        return "#1f77b4"
    return "#7f7f7f"


def collect_segments(rule):
    segments = []
    colors = []
    widths = []
    for name, sec in rule.get("secs", {}).items():
        pts = sec.get("geom", {}).get("pt3d", [])
        if len(pts) < 2:
            continue
        xyz = [(float(p[0]), float(p[1]), float(p[2])) for p in pts]
        diam = [float(p[3]) for p in pts]
        for i in range(len(xyz) - 1):
            segments.append([xyz[i], xyz[i + 1]])
            colors.append(sec_color(name))
            widths.append(max(0.3, min(4.0, (diam[i] + diam[i + 1]) * 0.35)))
    return segments, colors, widths


def set_equal_3d(ax, segments):
    xs = [p[0] for seg in segments for p in seg]
    ys = [p[1] for seg in segments for p in seg]
    zs = [p[2] for seg in segments for p in seg]
    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    cz = (min(zs) + max(zs)) / 2
    radius = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)) / 2
    ax.set_xlim(cx - radius, cx + radius)
    ax.set_ylim(cy - radius, cy + radius)
    ax.set_zlim(cz - radius, cz + radius)


def plot_cell(path: Path, out: Path, title: str | None = None) -> None:
    rule = load_pickle(path)
    segments, colors, widths = collect_segments(rule)

    fig = plt.figure(figsize=(13, 6), constrained_layout=True)
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)

    lc3d = Line3DCollection(segments, colors=colors, linewidths=widths, alpha=0.95)
    ax3d.add_collection3d(lc3d)
    set_equal_3d(ax3d, segments)
    ax3d.set_xlabel("x (um)")
    ax3d.set_ylabel("y (um)")
    ax3d.set_zlabel("z (um)")
    ax3d.view_init(elev=14, azim=-68)
    ax3d.set_title("3D morphology")

    for seg, color, width in zip(segments, colors, widths):
        ax2d.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], color=color, lw=width)
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.set_xlabel("x (um)")
    ax2d.set_ylabel("y (um)")
    ax2d.invert_yaxis()
    ax2d.set_title("x-y projection")

    fig.suptitle(title or path.name)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cell_params", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--title")
    args = parser.parse_args()
    plot_cell(args.cell_params, args.out, args.title)


if __name__ == "__main__":
    main()
