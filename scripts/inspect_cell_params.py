from __future__ import annotations

import argparse
import pickle
from pathlib import Path


def load_pickle(path: Path):
    with path.open("rb") as f:
        try:
            return pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)
            return pickle.load(f, encoding="latin1")


def summarize_cell_rule(path: Path) -> None:
    rule = load_pickle(path)
    secs = rule.get("secs", {})
    sec_lists = rule.get("secLists", {})

    n3d = 0
    npts3d = 0
    mech_names: set[str] = set()
    sec_type_counts: dict[str, int] = {}
    total_l = 0.0
    xyz_by_type: dict[str, list[tuple[float, float, float]]] = {}

    for name, sec in secs.items():
        prefix = "".join(ch for ch in name if not ch.isdigit() and ch not in "[]_")
        sec_type_counts[prefix or "other"] = sec_type_counts.get(prefix or "other", 0) + 1

        geom = sec.get("geom", {})
        if "pt3d" in geom:
            n3d += 1
            npts3d += len(geom["pt3d"])
            xyz_by_type.setdefault(prefix or "other", []).extend(
                (float(p[0]), float(p[1]), float(p[2])) for p in geom["pt3d"]
            )
        total_l += float(geom.get("L", 0.0) or 0.0)

        for mech_name in sec.get("mechs", {}).keys():
            mech_names.add(mech_name)

    print(f"file: {path}")
    print(f"top_level_keys: {sorted(rule.keys())}")
    print(f"sections: {len(secs)}")
    print(f"section_type_counts: {dict(sorted(sec_type_counts.items()))}")
    print(f"secLists: {sorted(sec_lists.keys())}")
    print(f"sections_with_pt3d: {n3d}")
    print(f"pt3d_points: {npts3d}")
    print(f"sum_section_L_um: {total_l:.2f}")
    print("pt3d_bbox_by_type:")
    for sec_type, pts in sorted(xyz_by_type.items()):
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        zs = [p[2] for p in pts]
        print(
            f"  {sec_type}: "
            f"x=[{min(xs):.1f}, {max(xs):.1f}], "
            f"y=[{min(ys):.1f}, {max(ys):.1f}], "
            f"z=[{min(zs):.1f}, {max(zs):.1f}]"
        )
    print(f"mechanisms: {sorted(mech_names)}")

    example_names = list(secs.keys())[:5]
    print("example_sections:")
    for name in example_names:
        sec = secs[name]
        geom = sec.get("geom", {})
        print(
            f"  {name}: geom_keys={sorted(geom.keys())}, "
            f"mechs={sorted(sec.get('mechs', {}).keys())}"
        )
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()

    for path in args.files:
        summarize_cell_rule(path)


if __name__ == "__main__":
    main()
