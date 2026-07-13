#!/usr/bin/env python3
"""Build the upstream M1 mechanisms for both NEURON and CoreNEURON."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_UPSTREAM = ROOT / "external" / "M1_NetPyNE_CellReports_2023"
EXCLUDED = {"dipole.mod", "dipole_pp.mod", "gabab.mod", "my_exp2syn.mod"}
STDLIB_VERBATIM_FILES = {
    "ch_CavL.mod", "ch_CavN.mod", "ch_KCaS.mod", "ch_Kdrfastngf.mod",
    "ch_KvAngf.mod", "ch_KvCaB.mod", "ch_leak.mod", "iconc_Ca.mod",
}
REPLACEMENTS = {
    "Nca.mod": [
        ("GLOBAL q10, temp, tadj, vmin, vmax, vshift", "GLOBAL q10, temp, vmin, vmax, vshift\n  RANGE tadj"),
    ],
    "cancr.mod": [("GLOBAL hinf, minf, s_inf", "RANGE hinf, minf, s_inf")],
    "catcb.mod": [("GLOBAL hinf, minf", "RANGE hinf, minf")],
    "kapcb.mod": [
        ("RANGE gkabar, ik", "RANGE gkabar, ik, qt"),
        ("        taun            (ms)", "        taun            (ms)\n        qt              (1)"),
        ("LOCAL qt", ""),
    ],
    "kapin.mod": [
        ("RANGE gkabar, ik", "RANGE gkabar, ik, qt"),
        ("        taun            (ms)", "        taun            (ms)\n        qt              (1)"),
        ("LOCAL qt", ""),
    ],
    "ch_CavL.mod": [("GLOBAL minf,mtau", "RANGE minf,mtau")],
    "ch_KvAngf.mod": [("GLOBAL ninf, linf, taul, taun", "RANGE ninf, linf, taul, taun")],
    "ch_KvCaB.mod": [("GLOBAL oinf, otau", "RANGE oinf, otau")],
    "naz.mod": [
        ("GLOBAL minf, hinf, mtau, htau", "RANGE minf, hinf, mtau, htau"),
        ("GLOBAL q10, temp, tadj, vmin, vmax, vshift", "GLOBAL q10, temp, vmin, vmax, vshift\n  RANGE tadj"),
    ],
}


def replace_once(text: str, old: str, new: str, path: Path) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"Expected one occurrence in {path}, found {count}: {old!r}")
    return text.replace(old, new)


def prepare_sources(upstream: Path, destination: Path) -> None:
    source = upstream / "sim" / "mod"
    if not source.is_dir():
        raise FileNotFoundError(f"Upstream MOD directory not found: {source}")
    shutil.rmtree(destination, ignore_errors=True)
    shutil.copytree(source, destination)
    for name in EXCLUDED:
        (destination / name).unlink()

    stdlib_block = re.compile(r"\nVERBATIM\n#include <stdlib\.h>.*?ENDVERBATIM\n", re.DOTALL)
    for name in STDLIB_VERBATIM_FILES:
        path = destination / name
        text = path.read_text(encoding="utf-8").replace("\r\n", "\n")
        text, count = stdlib_block.subn("\n", text, count=1)
        if count != 1:
            raise RuntimeError(f"Expected legacy stdlib VERBATIM block in {path}")
        path.write_text(text, encoding="utf-8")

    cavl = destination / "ch_CavL.mod"
    text = cavl.read_text(encoding="utf-8")
    text, count = re.subn(r"\n\s*VERBATIM\n\s*cai=_ion_cai;\n\s*ENDVERBATIM", "", text, count=1)
    if count != 1:
        raise RuntimeError(f"Expected explicit cai update in {cavl}")
    cavl.write_text(text, encoding="utf-8")

    for name, changes in REPLACEMENTS.items():
        path = destination / name
        text = path.read_text(encoding="utf-8")
        for old, new in changes:
            text = replace_once(text, old, new, path)
        path.write_text(text, encoding="utf-8")

    shutil.copy2(ROOT / "coreneuron" / "vecstim.mod", destination / "vecstim.mod")
    for custom_mod in (ROOT / "mod").glob("*.mod"):
        shutil.copy2(custom_mod, destination / custom_mod.name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream-dir", type=Path, default=DEFAULT_UPSTREAM)
    parser.add_argument("--prepared-dir", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    upstream = args.upstream_dir.resolve()
    sim_dir = upstream / "sim"
    prepared = (args.prepared_dir or sim_dir / ".coreneuron-mod").resolve()
    prepare_sources(upstream, prepared)
    print(f"Prepared CoreNEURON-compatible MOD sources in {prepared}", flush=True)
    if not args.prepare_only:
        subprocess.run(["nrnivmodl", "-coreneuron", str(prepared)], cwd=sim_dir, check=True)


if __name__ == "__main__":
    main()
