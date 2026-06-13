#!/usr/bin/env python3
"""Create a compact Markdown summary from a completed quiet run."""

import argparse
import json
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    manifest = json.loads((args.run_dir / "manifest.json").read_text())
    spikes = np.load(args.run_dir / "spikes.npz")
    vm = np.load(args.run_dir / "target_vm.npz")
    target_gid = int(vm["target_gid"][0])
    target_spikes = int(np.count_nonzero(spikes["spkid"] == target_gid))
    voltage = vm["v_mV"]
    lines = ["# Quiet PT5B Run Summary", "", f"- Trial: `{manifest['trial']}`", f"- Duration: `{manifest['duration_ms']} ms`", f"- MPI ranks: `{manifest['nhost']}`", f"- Wall time: `{manifest['elapsed_seconds']:.1f} s`", f"- Network spikes: `{spikes['spkt'].size}`", f"- Target PT5B gid: `{target_gid}`", f"- Target spikes: `{target_spikes}`", f"- Target Vm samples: `{voltage.size}`"]
    if voltage.size:
        lines += [f"- Target Vm mean: `{float(np.mean(voltage)):.3f} mV`", f"- Target Vm SD: `{float(np.std(voltage)):.3f} mV`", f"- Target Vm range: `{float(np.min(voltage)):.3f}` to `{float(np.max(voltage)):.3f} mV`"]
    (args.run_dir / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
