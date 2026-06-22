#!/usr/bin/env python3
"""Serve an interactive local viewer for recorded PT5B soma Vm traces."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np


HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>PT5B Vm Viewer</title>
  <style>
    body { font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; color: #202124; }
    h1 { font-size: 22px; margin: 0 0 8px; }
    .meta { color: #5f6368; margin-bottom: 18px; }
    .controls { display: flex; flex-wrap: wrap; gap: 12px; align-items: end; margin-bottom: 14px; }
    label { display: grid; gap: 4px; font-size: 13px; color: #3c4043; }
    select, input, button { font: inherit; padding: 6px 8px; }
    button { cursor: pointer; }
    canvas { width: 100%; height: 460px; border: 1px solid #dadce0; border-radius: 8px; display: block; }
    .stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 8px; margin: 12px 0 18px; }
    .card { border: 1px solid #dadce0; border-radius: 8px; padding: 10px; background: #fff; }
    .card .label { color: #5f6368; font-size: 12px; }
    .card .value { font-size: 18px; margin-top: 2px; }
    table { border-collapse: collapse; width: 100%; margin-top: 18px; font-size: 13px; }
    th, td { border-bottom: 1px solid #e0e0e0; padding: 6px 8px; text-align: right; }
    th:first-child, td:first-child { text-align: left; }
    tbody tr { cursor: pointer; }
    tbody tr:hover { background: #f1f3f4; }
    .hint { color: #5f6368; font-size: 13px; margin-top: 8px; }
  </style>
</head>
<body>
  <h1>PT5B Soma Membrane Potential Viewer</h1>
  <div id="meta" class="meta">Loading...</div>
  <div class="controls">
    <label>Cell GID
      <select id="gidSelect"></select>
    </label>
    <label>Type GID
      <input id="gidInput" type="number" min="0" step="1">
    </label>
    <label>Sort cells
      <select id="sortSelect">
        <option value="gid">GID</option>
        <option value="spikes_desc">Spike count desc</option>
        <option value="mean_desc">Mean Vm desc</option>
        <option value="mean_asc">Mean Vm asc</option>
        <option value="sd_desc">SD desc</option>
        <option value="max_desc">Max Vm desc</option>
      </select>
    </label>
    <label>Max plotted points
      <input id="maxPoints" type="number" min="500" max="20000" step="500" value="5000">
    </label>
    <button id="loadButton">Load Trace</button>
  </div>
  <div id="stats" class="stats"></div>
  <canvas id="plot" width="1400" height="520"></canvas>
  <div class="hint">Click a row below to load that cell. Traces are downsampled on demand from the NPZ file.</div>
  <table>
    <thead>
      <tr><th>GID</th><th>Spikes</th><th>Mean Vm</th><th>SD</th><th>Min Vm</th><th>Max Vm</th></tr>
    </thead>
    <tbody id="cellTable"></tbody>
  </table>

<script>
let cells = [];
let meta = {};

function fmt(x, digits = 2) {
  return Number.isFinite(x) ? x.toFixed(digits) : "n/a";
}

function sortedCells() {
  const mode = document.getElementById("sortSelect").value;
  const arr = [...cells];
  if (mode === "spikes_desc") arr.sort((a, b) => b.spike_count - a.spike_count || a.gid - b.gid);
  else if (mode === "mean_desc") arr.sort((a, b) => b.mean_mV - a.mean_mV || a.gid - b.gid);
  else if (mode === "mean_asc") arr.sort((a, b) => a.mean_mV - b.mean_mV || a.gid - b.gid);
  else if (mode === "sd_desc") arr.sort((a, b) => b.sd_mV - a.sd_mV || a.gid - b.gid);
  else if (mode === "max_desc") arr.sort((a, b) => b.max_mV - a.max_mV || a.gid - b.gid);
  else arr.sort((a, b) => a.gid - b.gid);
  return arr;
}

function rebuildCellLists() {
  const select = document.getElementById("gidSelect");
  const table = document.getElementById("cellTable");
  select.innerHTML = "";
  table.innerHTML = "";
  for (const cell of sortedCells()) {
    const opt = document.createElement("option");
    opt.value = cell.gid;
    opt.textContent = `${cell.gid} | spikes ${cell.spike_count} | mean ${fmt(cell.mean_mV, 1)} mV`;
    select.appendChild(opt);

    const tr = document.createElement("tr");
    tr.innerHTML = `<td>${cell.gid}</td><td>${cell.spike_count}</td><td>${fmt(cell.mean_mV)}</td><td>${fmt(cell.sd_mV)}</td><td>${fmt(cell.min_mV)}</td><td>${fmt(cell.max_mV)}</td>`;
    tr.addEventListener("click", () => {
      select.value = String(cell.gid);
      document.getElementById("gidInput").value = String(cell.gid);
      loadTrace(cell.gid);
    });
    table.appendChild(tr);
  }
}

function drawTrace(trace) {
  const canvas = document.getElementById("plot");
  const ctx = canvas.getContext("2d");
  const w = canvas.width, h = canvas.height;
  ctx.clearRect(0, 0, w, h);
  const pad = {left: 70, right: 20, top: 28, bottom: 48};
  const xs = trace.t_ms, ys = trace.v_mV;
  const xmin = xs[0], xmax = xs[xs.length - 1];
  const ymin = Math.min(-80, Math.min(...ys));
  const ymax = Math.max(30, Math.max(...ys));
  const xmap = x => pad.left + (x - xmin) / (xmax - xmin) * (w - pad.left - pad.right);
  const ymap = y => h - pad.bottom - (y - ymin) / (ymax - ymin) * (h - pad.top - pad.bottom);

  ctx.strokeStyle = "#e0e0e0";
  ctx.lineWidth = 1;
  ctx.font = "13px system-ui";
  ctx.fillStyle = "#5f6368";
  for (let i = 0; i <= 6; i++) {
    const yv = ymin + i * (ymax - ymin) / 6;
    const y = ymap(yv);
    ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(w - pad.right, y); ctx.stroke();
    ctx.fillText(`${fmt(yv, 0)} mV`, 8, y + 4);
  }
  for (let i = 0; i <= 6; i++) {
    const xv = xmin + i * (xmax - xmin) / 6;
    const x = xmap(xv);
    ctx.beginPath(); ctx.moveTo(x, h - pad.bottom); ctx.lineTo(x, h - pad.bottom + 5); ctx.stroke();
    ctx.fillText(`${fmt(xv, 0)} ms`, x - 22, h - 18);
  }

  ctx.strokeStyle = "#1769aa";
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  for (let i = 0; i < xs.length; i++) {
    const x = xmap(xs[i]), y = ymap(ys[i]);
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.strokeStyle = "#c62828";
  ctx.setLineDash([8, 5]);
  const meanY = ymap(trace.stats.mean_mV);
  ctx.beginPath(); ctx.moveTo(pad.left, meanY); ctx.lineTo(w - pad.right, meanY); ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "#202124";
  ctx.font = "16px system-ui";
  ctx.fillText(`GID ${trace.gid} soma Vm`, pad.left, 20);
}

function renderStats(trace) {
  const s = trace.stats;
  const items = [
    ["Spikes", s.spike_count],
    ["Mean Vm", `${fmt(s.mean_mV)} mV`],
    ["SD", `${fmt(s.sd_mV)} mV`],
    ["Min Vm", `${fmt(s.min_mV)} mV`],
    ["Max Vm", `${fmt(s.max_mV)} mV`],
    ["Returned points", trace.t_ms.length]
  ];
  document.getElementById("stats").innerHTML = items.map(([label, value]) =>
    `<div class="card"><div class="label">${label}</div><div class="value">${value}</div></div>`
  ).join("");
}

async function loadTrace(gidArg) {
  const gid = gidArg || document.getElementById("gidInput").value || document.getElementById("gidSelect").value;
  const maxPoints = document.getElementById("maxPoints").value;
  const res = await fetch(`/api/trace?gid=${encodeURIComponent(gid)}&max_points=${encodeURIComponent(maxPoints)}`);
  if (!res.ok) {
    alert(await res.text());
    return;
  }
  const trace = await res.json();
  document.getElementById("gidSelect").value = String(trace.gid);
  document.getElementById("gidInput").value = String(trace.gid);
  renderStats(trace);
  drawTrace(trace);
}

async function main() {
  const res = await fetch("/api/summary");
  const data = await res.json();
  cells = data.cells;
  meta = data.meta;
  document.getElementById("meta").textContent =
    `condition=${meta.condition}, trial=${meta.trial}, ihGbar=${meta.ihGbar}, cells=${data.cell_count}, duration=${data.duration_ms} ms, analysis starts at ${data.analysis_start_ms} ms`;
  rebuildCellLists();
  const target = meta.target_gid || cells[0].gid;
  document.getElementById("gidSelect").value = String(target);
  document.getElementById("gidInput").value = String(target);
  await loadTrace(target);
}

document.getElementById("sortSelect").addEventListener("change", rebuildCellLists);
document.getElementById("loadButton").addEventListener("click", () => loadTrace());
document.getElementById("gidSelect").addEventListener("change", e => {
  document.getElementById("gidInput").value = e.target.value;
  loadTrace(e.target.value);
});
main();
</script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="Run directory containing pt5b_vm.npz")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--analysis-start-ms", type=float, default=1000.0)
    return parser.parse_args()


class VmStore:
    def __init__(self, run_dir: Path, analysis_start_ms: float) -> None:
        self.run_dir = run_dir
        self.analysis_start_ms = analysis_start_ms
        with np.load(run_dir / "pt5b_vm.npz", allow_pickle=False) as data:
            self.t_ms = np.asarray(data["t_ms"], dtype=float)
            self.gids = np.asarray(data["gids"], dtype=int)
            self.v_mV = np.asarray(data["v_mV"], dtype=float)
        manifest_path = run_dir / "manifest.json"
        self.manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        self.gid_to_index = {int(gid): i for i, gid in enumerate(self.gids)}
        self.spike_counts = self._load_spike_counts()
        self.summary = self._compute_summary()

    def _load_spike_counts(self) -> np.ndarray:
        spikes_path = self.run_dir / "spikes.npz"
        if not spikes_path.exists():
            return np.zeros(self.gids.shape, dtype=int)
        with np.load(spikes_path, allow_pickle=False) as spikes:
            counts = Counter(spikes["spkid"].astype(int))
        return np.asarray([counts.get(int(gid), 0) for gid in self.gids], dtype=int)

    def _compute_summary(self) -> list[dict[str, float | int]]:
        mask = self.t_ms >= self.analysis_start_ms
        if not np.any(mask):
            raise ValueError("analysis start is later than the recorded trace")
        post_v = self.v_mV[:, mask]
        mean_v = post_v.mean(axis=1)
        sd_v = post_v.std(axis=1)
        min_v = self.v_mV.min(axis=1)
        max_v = self.v_mV.max(axis=1)
        return [
            {
                "gid": int(gid),
                "spike_count": int(spikes),
                "mean_mV": float(mean),
                "sd_mV": float(sd),
                "min_mV": float(minimum),
                "max_mV": float(maximum),
            }
            for gid, spikes, mean, sd, minimum, maximum in zip(
                self.gids, self.spike_counts, mean_v, sd_v, min_v, max_v
            )
        ]

    def trace_payload(self, gid: int, max_points: int) -> dict[str, object]:
        if gid not in self.gid_to_index:
            raise KeyError(f"GID {gid} was not recorded")
        idx = self.gid_to_index[gid]
        t_ds, v_ds = minmax_downsample(self.t_ms, self.v_mV[idx], max_points)
        return {
            "gid": gid,
            "t_ms": t_ds.tolist(),
            "v_mV": v_ds.tolist(),
            "stats": self.summary[idx],
        }


def minmax_downsample(t_ms: np.ndarray, v_mV: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    max_points = max(100, int(max_points))
    if t_ms.size <= max_points:
        return t_ms, v_mV
    bins = max(1, max_points // 2)
    edges = np.linspace(0, t_ms.size, bins + 1, dtype=int)
    out_t: list[float] = []
    out_v: list[float] = []
    for start, stop in zip(edges[:-1], edges[1:]):
        if stop <= start:
            continue
        segment = v_mV[start:stop]
        local_min = int(np.argmin(segment)) + start
        local_max = int(np.argmax(segment)) + start
        for idx in sorted({local_min, local_max}):
            out_t.append(float(t_ms[idx]))
            out_v.append(float(v_mV[idx]))
    return np.asarray(out_t), np.asarray(out_v)


def make_handler(store: VmStore) -> type[BaseHTTPRequestHandler]:
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: object) -> None:
            return

        def send_json(self, payload: object, status: HTTPStatus = HTTPStatus.OK) -> None:
            body = json.dumps(payload, separators=(",", ":")).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def send_text(self, text: str, status: HTTPStatus = HTTPStatus.OK, content_type: str = "text/plain") -> None:
            body = text.encode()
            self.send_response(status)
            self.send_header("Content-Type", f"{content_type}; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path == "/":
                self.send_text(HTML, content_type="text/html")
                return
            if parsed.path == "/api/summary":
                dt = float(np.median(np.diff(store.t_ms))) if store.t_ms.size > 1 else 0.0
                self.send_json(
                    {
                        "meta": store.manifest,
                        "cell_count": int(store.gids.size),
                        "duration_ms": float(store.t_ms[-1] - store.t_ms[0] + dt),
                        "analysis_start_ms": store.analysis_start_ms,
                        "cells": store.summary,
                    }
                )
                return
            if parsed.path == "/api/trace":
                query = parse_qs(parsed.query)
                try:
                    gid = int(query.get("gid", [store.gids[0]])[0])
                    max_points = int(query.get("max_points", [5000])[0])
                    self.send_json(store.trace_payload(gid, max_points))
                except (ValueError, KeyError) as exc:
                    self.send_text(str(exc), status=HTTPStatus.BAD_REQUEST)
                return
            self.send_text("not found", status=HTTPStatus.NOT_FOUND)

    return Handler


def main() -> None:
    args = parse_args()
    store = VmStore(args.run_dir, args.analysis_start_ms)
    server = ThreadingHTTPServer((args.host, args.port), make_handler(store))
    print(f"serving PT5B Vm viewer at http://{args.host}:{args.port}")
    print(f"run_dir={args.run_dir}")
    server.serve_forever()


if __name__ == "__main__":
    main()
