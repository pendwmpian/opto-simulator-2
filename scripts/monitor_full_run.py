#!/usr/bin/env python3
"""Sample aggregate memory for a running full-model job."""

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
import time
import psutil


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--match", default="run_full_quiet_pt5b.py")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval-seconds", type=float, default=30.0)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["time_utc", "process_count", "total_rss_bytes", "max_rss_bytes", "system_available_bytes", "system_used_percent"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        while True:
            rss = []
            for process in psutil.process_iter(["cmdline", "memory_info"]):
                try:
                    command = " ".join(process.info["cmdline"] or [])
                    if args.match in command and process.pid != psutil.Process().pid:
                        rss.append(process.info["memory_info"].rss)
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            memory = psutil.virtual_memory()
            writer.writerow({"time_utc": datetime.now(timezone.utc).isoformat(), "process_count": len(rss), "total_rss_bytes": sum(rss), "max_rss_bytes": max(rss, default=0), "system_available_bytes": memory.available, "system_used_percent": memory.percent})
            handle.flush()
            if not rss:
                break
            time.sleep(args.interval_seconds)


if __name__ == "__main__":
    main()
