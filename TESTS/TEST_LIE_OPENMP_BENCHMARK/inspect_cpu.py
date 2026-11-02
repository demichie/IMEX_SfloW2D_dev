#!/usr/bin/env python3
"""Report visible CPUs, process affinity, cgroup CPU quotas and placement."""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess


def read(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return None


def cpu_list_count(value):
    if not value:
        return None
    total = 0
    try:
        for item in value.split(","):
            parts = item.split("-")
            total += int(parts[-1]) - int(parts[0]) + 1
        return total
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    affinity = sorted(os.sched_getaffinity(0))
    cpu_max = read("/sys/fs/cgroup/cpu.max")
    quota = None
    if cpu_max:
        fields = cpu_max.split()
        if len(fields) == 2 and fields[0] != "max":
            quota = int(fields[0]) / int(fields[1])
    else:
        q = read("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        period = read("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        if q and period and int(q) >= 0:
            quota = int(q) / int(period)
    cpuset = read("/sys/fs/cgroup/cpuset.cpus.effective")
    if cpuset is None:
        cpuset = read("/sys/fs/cgroup/cpuset/cpuset.cpus")
    stats = read("/sys/fs/cgroup/cpu.stat")
    try:
        topology = subprocess.run(["lscpu", "-e=CPU,CORE,SOCKET,NODE,ONLINE"],
                                  capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        topology = ""
    result = {
        "machine": platform.platform(),
        "logical_cpus_visible": os.cpu_count(),
        "affinity_cpu_count": len(affinity),
        "affinity_cpu_ids": affinity,
        "cpuset_effective": cpuset,
        "cpuset_cpu_count": cpu_list_count(cpuset),
        "cpu_quota_cores": quota,
        "cpu_stat_initial": stats,
        "lscpu_topology": topology,
        "effective_parallel_capacity_upper_bound": min(
            len(affinity), quota if quota is not None else float("inf")),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    quota_text = "unlimited / not detected" if quota is None else f"{quota:.2f}"
    lines = [
        "## Runner CPU resources",
        "",
        "| Quantity | Observed |",
        "|---|---|",
        f"| Logical CPUs visible | {os.cpu_count()} |",
        f"| CPUs allowed by process affinity | {len(affinity)} |",
        f"| cgroup effective CPU set | {cpuset or 'not reported'} |",
        f"| CPU quota in logical cores | {quota_text} |",
        f"| Upper bound of concurrent CPU capacity | {result['effective_parallel_capacity_upper_bound']:.2f} |",
        "",
        "The capacity estimate is not a physical-core count. SMT, shared",
        "hardware and cgroup throttling can further limit performance.",
        "",
    ]
    args.markdown.write_text("\n".join(lines))
    print("\n".join(lines), flush=True)
    if topology:
        print("CPU topology:\n" + topology, flush=True)
    if stats:
        print("Initial cgroup CPU statistics:\n" + stats, flush=True)


if __name__ == "__main__":
    main()
