#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
compiler="${FC:-gfortran}"
workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

flags=(-std=f2008 -O0 -g -fopenmp -fcheck=all -fbacktrace -Wall -Wextra -Werror)
if [[ "${1:-}" == --benchmark || "${1:-}" == --benchmark-large ]]; then
    # Production-like optimization for timings; correctness runs first in the same executable.
    flags=(-std=f2008 -O3 -fopenmp -Wall -Wextra -Werror)
fi

"$compiler" "${flags[@]}" -J"$workdir" \
    compaction_kernels.f90 test_domain_compaction.f90 -o "$workdir/test_domain_compaction"

"$workdir/test_domain_compaction" "$@"

"$compiler" "${flags[@]}" -J"$workdir" \
    test_domain_halo.f90 -o "$workdir/test_domain_halo"
"$workdir/test_domain_halo"
