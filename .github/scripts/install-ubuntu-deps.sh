#!/usr/bin/env bash
# Skip APT when packages are present; bound mirror/network waits otherwise.
set -euo pipefail
if (( $# == 0 )); then
    echo "usage: $0 package [package ...]" >&2
    exit 2
fi
missing=()
for package in "$@"; do
    if dpkg-query -W -f='${Status}' "$package" 2>/dev/null | grep -Fxq 'install ok installed'; then
        printf 'APT_CACHE_HIT %s\n' "$package"
    else
        missing+=("$package")
    fi
done
if (( ${#missing[@]} == 0 )); then
    echo 'APT_SKIPPED all requested packages are installed'
    exit 0
fi
printf 'APT_MISSING %s\n' "${missing[*]}"
apt_options=(
    -o Acquire::Retries=2
    -o Acquire::http::Timeout=15
    -o Acquire::https::Timeout=15
    -o Dpkg::Use-Pty=0
)
sudo env DEBIAN_FRONTEND=noninteractive timeout --signal=TERM --kill-after=10s 150s \
    apt-get "${apt_options[@]}" update
sudo env DEBIAN_FRONTEND=noninteractive timeout --signal=TERM --kill-after=10s 240s \
    apt-get "${apt_options[@]}" install --yes --no-install-recommends "${missing[@]}"
for package in "${missing[@]}"; do
    dpkg-query -W -f='${Status}' "$package" | grep -Fxq 'install ok installed'
done
echo 'APT_INSTALLED all missing dependencies are available'
