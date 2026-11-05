#!/usr/bin/env bash
# Prefer the runner's existing APT indexes; refresh only on install failure.
# This avoids a separate network-dependent apt-get update in every matrix job.
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
    -o Acquire::Retries=3
    -o Acquire::http::Timeout=20
    -o Acquire::https::Timeout=20
    -o Dpkg::Use-Pty=0
)
apt_install() {
    sudo env DEBIAN_FRONTEND=noninteractive timeout --signal=TERM --kill-after=10s 300s \
        apt-get "${apt_options[@]}" install --yes --no-install-recommends "$@"
}
echo 'APT_INSTALL_FIRST attempting with runner package indexes'
if ! apt_install "${missing[@]}"; then
    echo 'APT_INSTALL_RETRY first attempt failed; refreshing indexes'
    # APT update can return success with some failed indexes. Require a
    # successful installation afterwards rather than treating update as proof.
    sudo env DEBIAN_FRONTEND=noninteractive timeout --signal=TERM --kill-after=10s 180s \
        apt-get "${apt_options[@]}" update
    apt_install "${missing[@]}"
fi
for package in "${missing[@]}"; do
    dpkg-query -W -f='${Status}' "$package" | grep -Fxq 'install ok installed'
done
echo 'APT_INSTALLED all requested dependencies are available'
