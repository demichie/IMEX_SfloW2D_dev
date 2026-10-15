#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-stages.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { in_sources=1; next }
    in_sources && /^[[:space:]]*$/ { exit }
    in_sources { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
baseline_revision=c648e9c5356c65d68d0624369d408cec781f0816
baseline_sha256=90d5d50923069798cd6135f18a940cdaa423a7d26644044ca9989e483d6bfd43
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
for profile in strict optimized; do
    mkdir "$build_dir/$profile"
    cd "$build_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    objects=''
    reference_objects=''
    for source in $sources; do
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            -I"$repo_dir/src/state" \
            $netcdf_fflags -c "$repo_dir/src/$source"
        object=${source##*/}
        objects="$objects ${object%.f90}.o"
        if [ "$source" != "spatial/reconstruction_2d.f90" ]; then
            reference_objects="$reference_objects $build_dir/$profile/${object%.f90}.o"
        fi
    done
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$test_dir/test_imex_stages.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_imex_stages.o -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages
    reference_dir="$build_dir/$profile/reference"
    mkdir "$reference_dir"
    git -C "$repo_dir" show \
        "$baseline_revision:src/spatial/reconstruction_2d.f90" \
        > "$reference_dir/reconstruction_2d.f90"
    actual_baseline_sha256=$(shasum -a 256 "$reference_dir/reconstruction_2d.f90" | awk '{print $1}')
    if [ "$actual_baseline_sha256" != "$baseline_sha256" ]; then
        echo "Pinned XY-5 reconstruction baseline checksum mismatch" >&2
        exit 1
    fi
    (
        cd "$reference_dir"
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            -I"$build_dir/$profile" -I"$repo_dir/src/state" -J"$reference_dir" \
            $netcdf_fflags -c reconstruction_2d.f90
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $reference_objects reconstruction_2d.o "$build_dir/$profile/test_imex_stages.o" \
            -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages_reference
    )
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$test_dir/test_stage_count.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_stage_count.o -llapack $netcdf_flibs $netcdf_clibs -o test_stage_count
    PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/compare_stages.py" \
        "$build_dir/$profile/test_imex_stages" "$profile" \
        --stage-count-executable "$build_dir/$profile/test_stage_count"
    PYTHONDONTWRITEBYTECODE=1 python3 "$test_dir/test_reconstruction_boundaries.py" \
        "$build_dir/$profile/test_imex_stages" "$profile" \
        --baseline-executable "$reference_dir/test_imex_stages_reference"
done
