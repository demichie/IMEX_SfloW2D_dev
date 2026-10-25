#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-curvature-source.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi

cd "$build_dir"

sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { inside=1; next }
    inside && /^[[:space:]]*$/ { exit }
    inside { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
    objects=''
    for source in $sources; do
        object=${source##*/}
        case "$source" in
            *.c) ${CC:-cc} -c "$repo_dir/src/$source" -o "${object%.c}.o" ;;
            *) gfortran ${profile_flags:--O0 -g -fcheck=all} -fopenmp -fbacktrace \
                 -ffpe-trap=invalid,zero,overflow -I"$repo_dir/src/state" \
                 $netcdf_fflags -c "$repo_dir/src/$source" ;;
        esac
        objects="$objects ${object%.*}.o"
    done
gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
    $netcdf_fflags -c "$test_dir/test_curvature_source.f90"
gfortran -O0 -g -fopenmp $objects test_curvature_source.o -llapack \
    $netcdf_flibs $netcdf_clibs -o test_curvature_source
./test_curvature_source
