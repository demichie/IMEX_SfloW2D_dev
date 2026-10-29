#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-test-state-conversion.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
cd "$build_dir"

# Compile exactly the same current production modules and dependency order
# as the N8 mass-exchange gate. No retired source paths or test-only stubs.
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { inside=1; next }
    inside && /^[[:space:]]*$/ { exit }
    inside { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
for source in $sources; do
    case "$source" in
        *.c)
            cc -c "$repo_dir/src/$source" -o "${source##*/}.o"
            ;;
        *.f90)
            gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace \
                -I"$repo_dir/src/state" $netcdf_fflags \
                -c "$repo_dir/src/$source"
            ;;
    esac
done
gfortran -O0 -g -Wall -Wextra -fopenmp -fcheck=all -fbacktrace \
    -I"$repo_dir/src/state" $netcdf_fflags \
    -c "$test_dir/test_state_conversion.f90"
gfortran -O0 -g -fopenmp -fcheck=all -fbacktrace ./*.o \
    -llapack -lblas $netcdf_flibs $netcdf_clibs -lm -o test_state_conversion
./test_state_conversion
