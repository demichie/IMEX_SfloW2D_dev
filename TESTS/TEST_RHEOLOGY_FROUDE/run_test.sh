#!/bin/sh
set -eu

test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
build_dir=$(mktemp -d /tmp/imex-rheology-froude.XXXXXX)
trap 'rm -rf "$build_dir"' EXIT HUP INT TERM

cd "$build_dir"
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { inside=1; next }
    inside && /^[[:space:]]*$/ { exit }
    inside { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)

# Trap invalid arithmetic in both profiles: optimized success alone must not
# conceal a zero/negative Froude denominator. No scientific tolerance is relaxed.
for profile in strict optimized; do
    case "$profile" in
        strict) profile_flags='-O0 -g -fcheck=all' ;;
        optimized) profile_flags='-Ofast -funroll-all-loops' ;;
    esac
    mkdir "$profile"
    cd "$profile"
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

    # Keep assertion/IEEE checks strict even when the numerical modules use
    # -Ofast, which otherwise permits optimizing away NaN/Inf checks.
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace \
        -ffpe-trap=invalid,zero,overflow -c "$test_dir/test_rheology_froude.f90"

    gfortran $profile_flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_rheology_froude.o -llapack $netcdf_flibs $netcdf_clibs -o test_rheology_froude
    ./test_rheology_froude
    cd ..
done
