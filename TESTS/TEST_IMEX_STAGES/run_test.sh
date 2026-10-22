#!/bin/sh
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
python=${PYTHON:-python3}
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
profiles=${IMEX_TEST_PROFILES:-"strict optimized"}
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)
for profile in $profiles; do
    mkdir "$build_dir/$profile"
    cd "$build_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    objects=''
    for source in $sources; do
        object=${source##*/}
        case "$source" in
            *.c)
                cc -c "$repo_dir/src/$source" -o "$object.o"
                ;;
            *)
                gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
                    -I"$repo_dir/src/state" \
                    $netcdf_fflags -c "$repo_dir/src/$source"
                ;;
        esac
        objects="$objects ${object%.f90}.o"
    done
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$test_dir/test_imex_stages.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_imex_stages.o -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages
    section=${IMEX_TEST_SECTION:-}
    if [ -z "$section" ]; then
        gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $netcdf_fflags -c "$test_dir/test_stage_count.f90"
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $objects test_stage_count.o -llapack $netcdf_flibs $netcdf_clibs -o test_stage_count
        PYTHONDONTWRITEBYTECODE=1 "$python" "$test_dir/compare_stages.py" \
            "$build_dir/$profile/test_imex_stages" "$profile" \
            --stage-count-executable "$build_dir/$profile/test_stage_count"
    else
        PYTHONDONTWRITEBYTECODE=1 "$python" "$test_dir/compare_stages.py" \
            "$build_dir/$profile/test_imex_stages" "$profile" --section "$section"
    fi
    if [ -z "$section" ] || [ "$section" = stages ]; then
        # The boundary driver asserts analytic trace values; historical spatial
        # modules cannot be linked to the current driver's changed workspace types.
        PYTHONDONTWRITEBYTECODE=1 "$python" "$test_dir/test_reconstruction_boundaries.py" \
            "$build_dir/$profile/test_imex_stages" "$profile"
    fi
done
