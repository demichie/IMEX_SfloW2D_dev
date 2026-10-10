#!/bin/sh
# Verify stationary cancellation and preservation of the previous moving operator.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
composition_dir="$repo_dir/TESTS/TEST_N7_COMPOSITION"
PYTHONDONTWRITEBYTECODE=1 python3 "$composition_dir/test_checker.py"
build_dir=$(mktemp -d /tmp/imex-hydrostatic-roundoff.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" = 1 ]; then
    echo "Keeping test work directory: $build_dir"
else
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi
# Freeze only the changed module. Both operators use the SAME other current
# objects and observer, isolating this correction; no .git is needed in an
# exported acceptance snapshot. The reference checksum is part of the contract.
python3 -c 'import hashlib,json,sys; c=json.load(open(sys.argv[1])); assert hashlib.sha256(open(sys.argv[2],"rb").read()).hexdigest()==c["baseline_module_sha256"]' \
    "$composition_dir/roundoff_guard_contract.json" "$test_dir/reference/hyperbolic_2d.f90"
netcdf_fflags=$(nf-config --fflags)
netcdf_flibs=$(nf-config --flibs)
netcdf_clibs=$(nc-config --libs)

build_observer() {
    source_dir=$1
    variant_dir=$2
    cd "$variant_dir"
    sources=$(awk '
        /^IMEX_SfloW2D_SOURCES =/ { reading=1; next }
        reading && /^[[:space:]]*$/ { exit }
        reading { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
    ' "$source_dir/Makefile.am")
    objects=''
    for source in $sources; do
        gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
            $netcdf_fflags -c "$source_dir/$source"
        objects="$objects ${source%.f90}.o"
    done
    # Use the same current observer for both production revisions.
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $netcdf_fflags -c "$repo_dir/TESTS/TEST_IMEX_STAGES/test_imex_stages.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $objects test_imex_stages.o -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages
}

for profile in strict optimized; do
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all' ;;
        optimized) flags='-Ofast -funroll-all-loops' ;;
    esac
    mkdir "$build_dir/$profile" "$build_dir/$profile/baseline" "$build_dir/$profile/candidate"
    build_observer "$repo_dir/src" "$build_dir/$profile/candidate"
    shared_objects=''
    for source in $sources; do
        if [ "$source" != hyperbolic_2d.f90 ]; then
            shared_objects="$shared_objects $build_dir/$profile/candidate/${source%.f90}.o"
        fi
    done
    cd "$build_dir/$profile/baseline"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        -I"$build_dir/$profile/candidate" $netcdf_fflags -c "$test_dir/reference/hyperbolic_2d.f90"
    gfortran -O0 -g -fcheck=all -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        -I"$build_dir/$profile/candidate" $netcdf_fflags -c "$repo_dir/TESTS/TEST_IMEX_STAGES/test_imex_stages.f90"
    gfortran $flags -fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow \
        $shared_objects hyperbolic_2d.o test_imex_stages.o -llapack $netcdf_flibs $netcdf_clibs -o test_imex_stages
    cd "$build_dir/$profile"
    PYTHONDONTWRITEBYTECODE=1 python3 "$composition_dir/check_roundoff_guard.py" \
        "$build_dir/$profile/baseline/test_imex_stages" \
        "$build_dir/$profile/candidate/test_imex_stages" "$profile"
    PYTHONDONTWRITEBYTECODE=1 python3 "$composition_dir/check_cases.py" \
        "$build_dir/$profile/candidate/test_imex_stages" "$profile" --equilibrium-only
done
