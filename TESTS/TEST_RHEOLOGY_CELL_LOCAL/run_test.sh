#!/bin/sh
# D02 characterization and D05/D06 extraction regression; nonzero means a check failed.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
work_dir=${D02_WORK_DIR:-$(mktemp -d /tmp/imex-d02.XXXXXX)}
python=${PYTHON:-python3}
mkdir -p "$work_dir"
work_dir=$(CDPATH= cd -- "$work_dir" && pwd)
echo "D02 evidence directory: $work_dir"
compiler=${FC:-gfortran}
failed=0
for profile in strict optimized; do
    mkdir -p "$work_dir/$profile"
    cd "$work_dir/$profile"
    case "$profile" in
        strict) flags='-O0 -g -fcheck=all -finit-real=snan' ;;
        optimized) flags='-Ofast -g -funroll-all-loops' ;;
    esac
    common='-fopenmp -fbacktrace -ffpe-trap=invalid,zero,overflow -ffree-line-length-none'
    objects=''
    for source in core/parameters_2d core/diagnostics_2d rheology/rheology_2d \
        rheology/rheology_analytic_utils_2d rheology/rheology_parameters_2d \
        rheology/rheology_voellmy_2d rheology/rheology_plastic_2d \
        rheology/rheology_froude_2d rheology/rheology_lucas_2d \
        numerics/lapack_interfaces_2d \
        core/constitutive_parameters_2d numerics/equation_metadata_2d core/complexify \
        geometry/geometry_2d state/state_conversion_2d \
        rheology/rheology_lava_2d rheology/rheology_lahar_2d \
        rheology/rheology_empirical_2d rheology/rheology_quadratic_2d \
        rheology/rheology_mui_2d rheology/rheology_mui_voellmy_2d \
        physics/thermal_sources_2d physics/pore_pressure_2d \
        physics/equation_terms_2d \
        numerics/nonlinear_solver_2d; do
        "$compiler" $flags $common -I"$repo_dir/src/state" -c "$repo_dir/src/$source.f90"
        objects="$objects ${source##*/}.o"
    done
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_model_names.f90"
    "$compiler" $common parameters_2d.o rheology_2d.o test_model_names.o \
        -o test_model_names
    ./test_model_names
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_source_modules.f90"
    "$compiler" $common $objects test_source_modules.o -llapack \
        -o test_source_modules
    ./test_source_modules
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_analytic_velocity.f90"
    "$compiler" $common $objects test_analytic_velocity.o -llapack \
        -o test_analytic_velocity
    ./test_analytic_velocity
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_lava_cooling.f90"
    "$compiler" $common parameters_2d.o diagnostics_2d.o thermal_sources_2d.o \
        test_lava_cooling.o -o test_lava_cooling
    ./test_lava_cooling
    "$compiler" -O0 -g -fcheck=all $common -c \
        "$test_dir/test_pore_pressure_diffusion.f90"
    "$compiler" $common $objects test_pore_pressure_diffusion.o -llapack \
        -o test_pore_pressure_diffusion
    ./test_pore_pressure_diffusion
    "$compiler" -O0 -g -fcheck=all $common -I"$repo_dir/src/state" -c \
        "$test_dir/test_grouped_implicit_stage.f90"
    "$compiler" $common $objects test_grouped_implicit_stage.o -llapack \
        -o test_grouped_implicit_stage
    ./test_grouped_implicit_stage
    # IEEE checks and the test-side FD subtraction stay strict in both profiles.
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/probe_sources.f90"
    "$compiler" $common $objects probe_sources.o -llapack -o probe_sources
    if [ "$profile" = optimized ]; then
        # Distinguish an exception exposed by traps from an ordinary optimized
        # run. Production module objects are identical in both executables.
        "$compiler" -O0 -g -fcheck=all -fopenmp -ffree-line-length-none \
            -c "$test_dir/probe_sources.f90" -o probe_untrapped.o
        "$compiler" -fopenmp $objects probe_untrapped.o -llapack -o probe_untrapped
    fi
    PYTHONDONTWRITEBYTECODE=1 "$python" "$test_dir/characterize.py" \
        "$work_dir/$profile/probe_sources" "$work_dir/$profile" "$profile" || failed=1
done
exit "$failed"
