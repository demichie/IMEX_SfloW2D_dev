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
    sources=$(awk '
        /^IMEX_SfloW2D_SOURCES =/ { found=1; next }
        found && /^[[:space:]]*$/ { exit }
        found { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
    ' "$repo_dir/src/Makefile.am")
    netcdf_fflags=$(nf-config --fflags)
    netcdf_flibs=$(nf-config --flibs)
    netcdf_clibs=$(nc-config --libs)
    for source in $sources; do
        object=${source##*/}
        case "$source" in
            *.c) ${CC:-cc} -c "$repo_dir/src/$source" -o "${object%.c}.o" ;;
            *.f90) "$compiler" $flags $common -I"$repo_dir/src/state" -I"$repo_dir/src/physics/rheology_laws" -I"$repo_dir/src/physics/thermal_laws" \
                $netcdf_fflags -c "$repo_dir/src/$source" ;;
        esac
        objects="$objects ${object%.*}.o"
    done
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_model_names.f90"
    model_objects="$objects"
    "$compiler" $common $model_objects test_model_names.o -llapack $netcdf_flibs $netcdf_clibs -o test_model_names
    ./test_model_names
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_rheology_namelist_alias.f90"
    "$compiler" $common $objects test_rheology_namelist_alias.o -llapack $netcdf_flibs $netcdf_clibs \
        -o test_rheology_namelist_alias
    ./test_rheology_namelist_alias
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_source_modules.f90"
    "$compiler" $common $objects test_source_modules.o -llapack $netcdf_flibs $netcdf_clibs \
        -o test_source_modules
    ./test_source_modules
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_analytic_velocity.f90"
    "$compiler" $common $objects test_analytic_velocity.o -llapack $netcdf_flibs $netcdf_clibs \
        -o test_analytic_velocity
    ./test_analytic_velocity
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_shared_implicit_kernels.f90"
    "$compiler" $common $objects test_shared_implicit_kernels.o -llapack $netcdf_flibs $netcdf_clibs \
        -o test_shared_implicit_kernels
    ./test_shared_implicit_kernels
    "$compiler" -O0 -g -fcheck=all $common -c "$test_dir/test_lava_cooling.f90"
    "$compiler" $common parameters_2d.o diagnostics_2d.o thermal_cooling_2d.o \
        thermal_laws_2d.o test_lava_cooling.o -o test_lava_cooling
    ./test_lava_cooling
    "$compiler" -O0 -g -fcheck=all $common -c \
        "$test_dir/test_pore_pressure_diffusion.f90"
    "$compiler" $common $objects test_pore_pressure_diffusion.o -llapack $netcdf_flibs $netcdf_clibs \
        -o test_pore_pressure_diffusion
    ./test_pore_pressure_diffusion
    # The grouped IMEX-stage and probe characterization harnesses were
    # retired; Lie group-local Newton coverage is tracked in D19.
done
exit "$failed"
