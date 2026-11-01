#!/bin/sh
# Lie-compatible mass-exchange physics and geometric projection regressions.
set -eu
test_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$test_dir/../.." && pwd)
python3 -B "$test_dir/test_projection_checker.py"
python3 -B "$test_dir/check_coupled_limiters.py"
python3 -B "$test_dir/check_process_momentum.py"
build_dir=$(mktemp -d /tmp/lie-projection.XXXXXX)
if [ "${KEEP_TEST_WORKDIR:-0}" != 1 ]; then
    trap 'rm -rf "$build_dir"' EXIT HUP INT TERM
fi
sources=$(awk '
    /^IMEX_SfloW2D_SOURCES =/ { found=1; next }
    found && /^[[:space:]]*$/ { exit }
    found { gsub(/\\/, ""); if ($1 != "IMEX_SfloW2D.f90") print $1 }
' "$repo_dir/src/Makefile.am")
fflags=$(nf-config --fflags)
flibs=$(nf-config --flibs)
clibs=$(nc-config --libs)
for profile in strict optimized; do
    mkdir -p "$build_dir/$profile"
    cd "$build_dir/$profile"
    if [ "$profile" = strict ]; then flags='-O0 -g -fcheck=all'; else flags='-Ofast -funroll-all-loops'; fi
    for source in $sources; do
        case "$source" in
            *.c) ${CC:-cc} -c "$repo_dir/src/$source" ;;
            *.f90) ${FC:-gfortran} $flags -fopenmp -fbacktrace -ffree-line-length-none \
               -I"$repo_dir/src/state" $fflags -c "$repo_dir/src/$source" ;;
        esac
    done
    ${FC:-gfortran} -O0 -g -fcheck=all -fopenmp -ffree-line-length-none \
      $fflags -c "$test_dir/test_projection_stress.f90"
    ${FC:-gfortran} -fopenmp ./*.o -llapack $flibs $clibs -o test_projection_stress
    for threads in 1 4; do
        mkdir "threads-$threads"
        (cd "threads-$threads"; OMP_NUM_THREADS=$threads OMP_DYNAMIC=FALSE \
            ../test_projection_stress "$threads" > projection.log)
        python3 "$test_dir/check_projection_stress.py" \
            "threads-$threads/projection.log" "$profile" "$threads" > "threads-$threads/evidence.log"
    done
    cmp threads-1/projection_snapshot.bin threads-4/projection_snapshot.bin
done
