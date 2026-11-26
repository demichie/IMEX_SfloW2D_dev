# Build from source

The source repository requires a Fortran compiler, Autoconf/Automake, NetCDF-Fortran, LAPACK/BLAS and Make. Use dependency libraries compatible with your chosen Fortran compiler.

On Ubuntu 24.04, a representative installation is:

```bash
sudo apt-get update
sudo apt-get install -y gfortran autoconf automake make \
    libnetcdff-dev libnetcdf-dev liblapack-dev libblas-dev
```

From the **repository root**:

```bash
autoreconf -fi
./configure --with-netcdf=/usr FC=gfortran
make -j2
make check
```

The executable is `src/IMEX_SfloW2D`. For a strict/debug build, configure with `--enable-strict-debug`. Generated Autotools artifacts are deliberately not committed.

## Run an example

```bash
cd EXAMPLES/EXAMPLE_DAMBREAK_TWO_LAYER_1D
python3 prepare.py
../../src/IMEX_SfloW2D
```

The example provides the input file and creates its DEM and initial thickness raster. Other examples may require extracting a pre-existing DEM ZIP; always follow their individual instructions.

## Other platforms

Use `--with-netcdf=/path/to/prefix` for a compiler-compatible NetCDF-Fortran installation. Consult the [repository README](https://github.com/demichie/IMEX_SfloW2D_dev#build-from-source) for build options.

A [Docker image](docker.md) avoids installing the Fortran compilation toolchain on the host once the image has been published.
