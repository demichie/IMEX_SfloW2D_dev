# Building the solver

IMEX_SfloW2D uses a Fortran compiler, Autotools, NetCDF-Fortran and LAPACK/BLAS. See the [repository README](https://github.com/demichie/IMEX_SfloW2D_dev#readme) for the authoritative platform-specific instructions.

A typical Autotools workflow is:

```bash
autoreconf -fi
./configure --with-netcdf=/path/to/netcdf-fortran
make -j2
make check
```

Use the correct NetCDF-Fortran installation for your compiler. Once built, consult an individual [example](../examples/index.md) for the necessary preparation, inputs and execution command.
