# initial information about the project
AC_INIT([IMEX_SfloW2D], [1.0], [mattia.demichielivitturi@ingv.it])

# check if the source folder is available
AC_CONFIG_SRCDIR([src/IMEX_SfloW2D.f90])

# Preserve the historical optimized default unless a build profile or FCFLAGS
# was explicitly supplied.
AC_ARG_ENABLE([strict-debug],
  [AS_HELP_STRING([--enable-strict-debug], [build with runtime checks and floating-point traps])],
  [enable_strict_debug=$enableval],
  [enable_strict_debug=no]
)
if test -z "$FCFLAGS"; then
  if test "x$enable_strict_debug" = "xyes"; then
    FCFLAGS="-O0 -g -Wall -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow"
  else
    FCFLAGS="-Ofast -funroll-all-loops"
  fi
fi

# check for fortran compiler
AC_PROG_FC

# Compiler-specific OpenMP flags can be overridden for non-GNU toolchains.
AC_ARG_VAR([OPENMP_FCFLAGS], [Fortran OpenMP compile and link flags])
if test -z "$OPENMP_FCFLAGS"; then
  OPENMP_FCFLAGS="-fopenmp"
fi
AC_SUBST([OPENMP_FCFLAGS])

# automake initialization (mandatory) including a check for automake API version >= 1.9
AM_INIT_AUTOMAKE([1.9 foreign])

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_ARG_WITH([netcdf],
  AS_HELP_STRING([--with-netcdf=DIR], [location of NetCDF installation]),
  [NETCDF=$withval],
  [NETCDF="/usr/local/netcdf"]
)

AC_ARG_VAR([NETCDF_LIBS], [Fortran NetCDF libraries])
AC_ARG_VAR([LAPACK_LIBS], [LAPACK and BLAS libraries])
if test -z "$NETCDF_LIBS"; then
  NETCDF_LIBS="-lnetcdff"
fi
if test -z "$LAPACK_LIBS"; then
  LAPACK_LIBS="-llapack -lblas"
fi
if test "x$NETCDF" = "xno"; then
  NETCDF_CPPFLAGS=""
  NETCDF_LDFLAGS=""
else
  NETCDF_CPPFLAGS="-I$NETCDF/include"
  NETCDF_LDFLAGS="-L$NETCDF/lib"
fi

AC_SUBST([NETCDF])
AC_SUBST([NETCDF_CPPFLAGS])
AC_SUBST([NETCDF_LDFLAGS])
AC_SUBST([NETCDF_LIBS])
AC_SUBST([LAPACK_LIBS])

AC_PREFIX_DEFAULT([..])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
