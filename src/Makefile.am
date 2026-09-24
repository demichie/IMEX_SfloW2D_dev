
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEX_SfloW2D
bindir=@bindir@

# list of sources for the binary
IMEX_SfloW2D_SOURCES = \
                parameters_2d.f90 \
                diagnostics_2d.f90 \
                lapack_interfaces_2d.f90 \
                constitutive_parameters_2d.f90 \
                model_layout_2d.f90 \
                equation_metadata_2d.f90 \
                runtime_2d.f90 \
                complexify.f90 \
                geometry_2d.f90 \
                state_conversion_2d.f90 \
                equation_terms_2d.f90 \
                nonconservative_2d.f90 \
                domain_2d.f90 \
                constitutive_2d.f90 \
                state_2d.f90 \
                mass_exchange_2d.f90 \
                nonlinear_solver_2d.f90 \
                reconstruction_2d.f90 \
                hyperbolic_2d.f90 \
                time_integration_2d.f90 \
                stochastic_random_2d.f90 \
                stochastic_module.f90 \
                simulation_2d.f90 \
                init_2d.f90 \
                inpout_2d.f90 \
                IMEX_SfloW2D.f90

EXTRA_DIST = state_conversion_phys_var.inc mixture_heat_capacity.inc          \
             mixture_from_mass.inc

state_conversion_2d.$(OBJEXT): state_conversion_phys_var.inc                  \
                               mixture_heat_capacity.inc                      \
                               mixture_from_mass.inc

FCFLAGS         = -fopenmp -Ofast -funroll-all-loops -I@NETCDF@/include -lnetcdff
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack -L@NETCDF@/lib -lnetcdff 

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran
NETCDF = @NETCDF@
