
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
                pccu_2d.f90 \
                domain_2d.f90 \
                state_2d.f90 \
                mass_exchange_2d.f90 \
                nonlinear_solver_2d.f90 \
                hp_reconstruction_2d.f90 \
                reconstruction_2d.f90 \
                hyperbolic_2d.f90 \
                spatial_operator_2d.f90 \
                time_integration_2d.f90 \
                stochastic_random_2d.f90 \
                stochastic_module.f90 \
                simulation_2d.f90 \
                init_2d.f90 \
                inpout_2d.f90 \
                IMEX_SfloW2D.f90

EXTRA_DIST = state_conversion_phys_var.inc velocity_from_conservative.inc     \
             mixture_heat_capacity.inc mixture_from_mass.inc

state_conversion_2d.$(OBJEXT): state_conversion_phys_var.inc                  \
                               velocity_from_conservative.inc                 \
                               mixture_heat_capacity.inc                      \
                               mixture_from_mass.inc

lapack_interfaces_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
constitutive_parameters_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
runtime_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
complexify.$(OBJEXT): parameters_2d.$(OBJEXT)
geometry_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
state_conversion_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                         \
                               constitutive_parameters_2d.$(OBJEXT)
equation_terms_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                           \
                             constitutive_parameters_2d.$(OBJEXT)              \
                             diagnostics_2d.$(OBJEXT)                          \
                             equation_metadata_2d.$(OBJEXT)                    \
                             complexify.$(OBJEXT) geometry_2d.$(OBJEXT)        \
                             state_conversion_2d.$(OBJEXT)
nonconservative_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
pccu_2d.$(OBJEXT): parameters_2d.$(OBJEXT) nonconservative_2d.$(OBJEXT)
domain_2d.$(OBJEXT): parameters_2d.$(OBJEXT) geometry_2d.$(OBJEXT)
state_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                                   \
                    constitutive_parameters_2d.$(OBJEXT)                      \
                    geometry_2d.$(OBJEXT) model_layout_2d.$(OBJEXT)
mass_exchange_2d.$(OBJEXT): parameters_2d.$(OBJEXT) diagnostics_2d.$(OBJEXT)  \
                            domain_2d.$(OBJEXT) equation_terms_2d.$(OBJEXT)    \
                            geometry_2d.$(OBJEXT) state_conversion_2d.$(OBJEXT)
nonlinear_solver_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                        \
                               diagnostics_2d.$(OBJEXT)                       \
                               equation_metadata_2d.$(OBJEXT)                 \
                               equation_terms_2d.$(OBJEXT)                    \
                               lapack_interfaces_2d.$(OBJEXT)
hp_reconstruction_2d.$(OBJEXT): parameters_2d.$(OBJEXT) geometry_2d.$(OBJEXT)
reconstruction_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                          \
                              constitutive_parameters_2d.$(OBJEXT)             \
                              equation_terms_2d.$(OBJEXT) geometry_2d.$(OBJEXT)\
                              hp_reconstruction_2d.$(OBJEXT)                   \
                              state_conversion_2d.$(OBJEXT)
hyperbolic_2d.$(OBJEXT): parameters_2d.$(OBJEXT) equation_terms_2d.$(OBJEXT)  \
                          geometry_2d.$(OBJEXT) hp_reconstruction_2d.$(OBJEXT) \
                          nonconservative_2d.$(OBJEXT) pccu_2d.$(OBJEXT)       \
                          reconstruction_2d.$(OBJEXT)
spatial_operator_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                       \
                                constitutive_parameters_2d.$(OBJEXT)          \
                                domain_2d.$(OBJEXT) geometry_2d.$(OBJEXT)     \
                                hyperbolic_2d.$(OBJEXT)                       \
                                reconstruction_2d.$(OBJEXT)                  \
                                state_conversion_2d.$(OBJEXT)
time_integration_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                        \
                               constitutive_parameters_2d.$(OBJEXT)           \
                               diagnostics_2d.$(OBJEXT) domain_2d.$(OBJEXT)   \
                               equation_metadata_2d.$(OBJEXT)                 \
                               equation_terms_2d.$(OBJEXT) geometry_2d.$(OBJEXT)\
                               nonlinear_solver_2d.$(OBJEXT)                  \
                               spatial_operator_2d.$(OBJEXT)                  \
                               state_conversion_2d.$(OBJEXT)
stochastic_random_2d.$(OBJEXT): parameters_2d.$(OBJEXT)
stochastic_module.$(OBJEXT): parameters_2d.$(OBJEXT) domain_2d.$(OBJEXT)     \
                              geometry_2d.$(OBJEXT) state_2d.$(OBJEXT)        \
                              stochastic_random_2d.$(OBJEXT)
simulation_2d.$(OBJEXT): parameters_2d.$(OBJEXT) domain_2d.$(OBJEXT)          \
                          mass_exchange_2d.$(OBJEXT)                         \
                          equation_metadata_2d.$(OBJEXT)                      \
                          equation_terms_2d.$(OBJEXT)                         \
                          model_layout_2d.$(OBJEXT)                           \
                          nonlinear_solver_2d.$(OBJEXT)                       \
                          runtime_2d.$(OBJEXT)                                \
                          spatial_operator_2d.$(OBJEXT)                       \
                          state_2d.$(OBJEXT) stochastic_module.$(OBJEXT)      \
                          time_integration_2d.$(OBJEXT)
init_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                                   \
                    constitutive_parameters_2d.$(OBJEXT) geometry_2d.$(OBJEXT)\
                    state_2d.$(OBJEXT) state_conversion_2d.$(OBJEXT)
inpout_2d.$(OBJEXT): parameters_2d.$(OBJEXT)                                 \
                      time_integration_2d.$(OBJEXT)                           \
                      constitutive_parameters_2d.$(OBJEXT)                    \
                      diagnostics_2d.$(OBJEXT) domain_2d.$(OBJEXT)            \
                      geometry_2d.$(OBJEXT) init_2d.$(OBJEXT)                 \
                      runtime_2d.$(OBJEXT) state_2d.$(OBJEXT)                 \
                      state_conversion_2d.$(OBJEXT)                           \
                      stochastic_module.$(OBJEXT)                             \
                      stochastic_random_2d.$(OBJEXT)
IMEX_SfloW2D.$(OBJEXT): parameters_2d.$(OBJEXT)                              \
                         constitutive_parameters_2d.$(OBJEXT)                 \
                         geometry_2d.$(OBJEXT) init_2d.$(OBJEXT)              \
                         inpout_2d.$(OBJEXT) mass_exchange_2d.$(OBJEXT)       \
                         simulation_2d.$(OBJEXT) state_conversion_2d.$(OBJEXT)

AM_FCFLAGS      = -fopenmp -I@NETCDF@/include
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack -L@NETCDF@/lib -lnetcdff 

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran
NETCDF = @NETCDF@
