
AUTOMAKE_OPTIONS = foreign subdir-objects

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEX_SfloW2D
bindir=@bindir@

# list of sources for the binary
IMEX_SfloW2D_SOURCES = \
                core/parameters_2d.f90 \
                core/diagnostics_2d.f90 \
                numerics/lapack_interfaces_2d.f90 \
                core/constitutive_parameters_2d.f90 \
                core/model_layout_2d.f90 \
                numerics/equation_metadata_2d.f90 \
                core/runtime_2d.f90 \
                core/complexify.f90 \
                geometry/geometry_2d.f90 \
                state/state_conversion_2d.f90 \
                physics/equation_terms_2d.f90 \
                spatial/nonconservative_2d.f90 \
                spatial/pccu_2d.f90 \
                geometry/domain_2d.f90 \
                state/state_2d.f90 \
                physics/mass_exchange_2d.f90 \
                numerics/nonlinear_solver_2d.f90 \
                spatial/hp_reconstruction_2d.f90 \
                spatial/reconstruction_2d.f90 \
                spatial/hyperbolic_2d.f90 \
                spatial/spatial_operator_2d.f90 \
                numerics/time_integration_2d.f90 \
                stochastic/stochastic_random_2d.f90 \
                stochastic/stochastic_module.f90 \
                core/simulation_2d.f90 \
                core/init_2d.f90 \
                io/inpout_2d.f90 \
                IMEX_SfloW2D.f90

AM_FCFLAGS = @OPENMP_FCFLAGS@ @NETCDF_CPPFLAGS@ -I$(srcdir)/state
EXTRA_DIST = state/state_conversion_phys_var.inc state/velocity_from_conservative.inc \
             state/mixture_heat_capacity.inc state/mixture_from_mass.inc

state/state_conversion_2d.$(OBJEXT): state/state_conversion_phys_var.inc      \
                                     state/velocity_from_conservative.inc     \
                                     state/mixture_heat_capacity.inc          \
                                     state/mixture_from_mass.inc

numerics/lapack_interfaces_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/constitutive_parameters_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/runtime_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/complexify.$(OBJEXT): core/parameters_2d.$(OBJEXT)
geometry/geometry_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
state/state_conversion_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT)
physics/equation_terms_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)             \
                                     numerics/equation_metadata_2d.$(OBJEXT)   \
                                     core/complexify.$(OBJEXT)                 \
                                     geometry/geometry_2d.$(OBJEXT)            \
                                     state/state_conversion_2d.$(OBJEXT)
spatial/nonconservative_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
spatial/pccu_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                        \
                           spatial/nonconservative_2d.$(OBJEXT)
geometry/domain_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                     \
                              geometry/geometry_2d.$(OBJEXT)
state/state_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                         \
                          core/constitutive_parameters_2d.$(OBJEXT)            \
                          geometry/geometry_2d.$(OBJEXT)                       \
                          core/model_layout_2d.$(OBJEXT)
physics/mass_exchange_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)               \
                                    core/diagnostics_2d.$(OBJEXT)               \
                                    geometry/domain_2d.$(OBJEXT)                \
                                    physics/equation_terms_2d.$(OBJEXT)         \
                                    geometry/geometry_2d.$(OBJEXT)              \
                                    state/state_conversion_2d.$(OBJEXT)
numerics/nonlinear_solver_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)           \
                                        core/diagnostics_2d.$(OBJEXT)           \
                                        numerics/equation_metadata_2d.$(OBJEXT) \
                                        physics/equation_terms_2d.$(OBJEXT)     \
                                        numerics/lapack_interfaces_2d.$(OBJEXT)
spatial/hp_reconstruction_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)           \
                                        geometry/geometry_2d.$(OBJEXT)
spatial/reconstruction_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/equation_terms_2d.$(OBJEXT)        \
                                     geometry/geometry_2d.$(OBJEXT)             \
                                     spatial/hp_reconstruction_2d.$(OBJEXT)     \
                                     state/state_conversion_2d.$(OBJEXT)
spatial/hyperbolic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                  \
                                 physics/equation_terms_2d.$(OBJEXT)            \
                                 geometry/geometry_2d.$(OBJEXT)                 \
                                 spatial/hp_reconstruction_2d.$(OBJEXT)         \
                                 spatial/nonconservative_2d.$(OBJEXT)           \
                                 spatial/pccu_2d.$(OBJEXT)                      \
                                 spatial/reconstruction_2d.$(OBJEXT)
spatial/spatial_operator_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)            \
                                       core/diagnostics_2d.$(OBJEXT)            \
                                       core/constitutive_parameters_2d.$(OBJEXT) \
                                       geometry/domain_2d.$(OBJEXT)              \
                                       geometry/geometry_2d.$(OBJEXT)           \
                                       spatial/hyperbolic_2d.$(OBJEXT)          \
                                       spatial/reconstruction_2d.$(OBJEXT)      \
                                       state/state_conversion_2d.$(OBJEXT)
numerics/time_integration_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)          \
                                         core/constitutive_parameters_2d.$(OBJEXT) \
                                         core/diagnostics_2d.$(OBJEXT)          \
                                         geometry/domain_2d.$(OBJEXT)            \
                                         numerics/equation_metadata_2d.$(OBJEXT) \
                                         physics/equation_terms_2d.$(OBJEXT)     \
                                         geometry/geometry_2d.$(OBJEXT)          \
                                         numerics/nonlinear_solver_2d.$(OBJEXT)  \
                                         spatial/spatial_operator_2d.$(OBJEXT)  \
                                         state/state_conversion_2d.$(OBJEXT)
stochastic/stochastic_random_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
stochastic/stochastic_module.$(OBJEXT): core/parameters_2d.$(OBJEXT)          \
                                        geometry/domain_2d.$(OBJEXT)             \
                                        geometry/geometry_2d.$(OBJEXT)          \
                                        state/state_2d.$(OBJEXT)                 \
                                        stochastic/stochastic_random_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                     \
                              geometry/domain_2d.$(OBJEXT)                      \
                              physics/mass_exchange_2d.$(OBJEXT)                \
                              numerics/equation_metadata_2d.$(OBJEXT)           \
                              physics/equation_terms_2d.$(OBJEXT)               \
                              core/model_layout_2d.$(OBJEXT)                    \
                              numerics/nonlinear_solver_2d.$(OBJEXT)            \
                              core/runtime_2d.$(OBJEXT)                         \
                              spatial/spatial_operator_2d.$(OBJEXT)             \
                              state/state_2d.$(OBJEXT)                          \
                              stochastic/stochastic_module.$(OBJEXT)            \
                              numerics/time_integration_2d.$(OBJEXT)
core/init_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                           \
                        core/constitutive_parameters_2d.$(OBJEXT)              \
                        geometry/geometry_2d.$(OBJEXT)                         \
                        state/state_2d.$(OBJEXT)                               \
                        state/state_conversion_2d.$(OBJEXT)
io/inpout_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                           \
                        numerics/time_integration_2d.$(OBJEXT)                 \
                        core/constitutive_parameters_2d.$(OBJEXT)              \
                        core/diagnostics_2d.$(OBJEXT)                           \
                        geometry/domain_2d.$(OBJEXT)                           \
                        geometry/geometry_2d.$(OBJEXT)                         \
                        core/init_2d.$(OBJEXT)                                 \
                        core/runtime_2d.$(OBJEXT)                              \
                        state/state_2d.$(OBJEXT)                               \
                        state/state_conversion_2d.$(OBJEXT)                    \
                        stochastic/stochastic_module.$(OBJEXT)                 \
                        stochastic/stochastic_random_2d.$(OBJEXT)
IMEX_SfloW2D.$(OBJEXT): core/parameters_2d.$(OBJEXT)                         \
                         core/constitutive_parameters_2d.$(OBJEXT)             \
                         geometry/geometry_2d.$(OBJEXT)                        \
                         core/init_2d.$(OBJEXT) io/inpout_2d.$(OBJEXT)         \
                         physics/mass_exchange_2d.$(OBJEXT)                     \
                         core/simulation_2d.$(OBJEXT)                          \
                         state/state_conversion_2d.$(OBJEXT)

IMEX_SfloW2D_LDFLAGS = @OPENMP_FCFLAGS@ @NETCDF_LDFLAGS@
IMEX_SfloW2D_LDADD = @NETCDF_LIBS@ @LAPACK_LIBS@ -lm
