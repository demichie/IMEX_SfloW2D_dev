
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
                physics/rheology_laws/rheology_analytic_utils_2d.f90 \
                physics/rheology_laws/rheology_parameters_2d.f90 \
                physics/rheology_laws/rheology_voellmy_2d.f90 \
                physics/rheology_laws/rheology_plastic_2d.f90 \
                physics/rheology_laws/rheology_froude_2d.f90 \
                physics/rheology_laws/rheology_lucas_2d.f90 \
                numerics/lapack_interfaces_2d.f90 \
                core/constitutive_parameters_2d.f90 \
                core/model_layout_2d.f90 \
                numerics/equation_metadata_2d.f90 \
                numerics/source_catalog_2d.f90 \
                numerics/source_registry_2d.f90 \
                core/runtime_2d.f90 \
                core/complexify.f90 \
                geometry/geometry_2d.f90 \
                state/state_conversion_2d.f90 \
                physics/rheology_laws/rheology_lava_2d.f90 \
                physics/rheology_laws/rheology_lahar_2d.f90 \
                physics/rheology_laws/rheology_empirical_2d.f90 \
                physics/rheology_laws/rheology_quadratic_2d.f90 \
                physics/rheology_laws/rheology_newmodel_2d.f90 \
                physics/rheology_laws/rheology_mui_2d.f90 \
                physics/rheology_laws/rheology_mui_voellmy_2d.f90 \
                physics/rheology_laws/rheology_2d.f90 \
                physics/thermal_laws_2d.f90 \
                physics/pore_pressure_laws_2d.f90 \
                physics/mass_exchange/mass_exchange_carrier_loss_2d.f90 \
                physics/mass_exchange/mass_exchange_erosion_2d.f90 \
                physics/mass_exchange/mass_exchange_deposition_2d.f90 \
                physics/mass_exchange/mass_exchange_lava_deposition_2d.f90 \
                numerics/source_input_2d.f90 \
                physics/mass_exchange/mass_exchange_vaporization_2d.f90 \
                physics/mass_exchange/mass_exchange_entrainment_2d.f90 \
                physics/mass_exchange/mass_exchange_2d.f90 \
                physics/source_configuration_2d.f90 \
                physics/mass_injection_terms/mass_injection_coefficients_2d.f90 \
                physics/mass_injection_terms/mass_boundary_injection_2d.f90 \
                physics/mass_injection_terms/mass_bottom_circular_injection_2d.f90 \
                physics/mass_injection_terms/mass_bottom_fissural_injection_2d.f90 \
                physics/mass_injection_terms/mass_radial_injection_2d.f90 \
                physics/source_terms_2d.f90 \
                physics/equation_terms_2d.f90 \
                spatial/nonconservative_2d.f90 \
                spatial/pccu_2d.f90 \
                geometry/domain_2d.f90 \
                state/state_2d.f90 \
                numerics/nonlinear_solver_2d.f90 \
                spatial/reconstruction_2d.f90 \
                spatial/hyperbolic_2d.f90 \
                spatial/spatial_operator_2d.f90 \
                numerics/time_integration_2d.f90 \
                numerics/step_source_executor_2d.f90 \
                stochastic/stochastic_random_2d.f90 \
                stochastic/stochastic_module.f90 \
                core/simulation_2d.f90 \
                core/init_2d.f90 \
                io/output_paths.f90 \
                io/output_dirs.c \
                io/inpout_2d.f90 \
                IMEX_SfloW2D.f90

AM_FCFLAGS = @OPENMP_FCFLAGS@ @NETCDF_CPPFLAGS@ -I$(srcdir)/state
CLEANFILES = *.mod
EXTRA_DIST = state/state_conversion_phys_var.inc state/velocity_from_conservative.inc \
             state/mixture_heat_capacity.inc state/mixture_from_mass.inc \
             python_diagnostics/kt_ritter_1d.f90

state/state_conversion_2d.$(OBJEXT): state/state_conversion_phys_var.inc      \
                                     state/velocity_from_conservative.inc     \
                                     state/mixture_heat_capacity.inc          \
                                     state/mixture_from_mass.inc

numerics/lapack_interfaces_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
physics/rheology_laws/rheology_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                 \
                                core/diagnostics_2d.$(OBJEXT)                 \
                                physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT)       \
                                physics/rheology_laws/rheology_plastic_2d.$(OBJEXT)       \
                                physics/rheology_laws/rheology_lava_2d.$(OBJEXT)          \
                                physics/rheology_laws/rheology_lahar_2d.$(OBJEXT)         \
                                physics/rheology_laws/rheology_empirical_2d.$(OBJEXT)     \
                                physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT)     \
                                physics/rheology_laws/rheology_newmodel_2d.$(OBJEXT)      \
                                physics/rheology_laws/rheology_froude_2d.$(OBJEXT)        \
                                physics/rheology_laws/rheology_lucas_2d.$(OBJEXT)         \
                                physics/rheology_laws/rheology_mui_2d.$(OBJEXT)           \
                                physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT)
physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)  \
                                                core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_parameters_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)         \
                                        physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                        physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_plastic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)         \
                                        physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_froude_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)          \
                                       physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                       physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_lucas_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)            \
                                      physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
core/constitutive_parameters_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/runtime_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/complexify.$(OBJEXT): core/parameters_2d.$(OBJEXT)
geometry/geometry_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
state/state_conversion_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                    core/constitutive_parameters_2d.$(OBJEXT) \
                                    physics/rheology_laws/rheology_parameters_2d.$(OBJEXT)  \
                                    state/state_conversion_2d.$(OBJEXT)        \
                                    physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)      \
                                            core/constitutive_parameters_2d.$(OBJEXT) \
                                            physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                            state/state_conversion_2d.$(OBJEXT) \
                                            physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_lava_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                     physics/rheology_laws/rheology_parameters_2d.$(OBJEXT)  \
                                     physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_lahar_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)            \
                                      core/constitutive_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_empirical_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)       \
                                          physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_newmodel_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                                     core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)        \
                                          physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/thermal_laws_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)               \
                                      core/diagnostics_2d.$(OBJEXT)
physics/pore_pressure_laws_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)         \
                                    core/constitutive_parameters_2d.$(OBJEXT) \
                                    core/complexify.$(OBJEXT)                 \
                                    core/diagnostics_2d.$(OBJEXT)              \
                                    geometry/geometry_2d.$(OBJEXT)             \
                                    state/state_conversion_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_carrier_loss_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                      core/constitutive_parameters_2d.$(OBJEXT) \
                                      physics/pore_pressure_laws_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_deposition_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_lahar_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_lava_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_vaporization_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/pore_pressure_laws_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_carrier_loss_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_deposition_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_vaporization_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT)
physics/mass_injection_terms/mass_boundary_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure_laws_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure_laws_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure_laws_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure_laws_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_boundary_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/equation_terms_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/thermal_laws_2d.$(OBJEXT)         \
                                     physics/pore_pressure_laws_2d.$(OBJEXT)   \
                                     physics/source_terms_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_2d.$(OBJEXT)           \
                                     core/diagnostics_2d.$(OBJEXT)             \
                                     numerics/equation_metadata_2d.$(OBJEXT)   \
                                     numerics/source_catalog_2d.$(OBJEXT)      \
                                     numerics/source_registry_2d.$(OBJEXT)     \
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
                          core/model_layout_2d.$(OBJEXT)                       \
                          state/state_conversion_2d.$(OBJEXT)                 \
                          core/diagnostics_2d.$(OBJEXT)
numerics/nonlinear_solver_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)           \
                                        core/diagnostics_2d.$(OBJEXT)           \
                                        numerics/equation_metadata_2d.$(OBJEXT) \
                                        physics/equation_terms_2d.$(OBJEXT)     \
                                        numerics/lapack_interfaces_2d.$(OBJEXT) \
                                        numerics/source_registry_2d.$(OBJEXT)   \
                                        numerics/source_catalog_2d.$(OBJEXT)
spatial/reconstruction_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/equation_terms_2d.$(OBJEXT) \
                                     physics/source_terms_2d.$(OBJEXT)        \
                                     geometry/geometry_2d.$(OBJEXT)             \
                                     state/state_conversion_2d.$(OBJEXT)
spatial/hyperbolic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                  \
                                 physics/equation_terms_2d.$(OBJEXT)            \
                                 geometry/geometry_2d.$(OBJEXT)                 \
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
                                         numerics/source_registry_2d.$(OBJEXT)   \
                                         physics/equation_terms_2d.$(OBJEXT)     \
                                         geometry/geometry_2d.$(OBJEXT)          \
                                         numerics/nonlinear_solver_2d.$(OBJEXT)  \
                                         spatial/spatial_operator_2d.$(OBJEXT)  \
                                         state/state_conversion_2d.$(OBJEXT)  \
                                         state/state_2d.$(OBJEXT)
numerics/step_source_executor_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)       \
                                          core/constitutive_parameters_2d.$(OBJEXT) \
                                          core/diagnostics_2d.$(OBJEXT)       \
                                          geometry/domain_2d.$(OBJEXT)        \
                                          numerics/source_registry_2d.$(OBJEXT) \
                                          numerics/source_catalog_2d.$(OBJEXT) \
                                          physics/source_terms_2d.$(OBJEXT) \
                                          physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT) \
                                          physics/mass_exchange/mass_exchange_2d.$(OBJEXT) \
                                          numerics/nonlinear_solver_2d.$(OBJEXT) \
                                          numerics/equation_metadata_2d.$(OBJEXT) \
                                          state/state_2d.$(OBJEXT) \
                                          geometry/geometry_2d.$(OBJEXT)     \
                                          state/state_conversion_2d.$(OBJEXT)
stochastic/stochastic_random_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
stochastic/stochastic_module.$(OBJEXT): core/parameters_2d.$(OBJEXT)          \
                                        geometry/domain_2d.$(OBJEXT)             \
                                        geometry/geometry_2d.$(OBJEXT)          \
                                        state/state_2d.$(OBJEXT)                 \
                                        stochastic/stochastic_random_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                     \
                              geometry/geometry_2d.$(OBJEXT)                    \
                              geometry/domain_2d.$(OBJEXT)                      \
                              numerics/equation_metadata_2d.$(OBJEXT)           \
                              physics/equation_terms_2d.$(OBJEXT)               \
                              core/model_layout_2d.$(OBJEXT)                    \
                              numerics/nonlinear_solver_2d.$(OBJEXT)            \
                              core/runtime_2d.$(OBJEXT)                         \
                              spatial/spatial_operator_2d.$(OBJEXT)             \
                              state/state_2d.$(OBJEXT)                          \
                              stochastic/stochastic_module.$(OBJEXT)            \
                              numerics/time_integration_2d.$(OBJEXT)            \
                              numerics/step_source_executor_2d.$(OBJEXT) \
                              physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT)
core/init_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                           \
                        core/constitutive_parameters_2d.$(OBJEXT)              \
                        geometry/geometry_2d.$(OBJEXT)                         \
                        state/state_2d.$(OBJEXT)                                \
                        state/state_conversion_2d.$(OBJEXT)
io/output_paths.$(OBJEXT): core/parameters_2d.$(OBJEXT)                       \
                          core/diagnostics_2d.$(OBJEXT)
io/inpout_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                           \
                        physics/rheology_laws/rheology_2d.$(OBJEXT)                        \
                        physics/rheology_laws/rheology_parameters_2d.$(OBJEXT)             \
                        physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT)                \
                        physics/rheology_laws/rheology_plastic_2d.$(OBJEXT)                \
                        physics/rheology_laws/rheology_lava_2d.$(OBJEXT)                   \
                        physics/rheology_laws/rheology_lahar_2d.$(OBJEXT)                  \
                        physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT)              \
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
                        stochastic/stochastic_random_2d.$(OBJEXT)              \
                        io/output_paths.$(OBJEXT)                              \
                        numerics/source_input_2d.$(OBJEXT)
IMEX_SfloW2D.$(OBJEXT): core/parameters_2d.$(OBJEXT)                         \
                         core/constitutive_parameters_2d.$(OBJEXT)             \
                         geometry/geometry_2d.$(OBJEXT)                        \
                         core/init_2d.$(OBJEXT) io/output_paths.$(OBJEXT)      \
                         io/inpout_2d.$(OBJEXT)                               \
                         numerics/time_integration_2d.$(OBJEXT)                \
                         core/simulation_2d.$(OBJEXT)                          \
                         state/state_conversion_2d.$(OBJEXT)

IMEX_SfloW2D_LDFLAGS = @OPENMP_FCFLAGS@ @NETCDF_LDFLAGS@
IMEX_SfloW2D_LDADD = @NETCDF_LIBS@ @LAPACK_LIBS@ -lm

# D17 metadata is built once and consumed by stage dispatch through the analytic
# bridge; no registry validation or copying is performed inside cell loops.
numerics/source_catalog_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
numerics/source_input_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT)
numerics/source_registry_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT) \
                                      numerics/equation_metadata_2d.$(OBJEXT)
physics/source_configuration_2d.$(OBJEXT): numerics/source_registry_2d.$(OBJEXT) \
                                         numerics/source_input_2d.$(OBJEXT) \
                                         core/constitutive_parameters_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): physics/source_configuration_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                              numerics/source_input_2d.$(OBJEXT)
