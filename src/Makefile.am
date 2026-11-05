
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
                numerics/lapack_interfaces_2d.f90 \
                core/constitutive_parameters_2d.f90 \
                core/model_layout_2d.f90 \
                numerics/source_catalog_2d.f90 \
                numerics/source_registry_2d.f90 \
                core/runtime_2d.f90 \
                core/complexify.f90 \
                geometry/geometry_2d.f90 \
                state/state_conversion_2d.f90 \
                state_new/state_new_types.f90 \
                state_new/state_new_composition.f90 \
                state_new/state_new_mixture.f90 \
                state_new/state_new_conversions.f90 \
                state_new/state_new_physical.f90 \
                state_new/state_new.f90 \
                physics/rheology_laws/rheology_voellmy_2d.f90 \
                physics/rheology_laws/rheology_plastic_2d.f90 \
                physics/rheology_laws/rheology_froude_2d.f90 \
                physics/rheology_laws/rheology_lucas_2d.f90 \
                physics/rheology_laws/rheology_lava_2d.f90 \
                physics/rheology_laws/rheology_lahar_2d.f90 \
                physics/rheology_laws/rheology_empirical_2d.f90 \
                physics/rheology_laws/rheology_quadratic_2d.f90 \
                physics/rheology_laws/rheology_newmodel_2d.f90 \
                physics/rheology_laws/rheology_mui_2d.f90 \
                physics/rheology_laws/rheology_mui_voellmy_2d.f90 \
                physics/rheology_laws/rheology_2d.f90 \
                physics/thermal_laws/thermal_cooling_2d.f90 \
                physics/thermal_laws/thermal_laws_2d.f90 \
                physics/pore_pressure/pore_pressure_darcy_2d.f90 \
                 physics/pore_pressure/pore_pressure_2d.f90 \
                physics/mass_exchange/mass_exchange_carrier_liquid_loss_2d.f90 \
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
                physics/mass_injection_terms/mass_injection_2d.f90 \
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
                numerics/lie_source_executor_2d.f90 \
                numerics/time_integration_2d.f90 \
                stochastic/stochastic_random_2d.f90 \
                stochastic/stochastic_module.f90 \
                core/simulation_2d.f90 \
                core/init_2d.f90 \
                io/output_paths.f90 \
                io/output_dirs.c \
                io/inpout_2d.f90 \
                IMEX_SfloW2D.f90

AM_FCFLAGS = @OPENMP_FCFLAGS@ @NETCDF_CPPFLAGS@ -I$(srcdir)/state -I$(srcdir)/state_new \
             -I$(srcdir)/physics/rheology_laws -I$(srcdir)/physics/thermal_laws
CLEANFILES = *.mod
EXTRA_DIST = state/state_conversion_phys_var.inc state/velocity_from_conservative.inc \
             state/mixture_heat_capacity.inc state/mixture_from_mass.inc \
             state_new/mixture_heat_capacity.inc state_new/mixture_from_mass.inc \
             state_new/density_from_mass.inc state_new/kinematics.inc \
             state_new/richardson.inc \
             physics/rheology_laws/rheology_voellmy_turbulent.inc \
             physics/rheology_laws/rheology_lava_viscous.inc \
             physics/rheology_laws/rheology_plastic_yield.inc \
             physics/rheology_laws/rheology_quadratic_drag.inc \
             physics/rheology_laws/rheology_empirical_drag.inc \
             physics/rheology_laws/rheology_lahar_drag.inc \
             physics/rheology_laws/rheology_newmodel_linear.inc \
             physics/thermal_laws/thermal_cooling_rate.inc \
             python_diagnostics/kt_ritter_1d.f90

state/state_conversion_2d.$(OBJEXT): state/state_conversion_phys_var.inc      \
                                     state/velocity_from_conservative.inc     \
                                     state/mixture_heat_capacity.inc          \
                                     state/mixture_from_mass.inc

# Keep runtime-sized, cell-local phase scratch on each OpenMP thread's stack.
# In this hot state-conversion unit, automatic arrays otherwise trigger frequent
# malloc/free calls for n_solid and n_add_gas scratch inside PCCU and speed bounds.
# Scope this GNU Fortran optimization to this object: global -fstack-arrays would
# risk moving unrelated large grid temporaries onto thread stacks.
state/state_conversion_2d.$(OBJEXT): AM_FCFLAGS += -fstack-arrays
spatial/hyperbolic_2d.$(OBJEXT): AM_FCFLAGS += -fstack-arrays
# Newton's small runtime-sized cell work arrays remain thread-private on stack.
numerics/nonlinear_solver_2d.$(OBJEXT): AM_FCFLAGS += -fstack-arrays
# Explicit source cell kernels hold only bounded per-cell scratch, never grids.
physics/source_terms_2d.$(OBJEXT): AM_FCFLAGS += -fstack-arrays

# F5: physical source consumers must wait for the canonical .mod producers,
# including under parallel make (Fortran modules are not inferred by automake).
physics/rheology_laws/rheology_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_plastic_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_froude_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_lucas_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_lava_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_lahar_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_empirical_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_newmodel_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_mui_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT) \
    state_new/state_new_mixture.$(OBJEXT)
physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT) \
    state_new/state_new_mixture.$(OBJEXT)
physics/thermal_laws/thermal_cooling_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/thermal_laws/thermal_laws_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/pore_pressure/pore_pressure_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/pore_pressure/pore_pressure_darcy_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_carrier_liquid_loss_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_deposition_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_exchange/mass_exchange_vaporization_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/mass_injection_terms/mass_injection_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
physics/rheology_laws/rheology_2d.$(OBJEXT): state_new/state_new_physical.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): state_new/state_new_physical.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
numerics/nonlinear_solver_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)
numerics/lie_source_executor_2d.$(OBJEXT): state_new/state_new_types.$(OBJEXT)

# F5: new general state kernels are built in module dependency order.
state_new/state_new_types.$(OBJEXT): core/parameters_2d.$(OBJEXT)
state_new/state_new_composition.$(OBJEXT): core/parameters_2d.$(OBJEXT)
# F5: physical particle descriptors belong to mixture physics.
state_new/state_new_mixture.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state_new/mixture_heat_capacity.inc \
                                     state_new/mixture_from_mass.inc \
                                     state_new/density_from_mass.inc
state_new/state_new_conversions.$(OBJEXT): state_new/state_new_types.$(OBJEXT) \
                                         state_new/state_new_mixture.$(OBJEXT) \
                                         state_new/state_new_composition.$(OBJEXT)
state_new/state_new_physical.$(OBJEXT): state_new/state_new_types.$(OBJEXT) \
                                      state_new/state_new_mixture.$(OBJEXT) \
                                      state_new/state_new_composition.$(OBJEXT) \
                                      state_new/kinematics.inc state_new/richardson.inc
state_new/state_new.$(OBJEXT): state_new/state_new_types.$(OBJEXT) \
                             state_new/state_new_composition.$(OBJEXT) \
                             state_new/state_new_mixture.$(OBJEXT) \
                             state_new/state_new_conversions.$(OBJEXT) \
                             state_new/state_new_physical.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT): state_new/state_new_mixture.$(OBJEXT) \
    state_new/state_new_composition.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT): state_new/state_new_mixture.$(OBJEXT) \
    state_new/state_new_composition.$(OBJEXT)
physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT): state_new/state_new_mixture.$(OBJEXT) \
    state_new/state_new_composition.$(OBJEXT)

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
                                        physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_plastic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)         \
                                        physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_froude_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)          \
                                       physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                       physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_lucas_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)            \
                                      physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
core/constitutive_parameters_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/runtime_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
core/complexify.$(OBJEXT): core/parameters_2d.$(OBJEXT)
geometry/geometry_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
state/state_conversion_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                    core/constitutive_parameters_2d.$(OBJEXT) \
                                    physics/rheology_laws/rheology_parameters_2d.$(OBJEXT)  \
                                    state_new/state_new_mixture.$(OBJEXT)        \
                                    physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)      \
                                            core/constitutive_parameters_2d.$(OBJEXT) \
                                            physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                            state_new/state_new_mixture.$(OBJEXT) \
                                            physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT)
physics/rheology_laws/rheology_lava_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                     physics/rheology_laws/rheology_parameters_2d.$(OBJEXT)  \
                                     physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_lahar_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)            \
                                      core/constitutive_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_parameters_2d.$(OBJEXT) \
                                      physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_empirical_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)       \
                                          physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_newmodel_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                                     core/diagnostics_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)        \
                                          physics/rheology_laws/rheology_analytic_utils_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/thermal_laws/thermal_cooling_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/thermal_laws/thermal_laws_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)               \
                                     physics/thermal_laws/thermal_cooling_2d.$(OBJEXT) \
                                      core/diagnostics_2d.$(OBJEXT)
physics/pore_pressure/pore_pressure_darcy_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)         \
                                    core/constitutive_parameters_2d.$(OBJEXT) \
                                    core/complexify.$(OBJEXT)                 \
                                    core/diagnostics_2d.$(OBJEXT)              \
                                    geometry/geometry_2d.$(OBJEXT)             \
                                    state/state_conversion_2d.$(OBJEXT)
physics/pore_pressure/pore_pressure_2d.$(OBJEXT): physics/pore_pressure/pore_pressure_darcy_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_carrier_liquid_loss_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_deposition_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_lahar_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_lava_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_vaporization_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_carrier_liquid_loss_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_deposition_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_vaporization_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT)
physics/mass_injection_terms/mass_boundary_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/mass_injection_terms/mass_injection_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                      physics/mass_injection_terms/mass_injection_2d.$(OBJEXT) \
                                      physics/mass_exchange/mass_exchange_2d.$(OBJEXT) \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT) \
                                     geometry/geometry_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT) \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_boundary_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_circular_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_bottom_fissural_injection_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_radial_injection_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_2d.$(OBJEXT) \
                                     core/diagnostics_2d.$(OBJEXT)
physics/equation_terms_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)             \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/thermal_laws/thermal_laws_2d.$(OBJEXT)         \
                                     physics/pore_pressure/pore_pressure_2d.$(OBJEXT)   \
                                     physics/source_terms_2d.$(OBJEXT) \
                                     physics/rheology_laws/rheology_2d.$(OBJEXT)           \
                                     core/diagnostics_2d.$(OBJEXT)             \
                                     numerics/source_catalog_2d.$(OBJEXT)      \
                                     numerics/source_registry_2d.$(OBJEXT)     \
                                     core/complexify.$(OBJEXT)                 \
                                     geometry/geometry_2d.$(OBJEXT)            \
                                     state_new/state_new_mixture.$(OBJEXT) \
                                     state_new/state_new_composition.$(OBJEXT)
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
                                        physics/equation_terms_2d.$(OBJEXT)     \
                                        numerics/lapack_interfaces_2d.$(OBJEXT) \
                                        numerics/source_registry_2d.$(OBJEXT)   \
                                        numerics/source_catalog_2d.$(OBJEXT)
spatial/reconstruction_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)              \
                                     core/constitutive_parameters_2d.$(OBJEXT) \
                                     physics/equation_terms_2d.$(OBJEXT) \
                                     physics/source_terms_2d.$(OBJEXT)        \
                                     geometry/geometry_2d.$(OBJEXT)             \
                                     state_new/state_new_composition.$(OBJEXT)
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
                                         numerics/source_registry_2d.$(OBJEXT)   \
                                         physics/equation_terms_2d.$(OBJEXT)     \
                                         geometry/geometry_2d.$(OBJEXT)          \
                                         numerics/nonlinear_solver_2d.$(OBJEXT)  \
                                         spatial/spatial_operator_2d.$(OBJEXT)  \
                                         state/state_conversion_2d.$(OBJEXT)  \
                                         state/state_2d.$(OBJEXT) \
                                         numerics/lie_source_executor_2d.$(OBJEXT)
numerics/lie_source_executor_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)       \
                                          core/constitutive_parameters_2d.$(OBJEXT) \
                                          core/diagnostics_2d.$(OBJEXT)       \
                                          geometry/domain_2d.$(OBJEXT)        \
                                          numerics/source_registry_2d.$(OBJEXT) \
                                          numerics/source_catalog_2d.$(OBJEXT) \
                                          physics/source_terms_2d.$(OBJEXT) \
                                          physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT) \
                                          numerics/nonlinear_solver_2d.$(OBJEXT) \
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
                              physics/equation_terms_2d.$(OBJEXT)               \
                              core/model_layout_2d.$(OBJEXT)                    \
                              numerics/nonlinear_solver_2d.$(OBJEXT)            \
                              core/runtime_2d.$(OBJEXT)                         \
                              spatial/spatial_operator_2d.$(OBJEXT)             \
                              state/state_2d.$(OBJEXT)                          \
                              stochastic/stochastic_module.$(OBJEXT)            \
                              numerics/time_integration_2d.$(OBJEXT)            \
                              numerics/lie_source_executor_2d.$(OBJEXT) \
                              physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT)
core/init_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)                           \
                        core/constitutive_parameters_2d.$(OBJEXT)              \
                        geometry/geometry_2d.$(OBJEXT)                         \
                        state/state_2d.$(OBJEXT)                                \
                        state/state_conversion_2d.$(OBJEXT)
io/output_paths.$(OBJEXT): core/parameters_2d.$(OBJEXT)                       \
                          core/diagnostics_2d.$(OBJEXT)
io/inpout_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
    state_new/state_new_mixture.$(OBJEXT)                           \
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

# Explicit Fortran module dependencies for deterministic parallel builds.
numerics/source_catalog_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
numerics/source_input_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT) \
                                     physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT)
numerics/source_registry_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT)
physics/source_configuration_2d.$(OBJEXT): numerics/source_registry_2d.$(OBJEXT) \
                                         numerics/source_input_2d.$(OBJEXT) \
                                         core/constitutive_parameters_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): physics/source_configuration_2d.$(OBJEXT) \
                                     physics/mass_injection_terms/mass_injection_coefficients_2d.$(OBJEXT) \
                              numerics/source_input_2d.$(OBJEXT)

# Additional module USE edges introduced by the Lie source architecture.
physics/rheology_laws/rheology_2d.$(OBJEXT): state/state_conversion_2d.$(OBJEXT)
physics/thermal_laws/thermal_laws_2d.$(OBJEXT): state/state_conversion_2d.$(OBJEXT)
physics/pore_pressure/pore_pressure_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT) \
                                     state/state_conversion_2d.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT) \
                                  numerics/source_registry_2d.$(OBJEXT)
numerics/nonlinear_solver_2d.$(OBJEXT): physics/source_terms_2d.$(OBJEXT)
physics/source_configuration_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT) \
                                          physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): numerics/source_registry_2d.$(OBJEXT)

# Remaining direct USE dependencies, audited against all production .f90 sources.
IMEX_SfloW2D.$(OBJEXT): numerics/source_registry_2d.$(OBJEXT) \
    numerics/source_catalog_2d.$(OBJEXT)
core/simulation_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
io/inpout_2d.$(OBJEXT): numerics/source_catalog_2d.$(OBJEXT)
numerics/source_input_2d.$(OBJEXT): core/constitutive_parameters_2d.$(OBJEXT) \
    core/parameters_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_carrier_liquid_loss_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_entrainment_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT) \
    core/constitutive_parameters_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_erosion_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/mass_exchange/mass_exchange_lava_deposition_2d.$(OBJEXT): state/state_conversion_2d.$(OBJEXT) \
    core/constitutive_parameters_2d.$(OBJEXT)
physics/pore_pressure/pore_pressure_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_empirical_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_froude_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_lahar_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_lava_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_lucas_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_plastic_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)
physics/source_configuration_2d.$(OBJEXT): core/parameters_2d.$(OBJEXT)
physics/source_terms_2d.$(OBJEXT): core/complexify.$(OBJEXT)
spatial/hyperbolic_2d.$(OBJEXT): state/state_conversion_2d.$(OBJEXT)
spatial/reconstruction_2d.$(OBJEXT): physics/mass_injection_terms/mass_boundary_injection_2d.$(OBJEXT)
state/state_conversion_2d.$(OBJEXT): core/diagnostics_2d.$(OBJEXT)

physics/rheology_laws/rheology_voellmy_2d.$(OBJEXT): physics/rheology_laws/rheology_voellmy_turbulent.inc

physics/rheology_laws/rheology_mui_voellmy_2d.$(OBJEXT): physics/rheology_laws/rheology_voellmy_turbulent.inc

physics/rheology_laws/rheology_lava_2d.$(OBJEXT): physics/rheology_laws/rheology_lava_viscous.inc

physics/rheology_laws/rheology_plastic_2d.$(OBJEXT): physics/rheology_laws/rheology_plastic_yield.inc

physics/rheology_laws/rheology_quadratic_2d.$(OBJEXT): physics/rheology_laws/rheology_quadratic_drag.inc

physics/rheology_laws/rheology_empirical_2d.$(OBJEXT): physics/rheology_laws/rheology_empirical_drag.inc

physics/rheology_laws/rheology_lahar_2d.$(OBJEXT): physics/rheology_laws/rheology_lahar_drag.inc
