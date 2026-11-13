!********************************************************************************
!> \brief Simulation component aggregation and lifecycle
!
!> This module owns the state, runtime data and numerical workspaces that form
!> one simulation and coordinates their lifecycle.
!>
!> Owns state, runtime, domain, layout, source registry and numerical workspaces. Initialization
!> follows their dependency order; finalization releases their storage.
!********************************************************************************

MODULE simulation_2d

    USE mass_injection_coefficients_2d, ONLY : set_source_time_coefficients

    USE source_registry_2d, ONLY : source_registry_type
    USE source_configuration_2d, ONLY : initialize_configured_sources
    USE diagnostics_2d, ONLY : fatal_error

    USE rheology_parameters_2d, ONLY : rheology_l1
    USE thermal_parameters_2d, ONLY : thermal_l1
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_l1, pore_pressure_cache_l1
    USE model_layout_2d, ONLY : model_layout_type
    USE layer_physics_2d, ONLY : layer_physics_type, LAYER_PHASE_GAS, &
        LAYER_PHASE_LIQUID, LAYER_PHASE_LAVA

    USE parameters_2d, ONLY : wp, n_layers, n_vars, n_eqns, n_solid, n_add_gas,                                        &
            n_stoch_vars, n_pore_vars, liquid_flag, gas_flag, lava_flag, &
            rheology_flag, stochastic_flag,                                         &
            configure_canonical_state_indices


    USE runtime_2d, ONLY : runtime_state_type
    USE state_fields, ONLY : state_type
    USE domain_2d, ONLY : domain_type
    USE spatial_operator_2d, ONLY : spatial_operator_type
    USE time_integration_2d, ONLY : time_integration_workspace_type
    USE stochastic_module, ONLY : stochastic_workspace_type

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: simulation_context_type
    PUBLIC :: simulation
    PUBLIC :: TIMESTEP_EVENT_SOURCE_WINDOWS, TIMESTEP_EVENT_STOCHASTIC_PREPARED
    PUBLIC :: TIMESTEP_EVENT_LIE_COMPLETE
    PUBLIC :: TIMESTEP_EVENT_CLOCK_ADVANCED

    INTEGER, PARAMETER :: TIMESTEP_EVENT_SOURCE_WINDOWS=1
    INTEGER, PARAMETER :: TIMESTEP_EVENT_STOCHASTIC_PREPARED=2
    INTEGER, PARAMETER :: TIMESTEP_EVENT_LIE_COMPLETE=3
    INTEGER, PARAMETER :: TIMESTEP_EVENT_CLOCK_ADVANCED=5

    ABSTRACT INTERFACE
        !> \brief Report one completed production-timestep phase to an observer.
        !> \param[in] event Completed phase identifier from TIMESTEP_EVENT_*.
        !> \param[in] time Step-start time, or the advanced time for the clock event [s].
        !> \param[in] dt Timestep length used by the completed phase [s].
    SUBROUTINE timestep_event_observer(event,time,dt)
            IMPORT :: wp
            INTEGER, INTENT(IN) :: event
            REAL(wp), INTENT(IN) :: time,dt
    END SUBROUTINE timestep_event_observer
    END INTERFACE

    !> \brief Ownership root for runtime, states and reusable numerical workspaces.
    !> \details Initialization follows the model-layout dependencies; finalization
    !>          releases the owned storage without moving it into global facades.
    TYPE :: simulation_context_type

        TYPE(runtime_state_type) :: runtime
        TYPE(model_layout_type) :: model_layout
        TYPE(layer_physics_type), ALLOCATABLE :: layer_physics(:) !< Validated before worker allocation.
        TYPE(source_registry_type) :: source_registry !< Ordered method-capable Lie source groups.
        TYPE(state_type) :: state
        TYPE(domain_type) :: domain
        TYPE(spatial_operator_type) :: spatial_operator
        TYPE(time_integration_workspace_type) :: time_integration
        TYPE(stochastic_workspace_type) :: stochastic

    CONTAINS

        PROCEDURE :: initialize => initialize_simulation
        PROCEDURE :: finalize => finalize_simulation
        PROCEDURE :: advance_timestep => advance_simulation_timestep

    END TYPE simulation_context_type

    TYPE(simulation_context_type) :: simulation

CONTAINS

    !> \brief Initialize model metadata and simulation-owned numerical workspaces.
    !>
    !> \param[in,out] this Simulation-owned state, configuration descriptors and numerical
    !> workspaces.

    SUBROUTINE initialize_simulation(this)

        CLASS(simulation_context_type), INTENT(INOUT) :: this

        CHARACTER(LEN=256) :: source_error
        INTEGER :: carrier_phase

        ! Canonical state indices must be available before any initial-condition
        ! conversion, including collapsing-volume initialization.
        CALL configure_canonical_state_indices()

        CALL this%model_layout%initialize( n_layers, n_vars, n_eqns, n_solid,                                          &
            n_add_gas, n_stoch_vars, n_pore_vars )

        ! Validate phase/rheology and bind canonical thermal configuration before
        ! source, state and OpenMP workspace initialization. The thermal pointer
        ! is non-owning; this descriptor never caches another coefficient copy.
        IF (ALLOCATED(this%layer_physics)) DEALLOCATE(this%layer_physics)
        ALLOCATE(this%layer_physics(this%model_layout%n_layers))
        carrier_phase=LAYER_PHASE_GAS
        IF (liquid_flag) carrier_phase=LAYER_PHASE_LIQUID
        IF (lava_flag) carrier_phase=LAYER_PHASE_LAVA
        CALL this%layer_physics(1)%initialize(carrier_phase,rheology_flag,rheology_l1,thermal_l1, &
            pore_pressure_l1,source_error)
        IF (LEN_TRIM(source_error)>0) CALL fatal_error(TRIM(source_error))

        CALL initialize_configured_sources(this%source_registry,source_error)
        IF (LEN_TRIM(source_error)>0) CALL fatal_error(TRIM(source_error))

        CALL this%state%initialize( this%model_layout )

        CALL this%spatial_operator%initialize
        CALL this%domain%initialize

        CALL this%time_integration%initialize(this%source_registry)

        CALL this%stochastic%initialize

        WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'

    END SUBROUTINE initialize_simulation

    !> \brief Apply one production timestep after the driver has selected dt.
    !> \details Preserves the fractional-step ordering: source-window averages,
    !>          OU preparation, full Lie spatial/source advancement, then one clock advance.
    !> \param[in,out] this Initialized simulation context.
    !> \param[in] observer Optional test hook called after each completed phase.
    SUBROUTINE advance_simulation_timestep(this,observer,pre_spatial_primitives_ready)
        CLASS(simulation_context_type), INTENT(INOUT) :: this
        PROCEDURE(timestep_event_observer), OPTIONAL :: observer
        LOGICAL, INTENT(IN), OPTIONAL :: pre_spatial_primitives_ready
        REAL(wp) :: step_start,step_length

        step_start=this%runtime%t
        step_length=this%runtime%dt

        CALL set_source_time_coefficients(step_start,step_length)
        IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_SOURCE_WINDOWS,step_start,step_length)

        ! The stochastic OU stage belongs between spatial transport and
        ! Lie sources, inside the integrator's unpublished trial transaction.
        IF (PRESENT(observer)) THEN
            CALL this%time_integration%advance(this%state,step_start,step_length,                                      &
                this%stochastic,this%domain,this%spatial_operator,                                              &
                pre_spatial_primitives_ready=pre_spatial_primitives_ready,                                      &
                stochastic_prepared_hook=notify_stochastic_prepared,                                                   &
                thermal=this%layer_physics(1)%thermal_configuration(), &
                rheology=this%layer_physics(1)%rheology_configuration(), &
                pore_pressure=this%layer_physics(1)%pore_pressure_configuration(), &
                pore_cache=pore_pressure_cache_l1)
        ELSE
            CALL this%time_integration%advance(this%state,step_start,step_length,                                      &
                this%stochastic,this%domain,this%spatial_operator,                                              &
                pre_spatial_primitives_ready=pre_spatial_primitives_ready,                                            &
                thermal=this%layer_physics(1)%thermal_configuration(), &
                rheology=this%layer_physics(1)%rheology_configuration(), &
                pore_pressure=this%layer_physics(1)%pore_pressure_configuration(), &
                pore_cache=pore_pressure_cache_l1)
        END IF
        IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_LIE_COMPLETE,step_start,step_length)

        this%runtime%t=step_start+step_length
        IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_CLOCK_ADVANCED,                                            &
            this%runtime%t,step_length)
    CONTAINS
        SUBROUTINE notify_stochastic_prepared()
            ! Report the new phase at its actual point, after spatial.
            CALL observer(TIMESTEP_EVENT_STOCHASTIC_PREPARED,step_start,step_length)
        END SUBROUTINE notify_stochastic_prepared
    END SUBROUTINE advance_simulation_timestep

    !> \brief Release the simulation state, metadata and solver workspaces.
    !>
    !> \param[in,out] this Simulation-owned state, configuration descriptors and numerical
    !> workspaces.

    SUBROUTINE finalize_simulation(this)

        CLASS(simulation_context_type), INTENT(INOUT) :: this

        CALL this%domain%finalize
        CALL this%time_integration%finalize
        CALL this%spatial_operator%finalize

        CALL this%state%finalize


        CALL this%source_registry%finalize()
        IF (ALLOCATED(this%layer_physics)) DEALLOCATE(this%layer_physics)
        CALL this%model_layout%finalize
        CALL this%stochastic%finalize

    END SUBROUTINE finalize_simulation

END MODULE simulation_2d
