!********************************************************************************
!> \brief Simulation component aggregation and lifecycle
!
!> This module owns the state, runtime data and numerical workspaces that form
!> one simulation and coordinates their lifecycle.
!>
!> Owns state, runtime, domain, layout, equation partition, source registry and numerical workspaces. Initialization
!> follows their dependency order; finalization releases their storage.
!********************************************************************************

MODULE simulation_2d

  USE mass_injection_coefficients_2d, ONLY : set_source_time_coefficients
  USE step_source_executor_2d, ONLY : execute_step_source_groups

  USE source_registry_2d, ONLY : source_registry_type
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE diagnostics_2d, ONLY : fatal_error

  USE equation_metadata_2d, ONLY : equation_partition_type
  USE model_layout_2d, ONLY : model_layout_type

  USE parameters_2d, ONLY : wp, n_layers, n_vars, n_eqns, n_solid, n_add_gas, &
       n_stoch_vars, n_pore_vars, liquid_flag, gas_flag, stochastic_flag, &
       configure_canonical_state_indices

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                &
       finalize_nonlinear_solver

  USE geometry_2d, ONLY : release_topography_workspace

  USE runtime_2d, ONLY : runtime_state_type
  USE state_2d, ONLY : state_type
  USE domain_2d, ONLY : domain_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE stochastic_module, ONLY : stochastic_workspace_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: simulation_context_type
  PUBLIC :: simulation
  PUBLIC :: TIMESTEP_EVENT_SOURCE_WINDOWS, TIMESTEP_EVENT_STOCHASTIC_PREPARED
  PUBLIC :: TIMESTEP_EVENT_RK_COMPLETE, TIMESTEP_EVENT_STEP_SOURCES_COMPLETE
  PUBLIC :: TIMESTEP_EVENT_CLOCK_ADVANCED

  INTEGER, PARAMETER :: TIMESTEP_EVENT_SOURCE_WINDOWS=1
  INTEGER, PARAMETER :: TIMESTEP_EVENT_STOCHASTIC_PREPARED=2
  INTEGER, PARAMETER :: TIMESTEP_EVENT_RK_COMPLETE=3
  INTEGER, PARAMETER :: TIMESTEP_EVENT_STEP_SOURCES_COMPLETE=4
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
     TYPE(source_registry_type) :: source_registry !< Validated metadata for stage dispatch and source scheduling.
     TYPE(equation_partition_type) :: equation_partition
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
  !> \param[in,out] this Simulation-owned state, configuration descriptors and numerical workspaces.

  SUBROUTINE initialize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CHARACTER(LEN=256) :: source_error

    ! Canonical state indices must be available before any initial-condition
    ! conversion, including collapsing-volume initialization.
    CALL configure_canonical_state_indices()

    CALL this%model_layout%initialize( n_layers, n_vars, n_eqns, n_solid,   &
         n_add_gas, n_stoch_vars, n_pore_vars )

    CALL initialize_configured_sources(this%source_registry, &
         this%equation_partition,source_error)
    IF (LEN_TRIM(source_error)>0) CALL fatal_error(TRIM(source_error))

    CALL this%state%initialize( this%model_layout )

    CALL this%spatial_operator%initialize
    CALL this%domain%initialize

    CALL initialize_nonlinear_solver
    CALL this%time_integration%initialize

    CALL this%stochastic%initialize

    WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'

  END SUBROUTINE initialize_simulation

  !> \brief Apply one production timestep after the driver has selected dt.
  !> \details Preserves the fractional-step ordering: source-window averages,
  !>          OU preparation, RK-IMEX, post-RK step sources, then clock advance.
  !> \param[in,out] this Initialized simulation context.
  !> \param[in] observer Optional test hook called after each completed phase.
  SUBROUTINE advance_simulation_timestep(this,observer)
    CLASS(simulation_context_type), INTENT(INOUT) :: this
    PROCEDURE(timestep_event_observer), OPTIONAL :: observer
    REAL(wp) :: step_start,step_length

    step_start=this%runtime%t
    step_length=this%runtime%dt

    CALL set_source_time_coefficients(step_start,step_length)
    IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_SOURCE_WINDOWS,step_start,step_length)

    IF (stochastic_flag) CALL this%stochastic%prepare_timestep(this%state,step_length)
    IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_STOCHASTIC_PREPARED, &
         step_start,step_length)

    CALL this%time_integration%advance(this%state%q,this%state%qp,step_start, &
         step_length,this%stochastic%effective_Z,this%equation_partition, &
         this%source_registry,this%domain,this%spatial_operator)
    IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_RK_COMPLETE,step_start,step_length)

    CALL execute_step_source_groups(this%state%q,this%state%qp,step_length, &
         this%domain,this%source_registry,this%stochastic%effective_Z,step_start)
    IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_STEP_SOURCES_COMPLETE, &
         step_start,step_length)

    this%runtime%t=step_start+step_length
    IF (PRESENT(observer)) CALL observer(TIMESTEP_EVENT_CLOCK_ADVANCED, &
         this%runtime%t,step_length)
  END SUBROUTINE advance_simulation_timestep

  !> \brief Release the simulation state, metadata and solver workspaces.
  !>
  !> \param[in,out] this Simulation-owned state, configuration descriptors and numerical workspaces.

  SUBROUTINE finalize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CALL this%domain%finalize
    CALL this%time_integration%finalize
    CALL this%spatial_operator%finalize

    CALL this%state%finalize

    CALL finalize_nonlinear_solver

    CALL this%source_registry%finalize()
    CALL this%equation_partition%finalize
    CALL this%model_layout%finalize
    CALL this%stochastic%finalize
    CALL release_topography_workspace

  END SUBROUTINE finalize_simulation

END MODULE simulation_2d
