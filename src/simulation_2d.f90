!********************************************************************************
!> \brief Simulation component aggregation and lifecycle
!
!> This module owns the state, runtime data and numerical workspaces that form
!> one simulation and coordinates their lifecycle.
!>
!> Owns state, runtime, domain, layout, equation partition and numerical workspaces. Initialization
!> follows their dependency order; finalization releases their storage.
!********************************************************************************

MODULE simulation_2d

  USE equation_terms_2d, ONLY : init_problem_param

  USE equation_metadata_2d, ONLY : equation_partition_type
  USE model_layout_2d, ONLY : model_layout_type

  USE parameters_2d, ONLY : n_layers, n_vars, n_eqns, n_solid, n_add_gas,   &
       n_stoch_vars, n_pore_vars, liquid_flag, gas_flag

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                &
       finalize_nonlinear_solver

  USE mass_exchange_2d, ONLY : release_topography_workspace

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

  !> \brief Ownership root for runtime, states and reusable numerical workspaces.
  !> \details Initialization follows the model-layout dependencies; finalization
  !>          releases the owned storage without moving it into global facades.
  TYPE :: simulation_context_type

     TYPE(runtime_state_type) :: runtime
     TYPE(model_layout_type) :: model_layout
     TYPE(equation_partition_type) :: equation_partition
     TYPE(state_type) :: state
     TYPE(domain_type) :: domain
     TYPE(spatial_operator_type) :: spatial_operator
     TYPE(time_integration_workspace_type) :: time_integration
     TYPE(stochastic_workspace_type) :: stochastic

   CONTAINS

     PROCEDURE :: initialize => initialize_simulation
     PROCEDURE :: finalize => finalize_simulation

  END TYPE simulation_context_type

  TYPE(simulation_context_type) :: simulation

CONTAINS

  !> \brief Initialize model metadata and simulation-owned numerical workspaces.
  !>
  !> \param[in,out] this Simulation-owned state, configuration descriptors and numerical workspaces.

  SUBROUTINE initialize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CALL this%model_layout%initialize( n_layers, n_vars, n_eqns, n_solid,   &
         n_add_gas, n_stoch_vars, n_pore_vars, liquid_flag .AND. gas_flag )

    CALL init_problem_param( this%equation_partition )

    CALL this%state%initialize( this%model_layout )

    CALL this%spatial_operator%initialize
    CALL this%domain%initialize

    CALL initialize_nonlinear_solver
    CALL this%time_integration%initialize

    CALL this%stochastic%initialize

    WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'

  END SUBROUTINE initialize_simulation

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

    CALL this%equation_partition%finalize
    CALL this%model_layout%finalize
    CALL this%stochastic%finalize
    CALL release_topography_workspace

  END SUBROUTINE finalize_simulation

END MODULE simulation_2d
