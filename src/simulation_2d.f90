!********************************************************************************
!> \brief Simulation component aggregation and lifecycle
!
!> This module owns the state, runtime data and numerical workspaces that form
!> one simulation and coordinates their lifecycle.
!********************************************************************************
MODULE simulation_2d

  USE constitutive_2d, ONLY : init_problem_param, finalize_problem_param

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                &
       finalize_nonlinear_solver

  USE runtime_2d, ONLY : runtime_state_type
  USE state_2d, ONLY : state_type
  USE domain_2d, ONLY : domain_type
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE stochastic_module, ONLY : stochastic_workspace_type

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: simulation_context_type
  PUBLIC :: simulation

  TYPE :: simulation_context_type

     TYPE(runtime_state_type) :: runtime
     TYPE(state_type) :: state
     TYPE(domain_type) :: domain
     TYPE(reconstruction_workspace_type) :: reconstruction
     TYPE(hyperbolic_workspace_type) :: hyperbolic
     TYPE(time_integration_workspace_type) :: time_integration
     TYPE(stochastic_workspace_type) :: stochastic

   CONTAINS

     PROCEDURE :: initialize => initialize_simulation
     PROCEDURE :: finalize => finalize_simulation

  END TYPE simulation_context_type

  TYPE(simulation_context_type) :: simulation

CONTAINS

  SUBROUTINE initialize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CALL init_problem_param

    CALL this%state%initialize

    CALL this%reconstruction%initialize
    CALL this%hyperbolic%initialize
    CALL this%domain%initialize

    CALL initialize_nonlinear_solver
    CALL this%time_integration%initialize

    CALL this%stochastic%initialize

    WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'

  END SUBROUTINE initialize_simulation

  SUBROUTINE finalize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CALL this%domain%finalize
    CALL this%time_integration%finalize
    CALL this%hyperbolic%finalize
    CALL this%reconstruction%finalize

    CALL this%state%finalize

    CALL finalize_nonlinear_solver

    CALL finalize_problem_param
    CALL this%stochastic%finalize

  END SUBROUTINE finalize_simulation

END MODULE simulation_2d
