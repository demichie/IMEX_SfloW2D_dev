!********************************************************************************
!> \brief Simulation component aggregation and lifecycle
!
!> This module provides an explicit view of the state, runtime data and
!> numerical workspaces that form one simulation. The context refers to the
!> transitional module instances; numerical kernels can therefore be migrated
!> incrementally without duplicating their allocatable data.
!********************************************************************************
MODULE simulation_2d

  USE constitutive_2d, ONLY : init_problem_param, finalize_problem_param

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                &
       finalize_nonlinear_solver

  USE runtime_2d, ONLY : runtime_state_type, runtime_instance => runtime
  USE state_2d, ONLY : state_type, state_instance => state
  USE domain_2d, ONLY : domain_type, domain_instance => domain
  USE reconstruction_2d, ONLY : reconstruction_workspace_type,               &
       reconstruction_instance => reconstruction_workspace
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type,                        &
       hyperbolic_instance => hyperbolic_workspace
  USE time_integration_2d, ONLY : time_integration_workspace_type,            &
       time_integration_instance => time_integration_workspace
  USE stochastic_module, ONLY : stochastic_workspace_type,                   &
       stochastic_instance => stochastic_workspace

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: simulation_context_type
  PUBLIC :: simulation

  TYPE :: simulation_context_type

     TYPE(runtime_state_type), POINTER :: runtime => NULL()
     TYPE(state_type), POINTER :: state => NULL()
     TYPE(domain_type), POINTER :: domain => NULL()
     TYPE(reconstruction_workspace_type), POINTER :: reconstruction => NULL()
     TYPE(hyperbolic_workspace_type), POINTER :: hyperbolic => NULL()
     TYPE(time_integration_workspace_type), POINTER :: time_integration => NULL()
     TYPE(stochastic_workspace_type), POINTER :: stochastic => NULL()

   CONTAINS

     PROCEDURE, PRIVATE :: bind_components
     PROCEDURE :: initialize => initialize_simulation
     PROCEDURE :: finalize => finalize_simulation

  END TYPE simulation_context_type

  TYPE(simulation_context_type) :: simulation

CONTAINS

  SUBROUTINE bind_components(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    this%runtime => runtime_instance
    this%state => state_instance
    this%domain => domain_instance
    this%reconstruction => reconstruction_instance
    this%hyperbolic => hyperbolic_instance
    this%time_integration => time_integration_instance
    this%stochastic => stochastic_instance

  END SUBROUTINE bind_components

  SUBROUTINE initialize_simulation(this)

    CLASS(simulation_context_type), INTENT(INOUT) :: this

    CALL this%bind_components

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

    NULLIFY(this%runtime)
    NULLIFY(this%state)
    NULLIFY(this%domain)
    NULLIFY(this%reconstruction)
    NULLIFY(this%hyperbolic)
    NULLIFY(this%time_integration)
    NULLIFY(this%stochastic)

  END SUBROUTINE finalize_simulation

END MODULE simulation_2d
