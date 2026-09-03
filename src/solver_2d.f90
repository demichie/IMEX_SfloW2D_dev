!********************************************************************************
!> \brief Numerical component lifecycle orchestration
!
!> This module initializes and finalizes the components used by the solver.
!> Scientific state and transient workspaces are owned by those components.
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************
MODULE solver_2d

  USE constitutive_2d, ONLY : init_problem_param, finalize_problem_param

  USE state_2d, ONLY : state

  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver,                 &
       finalize_nonlinear_solver

  USE reconstruction_2d, ONLY : reconstruction_workspace

  USE hyperbolic_2d, ONLY : hyperbolic_workspace

  USE domain_2d, ONLY : domain
  USE time_integration_2d, ONLY : time_integration_workspace
  USE stochastic_module, ONLY : stochastic_workspace

  IMPLICIT none

CONTAINS

  !******************************************************************************
  !> \brief Memory allocation
  !
  !> This subroutine allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE allocate_solver_variables
    
    IMPLICIT NONE

    CALL init_problem_param

    CALL state%initialize

    CALL reconstruction_workspace%initialize
    CALL hyperbolic_workspace%initialize
    CALL domain%initialize

    CALL initialize_nonlinear_solver
    CALL time_integration_workspace%initialize

    CALL stochastic_workspace%initialize
    
    WRITE(*,*) 'ALLOCATION OF ARRAYS COMPLETED'
    
    RETURN
    
  END SUBROUTINE allocate_solver_variables

  !******************************************************************************
  !> \brief Memory deallocation
  !
  !> This subroutine de-allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE deallocate_solver_variables

    CALL domain%finalize
    CALL time_integration_workspace%finalize
    CALL hyperbolic_workspace%finalize
    CALL reconstruction_workspace%finalize

    CALL state%finalize

    CALL finalize_nonlinear_solver

    CALL finalize_problem_param
    CALL stochastic_workspace%finalize
    
    RETURN
    
  END SUBROUTINE deallocate_solver_variables


END MODULE solver_2d
