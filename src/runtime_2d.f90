!********************************************************************************
!> \brief Runtime state ownership
!
!> This module owns the mutable time state of a simulation.
!********************************************************************************
MODULE runtime_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: runtime_state_type
  PUBLIC :: runtime

  TYPE :: runtime_state_type

     !> Current simulation time
     REAL(wp) :: t = 0.0_wp

     !> Current time step
     REAL(wp) :: dt = 0.0_wp

  END TYPE runtime_state_type

  TYPE(runtime_state_type) :: runtime

END MODULE runtime_2d
