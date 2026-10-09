!********************************************************************************
!> \brief Runtime state ownership
!
!> This module owns the mutable time state of a simulation.
!>
!> Groups the current simulation time and timestep history independently of physical state arrays.
!> Binary restarts preserve this history together with the output schedule.
!********************************************************************************

MODULE runtime_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: runtime_state_type

  !> \brief Mutable clock, timestep history and runout bookkeeping for one run.
  TYPE :: runtime_state_type

     !> Current simulation time
     REAL(wp) :: t = 0.0_wp

     !> Current time step
     REAL(wp) :: dt = 0.0_wp

     !> Time steps used by the temporal smoothing limiter
     REAL(wp) :: dt_old = 0.0_wp
     REAL(wp) :: dt_old_old = 0.0_wp

  END TYPE runtime_state_type

END MODULE runtime_2d
