!> Compatibility API for the former specialized post-RK mass exchange entry point.
!> STEP dispatch and conservative updates are owned by step_source_executor_2d.
MODULE mass_exchange_step_2d
  USE step_source_executor_2d, ONLY : &
       apply_mass_exchange_step => execute_step_source_groups, &
       apply_mass_exchange_increment
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: apply_mass_exchange_step, apply_mass_exchange_increment
END MODULE mass_exchange_step_2d
