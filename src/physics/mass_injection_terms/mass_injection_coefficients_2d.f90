!> \brief Per-step source-window coefficients shared by injection laws.
!> \details Updated before the RK stages by the existing window initializer.
MODULE mass_injection_coefficients_2d
  USE parameters_2d, ONLY : wp
  IMPLICIT NONE
  REAL(wp) :: bottom_coeff_step = -1.0_wp
  REAL(wp) :: bottom_vel_step = -1.0_wp
  REAL(wp) :: lateral_coeff_step = -1.0_wp
  REAL(wp), ALLOCATABLE :: fissure_coeff_step(:)
END MODULE mass_injection_coefficients_2d
