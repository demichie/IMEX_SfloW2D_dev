PROGRAM test_curvature_source

  USE parameters_2d, ONLY : wp
  USE equation_terms_2d, ONLY : curvature_acceleration
  USE equation_terms_2d, ONLY : eval_curvature_momentum_source

  IMPLICIT NONE

  REAL(wp), PARAMETER :: Bx = 0.35_wp
  REAL(wp), PARAMETER :: By = -0.22_wp
  REAL(wp), PARAMETER :: Bxx = 0.08_wp
  REAL(wp), PARAMETER :: Bxy = -0.035_wp
  REAL(wp), PARAMETER :: Byy = 0.06_wp
  REAL(wp), PARAMETER :: G = 0.42_wp
  REAL(wp), PARAMETER :: mixture_mass = 1800.0_wp
  REAL(wp), PARAMETER :: u = 2.4_wp
  REAL(wp), PARAMETER :: v = -1.7_wp

  REAL(wp) :: acceleration, expected_acceleration
  REAL(wp) :: source_x, source_y, expected_x, expected_y
  REAL(wp) :: tolerance

  expected_acceleration = Bxx*u*u + 2.0_wp*Bxy*u*v + Byy*v*v
  CALL assert_small( 'curvature acceleration',                             &
       ABS(curvature_acceleration(Bxx,Bxy,Byy,u,v)-expected_acceleration),  &
       64.0_wp*EPSILON(1.0_wp)*ABS(expected_acceleration) )

  acceleration = curvature_acceleration(Bxx,Bxy,Byy,u,v)
  CALL eval_curvature_momentum_source( Bx, By, Bxx, Bxy, Byy, G,           &
       mixture_mass, u, v, source_x, source_y )

  expected_x = -G*mixture_mass*expected_acceleration*Bx
  expected_y = -G*mixture_mass*expected_acceleration*By
  tolerance = 128.0_wp*EPSILON(1.0_wp)                                    &
       * MAX(1.0_wp,ABS(expected_x),ABS(expected_y))

  CALL assert_small('returned acceleration',                              &
       ABS(acceleration-expected_acceleration),tolerance)
  CALL assert_small('x curvature source',ABS(source_x-expected_x),tolerance)
  CALL assert_small('y curvature source',ABS(source_y-expected_y),tolerance)

  ! The mixed-Hessian contribution must not be dropped in two dimensions.
  CALL assert_small( 'mixed Hessian contribution',                         &
       ABS(expected_acceleration-(Bxx*u*u+Byy*v*v)-2.0_wp*Bxy*u*v),         &
       64.0_wp*EPSILON(1.0_wp)*ABS(expected_acceleration) )

  WRITE(*,*) 'PASS: slope-weighted curvature source verified'

CONTAINS

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF ( value .GT. limit_value ) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,' > ',limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_curvature_source
