PROGRAM test_pccu

  USE parameters_2d, ONLY : wp
  USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y
  USE pccu_2d, ONLY : eval_hydrostatic_path
  USE pccu_2d, ONLY : eval_central_upwind_flux
  USE pccu_2d, ONLY : eval_oriented_pccu_pair

  IMPLICIT NONE

  CALL check_hydrostatic_path
  CALL check_oriented_identity
  CALL check_degenerate_pair

  WRITE(*,*) 'PASS: G=1 HP-PCCU algebra verified'

CONTAINS

  SUBROUTINE check_hydrostatic_path

    REAL(wp) :: path_x(6), path_y(6)
    REAL(wp) :: hL, hR, gammaL, gammaR, etaL, etaR
    REAL(wp) :: expected, mean_h_squared, tolerance

    hL = 0.7_wp
    hR = 1.6_wp
    gammaL = 8450.0_wp
    gammaR = 9320.0_wp
    etaL = 2.4_wp
    etaR = 2.9_wp

    CALL eval_hydrostatic_path(PATH_DIR_X,hL,gammaL,etaL,hR,gammaR,etaR,path_x)
    CALL eval_hydrostatic_path(PATH_DIR_Y,hL,gammaL,etaL,hR,gammaR,etaR,path_y)

    mean_h_squared = (hL*hL+hL*hR+hR*hR)/3.0_wp
    expected = -(etaR-etaL) * ( gammaL*hL                            &
         + 0.5_wp*(gammaL*(hR-hL)+(gammaR-gammaL)*hL)                &
         + (gammaR-gammaL)*(hR-hL)/3.0_wp )                         &
         - 0.5_wp*(gammaR-gammaL)*mean_h_squared

    ! The first parenthesis is the exact integral of Gamma(s)*h(s).
    tolerance = 2048.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,ABS(expected))
    CALL assert_small('x hydrostatic momentum',ABS(path_x(2)-expected),tolerance)
    CALL assert_small('y hydrostatic momentum',ABS(path_y(3)-expected),tolerance)
    CALL assert_small('x non-normal components',                            &
         MAXVAL(ABS(path_x([1,3,4,5,6]))),tolerance)
    CALL assert_small('y non-normal components',                            &
         MAXVAL(ABS(path_y([1,2,4,5,6]))),tolerance)
    CALL assert_small('thermal path component x',ABS(path_x(4)),0.0_wp)
    CALL assert_small('thermal path component y',ABS(path_y(4)),0.0_wp)

  END SUBROUTINE check_hydrostatic_path

  SUBROUTINE check_oriented_identity

    REAL(wp) :: fluxL(6), fluxR(6), qL(6), qR(6)
    REAL(wp) :: H(6), path(6), valueL(6), valueR(6), expected(6)
    REAL(wp) :: a_minus, a_plus, tolerance

    fluxL = [ 1.0_wp, -2.0_wp, 0.5_wp, 4.0_wp, -0.7_wp, 0.2_wp ]
    fluxR = [ 0.3_wp, 1.4_wp, -0.8_wp, 2.0_wp, 0.1_wp, -0.4_wp ]
    qL = [ 2.0_wp, -0.5_wp, 0.4_wp, 8.0_wp, 0.2_wp, 0.1_wp ]
    qR = [ 1.5_wp, 0.7_wp, -0.2_wp, 6.0_wp, 0.1_wp, 0.3_wp ]
    path = 0.0_wp
    path(2) = -3.75_wp
    a_minus = -1.2_wp
    a_plus = 2.8_wp

    CALL eval_central_upwind_flux(a_minus,a_plus,fluxL,fluxR,qL,qR,H)
    expected = (a_plus*fluxL-a_minus*fluxR+a_plus*a_minus*(qR-qL))           &
         / (a_plus-a_minus)
    tolerance = 512.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,MAXVAL(ABS(expected)))
    CALL assert_small('central-upwind flux',MAXVAL(ABS(H-expected)),tolerance)

    CALL eval_oriented_pccu_pair(H,path,a_minus,a_plus,valueL,valueR)
    CALL assert_small('oriented path identity',                              &
         MAXVAL(ABS(valueR-valueL-path)),tolerance)
    CALL assert_small('oriented equation 4',ABS(valueR(4)-valueL(4)),0.0_wp)

  END SUBROUTINE check_oriented_identity

  SUBROUTINE check_degenerate_pair

    REAL(wp) :: fluxL(4), fluxR(4), qL(4), qR(4)
    REAL(wp) :: H(4), path(4), valueL(4), valueR(4)

    fluxL = [ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp ]
    fluxR = [ 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp ]
    qL = 0.0_wp
    qR = 1.0_wp
    path = [ 0.0_wp, -2.0_wp, 0.0_wp, 0.0_wp ]

    CALL eval_central_upwind_flux(0.0_wp,0.0_wp,fluxL,fluxR,qL,qR,H)
    CALL assert_small('degenerate central flux',                            &
         MAXVAL(ABS(H-0.5_wp*(fluxL+fluxR))),0.0_wp)
    CALL eval_oriented_pccu_pair(H,path,0.0_wp,0.0_wp,valueL,valueR)
    CALL assert_small('degenerate oriented pair',MAXVAL(ABS(valueL-H)),0.0_wp)
    CALL assert_small('degenerate pair equality',MAXVAL(ABS(valueR-valueL)),0.0_wp)

  END SUBROUTINE check_degenerate_pair

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF ( value .GT. limit_value ) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,' > ',limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_pccu
