PROGRAM test_hp_reconstruction

  USE parameters_2d, ONLY : wp, theta
  USE geometry_2d, ONLY : limit
  USE hp_reconstruction_2d, ONLY : reconstruct_hp_line
  USE hp_reconstruction_2d, ONLY : hp_dry_tolerance

  IMPLICIT NONE

  theta = 1.3_wp

  CALL check_lake_at_rest
  CALL check_shoreline_positivity
  CALL check_momentum_admissibility
  CALL check_flat_bed_non_regression
  CALL check_direction_independence

  WRITE(*,*) 'PASS: HP reconstruction invariants verified'

CONTAINS

  SUBROUTINE check_lake_at_rest

    INTEGER, PARAMETER :: n = 7
    REAL(wp), PARAMETER :: H0 = 3.0_wp
    REAL(wp) :: B_face(n+1), Bm(n), Bp(n), h(n), u(n)
    REAL(wp) :: hm0(n), hp0(n), hum0(n), hup0(n), um0(n), up0(n)
    REAL(wp) :: hm(n), hp(n), hum(n), hup(n), etam(n), etap(n), w(n)
    REAL(wp) :: tolerance

    B_face = [ 0.0_wp, 0.2_wp, 0.45_wp, 1.35_wp, 1.35_wp,                  &
         0.85_wp, 0.65_wp, 0.55_wp ]
    Bm = B_face(1:n)
    Bp = B_face(2:n+1)
    h = H0 - 0.5_wp*(Bm+Bp)
    u = 0.0_wp

    CALL direct_candidates(h,hm0,hp0,3,1.0_wp)
    hum0 = 0.0_wp
    hup0 = 0.0_wp
    um0 = 0.0_wp
    up0 = 0.0_wp

    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm,hp,hum,hup,etam,etap,w)

    tolerance = 256.0_wp*EPSILON(1.0_wp)*H0
    CALL assert_small('lake-at-rest west eta',MAXVAL(ABS(etam-H0)),tolerance)
    CALL assert_small('lake-at-rest east eta',MAXVAL(ABS(etap-H0)),tolerance)
    CALL assert_small('lake-at-rest west h',MAXVAL(ABS(hm+Bm-H0)),tolerance)
    CALL assert_small('lake-at-rest east h',MAXVAL(ABS(hp+Bp-H0)),tolerance)
    CALL assert_small('lake-at-rest momentum',                              &
         MAX(MAXVAL(ABS(hum)),MAXVAL(ABS(hup))),tolerance)

  END SUBROUTINE check_lake_at_rest

  SUBROUTINE check_shoreline_positivity

    INTEGER, PARAMETER :: n = 6
    REAL(wp), PARAMETER :: H0 = 1.4_wp
    REAL(wp) :: B_face(n+1), Bm(n), Bp(n), h(n), u(n)
    REAL(wp) :: hm0(n), hp0(n), hum0(n), hup0(n), um0(n), up0(n)
    REAL(wp) :: hm(n), hp(n), hum(n), hup(n), etam(n), etap(n), w(n)

    B_face = [ 0.0_wp, 0.4_wp, 0.9_wp, 1.35_wp, 1.55_wp, 1.8_wp, 2.1_wp ]
    Bm = B_face(1:n)
    Bp = B_face(2:n+1)
    h = MAX(H0-0.5_wp*(Bm+Bp),0.0_wp)
    u = 0.0_wp

    CALL direct_candidates(h,hm0,hp0,3,1.0_wp)
    hum0 = 0.0_wp
    hup0 = 0.0_wp
    um0 = 0.0_wp
    up0 = 0.0_wp

    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm,hp,hum,hup,etam,etap,w)

    IF ( MIN(MINVAL(hm),MINVAL(hp)) .LT. 0.0_wp ) THEN
       WRITE(*,*) 'FAIL: shoreline reconstruction produced negative h'
       ERROR STOP 1
    END IF
    IF ( ANY( ( h .LE. hp_dry_tolerance ) .AND. ( w .NE. 1.0_wp ) ) ) THEN
       WRITE(*,*) 'FAIL: dry shoreline cells did not select eta reconstruction'
       ERROR STOP 1
    END IF

  END SUBROUTINE check_shoreline_positivity

  SUBROUTINE check_momentum_admissibility

    INTEGER, PARAMETER :: n = 6
    REAL(wp) :: Bm(n), Bp(n), h(n), u(n)
    REAL(wp) :: hm0(n), hp0(n), hum0(n), hup0(n), um0(n), up0(n)
    REAL(wp) :: hm(n), hp(n), hum(n), hup(n), etam(n), etap(n), w(n)
    REAL(wp) :: umin, umax, tolerance, mean_error, velocity_error
    INTEGER :: i

    Bm = [ 0.0_wp, 0.1_wp, 0.2_wp, 0.8_wp, 0.9_wp, 0.7_wp ]
    Bp = [ 0.1_wp, 0.2_wp, 0.8_wp, 0.9_wp, 0.7_wp, 0.6_wp ]
    h = [ 1.2_wp, 1.0_wp, 0.85_wp, 0.7_wp, 0.9_wp, 1.1_wp ]
    u = [ -0.4_wp, -0.1_wp, 0.6_wp, 1.1_wp, 0.3_wp, -0.2_wp ]
    hm0 = h - [ 0.0_wp, 0.1_wp, 0.15_wp, 0.2_wp, 0.1_wp, 0.0_wp ]
    hp0 = 2.0_wp*h-hm0
    hum0 = h*u - [ 0.0_wp, 0.8_wp, -0.9_wp, 1.1_wp, -0.7_wp, 0.0_wp ]
    hup0 = 2.0_wp*h*u-hum0
    um0 = u - 0.35_wp
    up0 = u + 0.35_wp

    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm,hp,hum,hup,etam,etap,w)

    tolerance = 1024.0_wp*EPSILON(1.0_wp)
    mean_error = MAXVAL(ABS(0.5_wp*(hum+hup)-h*u))
    CALL assert_small('normal momentum mean',mean_error,tolerance)
    CALL assert_small('thickness mean',MAXVAL(ABS(0.5_wp*(hm+hp)-h)),tolerance)

    velocity_error = 0.0_wp
    DO i = 1, n
       umin = MINVAL(u(MAX(1,i-1):MIN(n,i+1)))
       umax = MAXVAL(u(MAX(1,i-1):MIN(n,i+1)))
       IF (hm(i) .GT. hp_dry_tolerance) THEN
          velocity_error = MAX(velocity_error,umin-hum(i)/hm(i),             &
               hum(i)/hm(i)-umax)
       END IF
       IF (hp(i) .GT. hp_dry_tolerance) THEN
          velocity_error = MAX(velocity_error,umin-hup(i)/hp(i),             &
               hup(i)/hp(i)-umax)
       END IF
    END DO
    CALL assert_small('normal velocity interval',MAX(velocity_error,0.0_wp), &
         tolerance)

  END SUBROUTINE check_momentum_admissibility

  SUBROUTINE check_flat_bed_non_regression

    INTEGER, PARAMETER :: n = 5
    REAL(wp) :: Bm(n), Bp(n), h(n), u(n)
    REAL(wp) :: hm0(n), hp0(n), hum0(n), hup0(n), um0(n), up0(n)
    REAL(wp) :: hm(n), hp(n), hum(n), hup(n), etam(n), etap(n), w(n)
    REAL(wp) :: hu(n), tolerance

    Bm = 2.0_wp
    Bp = 2.0_wp
    h = 1.25_wp
    u = [ -0.4_wp, -0.2_wp, 0.0_wp, 0.2_wp, 0.4_wp ]
    hu = h*u
    CALL direct_candidates(h,hm0,hp0,3,1.0_wp)
    CALL direct_candidates(hu,hum0,hup0,3,1.0_wp)
    CALL direct_candidates(u,um0,up0,3,1.0_wp)

    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm,hp,hum,hup,etam,etap,w)

    tolerance = 512.0_wp*EPSILON(1.0_wp)
    CALL assert_small('flat-bed h west',MAXVAL(ABS(hm-hm0)),tolerance)
    CALL assert_small('flat-bed h east',MAXVAL(ABS(hp-hp0)),tolerance)
    CALL assert_small('flat-bed hu west',MAXVAL(ABS(hum-hum0)),tolerance)
    CALL assert_small('flat-bed hu east',MAXVAL(ABS(hup-hup0)),tolerance)

  END SUBROUTINE check_flat_bed_non_regression

  SUBROUTINE check_direction_independence

    INTEGER, PARAMETER :: n = 3
    REAL(wp) :: Bm(n), Bp(n), h(n), u(n), hm0(n), hp0(n)
    REAL(wp) :: hum0(n), hup0(n), um0(n), up0(n)
    REAL(wp) :: hm_x(n), hp_x(n), hum_x(n), hup_x(n), em_x(n), ep_x(n), w_x(n)
    REAL(wp) :: hm_y(n), hp_y(n), hum_y(n), hup_y(n), em_y(n), ep_y(n), w_y(n)

    Bm = [ 0.2_wp, 0.7_wp, 0.9_wp ]
    Bp = [ 0.7_wp, 0.9_wp, 0.4_wp ]
    h = [ 1.1_wp, 0.8_wp, 1.0_wp ]
    u = [ -0.2_wp, 0.4_wp, 0.1_wp ]
    CALL direct_candidates(h,hm0,hp0,3,1.0_wp)
    CALL direct_candidates(h*u,hum0,hup0,3,1.0_wp)
    CALL direct_candidates(u,um0,up0,3,1.0_wp)

    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm_x,hp_x,hum_x,hup_x,em_x,ep_x,w_x)
    CALL reconstruct_hp_line(h,u,Bm,Bp,hm0,hp0,hum0,hup0,um0,up0,3,        &
         1.0_wp,hm_y,hp_y,hum_y,hup_y,em_y,ep_y,w_y)

    CALL assert_small('row/column HP core',MAX(                             &
         MAXVAL(ABS(hm_x-hm_y)),MAXVAL(ABS(hp_x-hp_y)),                    &
         MAXVAL(ABS(hum_x-hum_y)),MAXVAL(ABS(hup_x-hup_y)),                &
         MAXVAL(ABS(em_x-em_y)),MAXVAL(ABS(ep_x-ep_y)),                    &
         MAXVAL(ABS(w_x-w_y))),0.0_wp)

  END SUBROUTINE check_direction_independence

  SUBROUTINE direct_candidates(values,value_minus,value_plus,limiter_id,coeff)

    REAL(wp), INTENT(IN) :: values(:), coeff
    REAL(wp), INTENT(OUT) :: value_minus(:), value_plus(:)
    INTEGER, INTENT(IN) :: limiter_id
    REAL(wp) :: stencil(3), coordinates(3), slope
    INTEGER :: i, n

    n = SIZE(values)
    coordinates = [ -1.0_wp, 0.0_wp, 1.0_wp ]
    value_minus = values
    value_plus = values
    DO i = 2, n-1
       stencil = values(i-1:i+1)
       CALL limit(stencil,coordinates,limiter_id,slope)
       value_minus(i) = values(i)-0.5_wp*coeff*slope
       value_plus(i) = values(i)+0.5_wp*coeff*slope
    END DO

  END SUBROUTINE direct_candidates

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF ( value .GT. limit_value ) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,' > ',limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_hp_reconstruction
