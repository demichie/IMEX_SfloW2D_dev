!********************************************************************************
!> \brief Hydrostatic/positivity-preserving reconstruction helpers
!>
!> This module contains the direction-independent scalar core of the HP
!> reconstruction.  The same routine is applied to grid rows and columns by
!> reconstruction_2d; fluxes and source terms deliberately remain outside this
!> module.
!********************************************************************************
MODULE hp_reconstruction_2d

  USE parameters_2d, ONLY : wp
  USE geometry_2d, ONLY : limit

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER, PUBLIC :: hp_dry_tolerance = 1.0E-10_wp
  REAL(wp), PARAMETER, PUBLIC :: hp_bed_step_tolerance = 1.0E-12_wp

  PUBLIC :: reconstruct_hp_line

CONTAINS

  !******************************************************************************
  !> \brief Apply the HP reconstruction to one row or column of cells
  !>
  !> Direct reconstructions of h, h*u_n and u_n are supplied by the caller.
  !> This routine reconstructs eta=h+B, enforces non-negative endpoint depths,
  !> computes the parameter-free continuity blend, and limits the final normal
  !> volumetric momentum while preserving its cell mean.
  !>
  !> Boundary eta slopes use the frozen-reference zero-gradient convention.
  !******************************************************************************
  SUBROUTINE reconstruct_hp_line( h_center, u_center, B_minus, B_plus,          &
       h_minus_direct, h_plus_direct, hu_minus_direct, hu_plus_direct,          &
       u_minus_candidate, u_plus_candidate, limiter_id,                        &
       reconstruction_coefficient, h_minus, h_plus, hu_minus, hu_plus,         &
       eta_minus, eta_plus, w_eta )

    REAL(wp), INTENT(IN) :: h_center(:)
    REAL(wp), INTENT(IN) :: u_center(:)
    REAL(wp), INTENT(IN) :: B_minus(:)
    REAL(wp), INTENT(IN) :: B_plus(:)
    REAL(wp), INTENT(IN) :: h_minus_direct(:)
    REAL(wp), INTENT(IN) :: h_plus_direct(:)
    REAL(wp), INTENT(IN) :: hu_minus_direct(:)
    REAL(wp), INTENT(IN) :: hu_plus_direct(:)
    REAL(wp), INTENT(IN) :: u_minus_candidate(:)
    REAL(wp), INTENT(IN) :: u_plus_candidate(:)
    INTEGER, INTENT(IN) :: limiter_id
    REAL(wp), INTENT(IN) :: reconstruction_coefficient

    REAL(wp), INTENT(OUT) :: h_minus(:)
    REAL(wp), INTENT(OUT) :: h_plus(:)
    REAL(wp), INTENT(OUT) :: hu_minus(:)
    REAL(wp), INTENT(OUT) :: hu_plus(:)
    REAL(wp), INTENT(OUT) :: eta_minus(:)
    REAL(wp), INTENT(OUT) :: eta_plus(:)
    REAL(wp), INTENT(OUT) :: w_eta(:)

    REAL(wp) :: eta_center(SIZE(h_center))
    REAL(wp) :: h_minus_eta(SIZE(h_center))
    REAL(wp) :: h_plus_eta(SIZE(h_center))
    REAL(wp) :: h_minus_orig(SIZE(h_center))
    REAL(wp) :: h_plus_orig(SIZE(h_center))
    REAL(wp) :: Eh(SIZE(h_center)), Eeta(SIZE(h_center))
    REAL(wp) :: eta_stencil(3), coordinate_stencil(3)
    REAL(wp) :: eta_slope, slope_min, slope_max
    REAL(wp) :: denominator, distribution_tolerance
    REAL(wp) :: momentum_weight, dm_target, dm_lower, dm_upper, dm
    REAL(wp) :: u_min, u_max
    INTEGER :: i, number_of_cells

    number_of_cells = SIZE(h_center)

    IF ( number_of_cells .LE. 0 ) RETURN

    IF ( ( SIZE(u_center) .NE. number_of_cells ) .OR.                         &
         ( SIZE(B_minus) .NE. number_of_cells ) .OR.                          &
         ( SIZE(B_plus) .NE. number_of_cells ) .OR.                           &
         ( SIZE(h_minus_direct) .NE. number_of_cells ) .OR.                   &
         ( SIZE(h_plus_direct) .NE. number_of_cells ) .OR.                    &
         ( SIZE(hu_minus_direct) .NE. number_of_cells ) .OR.                  &
         ( SIZE(hu_plus_direct) .NE. number_of_cells ) .OR.                   &
         ( SIZE(u_minus_candidate) .NE. number_of_cells ) .OR.                &
         ( SIZE(u_plus_candidate) .NE. number_of_cells ) .OR.                 &
         ( SIZE(h_minus) .NE. number_of_cells ) .OR.                          &
         ( SIZE(h_plus) .NE. number_of_cells ) .OR.                           &
         ( SIZE(hu_minus) .NE. number_of_cells ) .OR.                         &
         ( SIZE(hu_plus) .NE. number_of_cells ) .OR.                          &
         ( SIZE(eta_minus) .NE. number_of_cells ) .OR.                        &
         ( SIZE(eta_plus) .NE. number_of_cells ) .OR.                         &
         ( SIZE(w_eta) .NE. number_of_cells ) ) THEN
       ERROR STOP 'reconstruct_hp_line: inconsistent array sizes'
    END IF

    eta_center = h_center + 0.5_wp * ( B_minus + B_plus )
    coordinate_stencil = [ -1.0_wp, 0.0_wp, 1.0_wp ]

    DO i = 1, number_of_cells

       eta_slope = 0.0_wp
       IF ( ( i .GT. 1 ) .AND. ( i .LT. number_of_cells ) ) THEN
          eta_stencil = eta_center(i-1:i+1)
          CALL limit( eta_stencil, coordinate_stencil, limiter_id, eta_slope )
          eta_slope = reconstruction_coefficient * eta_slope
       END IF

       slope_min = 2.0_wp * ( B_plus(i) - eta_center(i) )
       slope_max = 2.0_wp * ( eta_center(i) - B_minus(i) )
       eta_slope = MIN( MAX( eta_slope, slope_min ), slope_max )

       eta_minus(i) = eta_center(i) - 0.5_wp * eta_slope
       eta_plus(i) = eta_center(i) + 0.5_wp * eta_slope
       h_minus_eta(i) = MAX( eta_minus(i) - B_minus(i), 0.0_wp )
       h_plus_eta(i) = MAX( eta_plus(i) - B_plus(i), 0.0_wp )

    END DO

    h_minus_orig = MAX( h_minus_direct, 0.0_wp )
    h_plus_orig = MAX( h_plus_direct, 0.0_wp )

    Eh = 0.0_wp
    Eeta = 0.0_wp
    DO i = 1, number_of_cells-1
       Eh(i) = Eh(i) + ABS( h_plus_orig(i) - h_minus_orig(i+1) )
       Eh(i+1) = Eh(i+1) + ABS( h_plus_orig(i) - h_minus_orig(i+1) )
       Eeta(i) = Eeta(i) + ABS( h_plus_eta(i) - h_minus_eta(i+1) )
       Eeta(i+1) = Eeta(i+1) + ABS( h_plus_eta(i) - h_minus_eta(i+1) )
    END DO

    DO i = 1, number_of_cells

       denominator = Eh(i) + Eeta(i)
       distribution_tolerance = 1.0E-14_wp * MAX( 1.0_wp, h_center(i) )
       IF ( denominator .GT. distribution_tolerance ) THEN
          w_eta(i) = Eh(i) / denominator
       ELSE
          w_eta(i) = 0.5_wp
       END IF

       IF ( ( h_minus_eta(i) .LE. hp_dry_tolerance ) .OR.                    &
            ( h_plus_eta(i) .LE. hp_dry_tolerance ) .OR.                     &
            ( h_center(i) .LE. hp_dry_tolerance ) ) w_eta(i) = 1.0_wp

       h_minus(i) = ( 1.0_wp-w_eta(i) ) * h_minus_orig(i)                    &
            + w_eta(i) * h_minus_eta(i)
       h_plus(i) = ( 1.0_wp-w_eta(i) ) * h_plus_orig(i)                      &
            + w_eta(i) * h_plus_eta(i)

       IF ( ABS( B_plus(i)-B_minus(i) ) .GT. hp_bed_step_tolerance ) THEN
          momentum_weight = w_eta(i)
       ELSE
          momentum_weight = 0.0_wp
       END IF

       hu_minus(i) = ( 1.0_wp-momentum_weight ) * hu_minus_direct(i)         &
            + momentum_weight * h_minus_eta(i) * u_minus_candidate(i)
       hu_plus(i) = ( 1.0_wp-momentum_weight ) * hu_plus_direct(i)           &
            + momentum_weight * h_plus_eta(i) * u_plus_candidate(i)

       u_min = MINVAL( u_center(MAX(1,i-1):MIN(number_of_cells,i+1)) )
       u_max = MAXVAL( u_center(MAX(1,i-1):MIN(number_of_cells,i+1)) )

       dm_target = 0.5_wp * ( hu_plus(i)-hu_minus(i) )                       &
            - u_center(i) * 0.5_wp * ( h_plus(i)-h_minus(i) )
       dm_lower = MAX( ( u_center(i)-u_max ) * h_minus(i),                   &
            ( u_min-u_center(i) ) * h_plus(i) )
       dm_upper = MIN( ( u_center(i)-u_min ) * h_minus(i),                   &
            ( u_max-u_center(i) ) * h_plus(i) )
       dm = MIN( MAX( dm_target, dm_lower ), dm_upper )

       hu_minus(i) = u_center(i) * h_minus(i) - dm
       hu_plus(i) = u_center(i) * h_plus(i) + dm

    END DO

  END SUBROUTINE reconstruct_hp_line

END MODULE hp_reconstruction_2d
