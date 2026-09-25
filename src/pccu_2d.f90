!********************************************************************************
!> \brief Path-conservative central-upwind algebra
!>
!> This module contains the model-independent G=1 PCCU building blocks.  The
!> thermodynamic endpoint coefficient Gamma is supplied by equation_terms_2d;
!> this module only integrates the hydrostatic one-form and builds the oriented
!> face pair.
!********************************************************************************
MODULE pccu_2d

  USE parameters_2d, ONLY : wp
  USE nonconservative_2d, ONLY : eval_path_contribution
  USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER, PUBLIC :: pccu_speed_tolerance = 1.0E-14_wp

  PUBLIC :: eval_hydrostatic_path
  PUBLIC :: eval_central_upwind_flux
  PUBLIC :: eval_oriented_pccu_pair

CONTAINS

  !******************************************************************************
  !> \brief Integrate the G=1 hydrostatic path between two endpoint states
  !>
  !> The returned vector is zero except in the momentum component normal to the
  !> selected direction.  In particular, the thermal component is identically
  !> zero by construction.
  !******************************************************************************
  SUBROUTINE eval_hydrostatic_path( direction, h_left, gamma_left, eta_left,   &
       h_right, gamma_right, eta_right, path_contribution )

    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: h_left, gamma_left, eta_left
    REAL(wp), INTENT(IN) :: h_right, gamma_right, eta_right
    REAL(wp), INTENT(OUT) :: path_contribution(:)

    REAL(wp) :: state_left(3), state_right(3)

    state_left = [ h_left, gamma_left, eta_left ]
    state_right = [ h_right, gamma_right, eta_right ]

    CALL eval_path_contribution( direction, state_left, state_right,          &
         hydrostatic_integrand, path_contribution, n_quad=3 )

  END SUBROUTINE eval_hydrostatic_path

  !******************************************************************************
  !> \brief Central-upwind flux built from inertial endpoint fluxes
  !******************************************************************************
  SUBROUTINE eval_central_upwind_flux( a_minus, a_plus, flux_left, flux_right, &
       state_left, state_right, numerical_flux )

    REAL(wp), INTENT(IN) :: a_minus, a_plus
    REAL(wp), INTENT(IN) :: flux_left(:), flux_right(:)
    REAL(wp), INTENT(IN) :: state_left(:), state_right(:)
    REAL(wp), INTENT(OUT) :: numerical_flux(:)

    REAL(wp) :: denominator

    IF ( ( SIZE(flux_right) .NE. SIZE(flux_left) ) .OR.                       &
         ( SIZE(state_left) .NE. SIZE(flux_left) ) .OR.                       &
         ( SIZE(state_right) .NE. SIZE(flux_left) ) .OR.                      &
         ( SIZE(numerical_flux) .NE. SIZE(flux_left) ) ) THEN
       ERROR STOP 'eval_central_upwind_flux: inconsistent array sizes'
    END IF

    denominator = a_plus - a_minus
    IF ( denominator .GT. pccu_speed_tolerance ) THEN
       numerical_flux = ( a_plus*flux_left - a_minus*flux_right              &
            + a_plus*a_minus*(state_right-state_left) ) / denominator
    ELSE
       numerical_flux = 0.5_wp * ( flux_left + flux_right )
    END IF

  END SUBROUTINE eval_central_upwind_flux

  !******************************************************************************
  !> \brief Form the two oriented PCCU values seen by adjacent cells
  !******************************************************************************
  SUBROUTINE eval_oriented_pccu_pair( numerical_flux, path_contribution,       &
       a_minus, a_plus, value_left, value_right )

    REAL(wp), INTENT(IN) :: numerical_flux(:), path_contribution(:)
    REAL(wp), INTENT(IN) :: a_minus, a_plus
    REAL(wp), INTENT(OUT) :: value_left(:), value_right(:)

    REAL(wp) :: denominator

    IF ( ( SIZE(path_contribution) .NE. SIZE(numerical_flux) ) .OR.           &
         ( SIZE(value_left) .NE. SIZE(numerical_flux) ) .OR.                  &
         ( SIZE(value_right) .NE. SIZE(numerical_flux) ) ) THEN
       ERROR STOP 'eval_oriented_pccu_pair: inconsistent array sizes'
    END IF

    denominator = a_plus - a_minus
    IF ( denominator .GT. pccu_speed_tolerance ) THEN
       value_left = numerical_flux + a_minus/denominator * path_contribution
       value_right = numerical_flux + a_plus/denominator * path_contribution
    ELSE
       value_left = numerical_flux
       value_right = numerical_flux
    END IF

  END SUBROUTINE eval_oriented_pccu_pair

  SUBROUTINE hydrostatic_integrand( direction, path_state, dstate_ds,         &
       integrand )

    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)

    REAL(wp) :: h, gamma
    INTEGER :: normal_momentum

    IF ( SIZE(path_state) .NE. 3 .OR. SIZE(dstate_ds) .NE. 3 ) THEN
       ERROR STOP 'hydrostatic_integrand: expected state (h,Gamma,eta)'
    END IF

    SELECT CASE ( direction )
    CASE ( PATH_DIR_X )
       normal_momentum = 2
    CASE ( PATH_DIR_Y )
       normal_momentum = 3
    CASE DEFAULT
       ERROR STOP 'hydrostatic_integrand: invalid direction'
    END SELECT

    IF ( SIZE(integrand) .LT. normal_momentum ) THEN
       ERROR STOP 'hydrostatic_integrand: output has no normal momentum slot'
    END IF

    h = path_state(1)
    gamma = path_state(2)
    integrand = 0.0_wp
    integrand(normal_momentum) = -( gamma*h*dstate_ds(3)                     &
         + 0.5_wp*h*h*dstate_ds(2) )

  END SUBROUTINE hydrostatic_integrand

END MODULE pccu_2d
