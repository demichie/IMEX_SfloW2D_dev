!********************************************************************************
!> \brief Path-conservative central-upwind algebra
!>
!> This module contains the model-independent PCCU building blocks.  The
!> thermodynamic endpoint coefficient Gamma and the slope coefficient G are
!> supplied by the caller; this module integrates the slope-corrected
!> hydrostatic one-form and builds the oriented face pair.
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
  !> \brief Integrate the slope-corrected hydrostatic one-form between two endpoints.
  !>
  !> The returned vector is zero except in the momentum component normal to the
  !> selected direction.  In particular, the thermal component is identically
  !> zero by construction.  A face path passes the same shared G at both
  !> endpoints; a cell-internal path passes its two shared face values so that
  !> G varies linearly along the quadrature path.  The integrated one-form is
  !> -G*(Gamma*h*deta + 0.5*h**2*dGamma).
  !>
  !> \param[in] direction Cartesian path direction: PATH_DIR_X or PATH_DIR_Y.
  !> \param[in] h_left Depth at the negative/left path endpoint [m].
  !> \param[in] gamma_left Hydrostatic coefficient Gamma at the negative/left endpoint.
  !> \param[in] eta_left Free-surface elevation at the negative/left endpoint or stencil neighbour
  !>                     [m].
  !> \param[in] grav_coeff_left Large-slope factor G at the negative/left endpoint.
  !> \param[in] h_right Depth at the positive/right path endpoint [m].
  !> \param[in] gamma_right Hydrostatic coefficient Gamma at the positive/right endpoint.
  !> \param[in] eta_right Free-surface elevation at the positive/right endpoint or stencil neighbour
  !>                      [m].
  !> \param[in] grav_coeff_right Large-slope factor G at the positive/right endpoint.
  !> \param[out] path_contribution Integrated nonconservative path contribution, one entry per
  !>                               equation.
  !******************************************************************************

  SUBROUTINE eval_hydrostatic_path( direction, h_left, gamma_left, eta_left,   &
       grav_coeff_left, h_right, gamma_right, eta_right, grav_coeff_right,    &
       path_contribution )

    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: h_left, gamma_left, eta_left, grav_coeff_left
    REAL(wp), INTENT(IN) :: h_right, gamma_right, eta_right, grav_coeff_right
    REAL(wp), INTENT(OUT) :: path_contribution(:)

    REAL(wp) :: state_left(4), state_right(4)

    state_left = [ h_left, gamma_left, eta_left, grav_coeff_left ]
    state_right = [ h_right, gamma_right, eta_right, grav_coeff_right ]

    CALL eval_path_contribution( direction, state_left, state_right,          &
         hydrostatic_integrand, path_contribution, n_quad=3 )

  END SUBROUTINE eval_hydrostatic_path

  !******************************************************************************
  !> \brief Combine inertial endpoint fluxes with central-upwind dissipation.
  !>
  !> \param[in] a_minus Nonpositive lower characteristic bound at the face [m s^-1].
  !> \param[in] a_plus Nonnegative upper characteristic bound at the face [m s^-1].
  !> \param[in] flux_left Inertial conservative flux of the negative/left endpoint state.
  !> \param[in] flux_right Inertial conservative flux of the positive/right endpoint state.
  !> \param[in] state_left Conservative state at the negative/left endpoint.
  !> \param[in] state_right Conservative state at the positive/right endpoint.
  !> \param[out] numerical_flux Central-upwind conservative flux before the path split.
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
  !> \brief Split a numerical flux and path jump into the two cell-facing values.
  !>
  !> \param[in] numerical_flux Central-upwind conservative flux before the path split.
  !> \param[in] path_contribution Integrated nonconservative path contribution, one entry per
  !>                              equation.
  !> \param[in] a_minus Nonpositive lower characteristic bound at the face [m s^-1].
  !> \param[in] a_plus Nonnegative upper characteristic bound at the face [m s^-1].
  !> \param[out] value_left Oriented face value consumed by the negative/left adjacent cell.
  !> \param[out] value_right Oriented face value consumed by the positive/right adjacent cell.
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

  !> \brief Evaluate the hydrostatic momentum one-form at a quadrature state.
  !>
  !> \param[in] direction Cartesian path direction: PATH_DIR_X or PATH_DIR_Y.
  !> \param[in] path_state Extended state evaluated at the current path quadrature coordinate.
  !> \param[in] dstate_ds Derivative of the extended path state with respect to its unit-interval
  !>                      coordinate.
  !> \param[out] integrand Model one-form contracted with the path tangent, one entry per equation.

  SUBROUTINE hydrostatic_integrand( direction, path_state, dstate_ds,         &
       integrand )

    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: path_state(:), dstate_ds(:)
    REAL(wp), INTENT(OUT) :: integrand(:)

    REAL(wp) :: h, gamma, grav_coeff
    INTEGER :: normal_momentum

    IF ( SIZE(path_state) .NE. 4 .OR. SIZE(dstate_ds) .NE. 4 ) THEN
       ERROR STOP 'hydrostatic_integrand: expected state (h,Gamma,eta,G)'
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
    grav_coeff = path_state(4)
    integrand = 0.0_wp
    integrand(normal_momentum) = -grav_coeff * ( gamma*h*dstate_ds(3)        &
         + 0.5_wp*h*h*dstate_ds(2) )

  END SUBROUTINE hydrostatic_integrand

END MODULE pccu_2d
