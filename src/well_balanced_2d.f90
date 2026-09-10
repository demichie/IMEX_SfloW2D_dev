!********************************************************************************
!> \brief Algebraic building blocks for the path-conservative topography scheme
!********************************************************************************
MODULE well_balanced_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_hydrostatic_path_integral

CONTAINS

  !******************************************************************************
  !> \brief Integrate the hydrostatic source along a linear projected path
  !>
  !> The path variables h, Gamma, G, B and the normal velocity are interpolated
  !> linearly between their left and right endpoints. Three-point
  !> Gauss-Legendre quadrature is exact for the resulting momentum and energy
  !> integrands, whose polynomial degrees are at most three and four.
  !******************************************************************************

  PURE SUBROUTINE eval_hydrostatic_path_integral( hL, hR, gammaL, gammaR,     &
       gravL, gravR, bedL, bedR, velL, velR, momentum_path, energy_path )

    REAL(wp), INTENT(IN) :: hL, hR
    REAL(wp), INTENT(IN) :: gammaL, gammaR
    REAL(wp), INTENT(IN) :: gravL, gravR
    REAL(wp), INTENT(IN) :: bedL, bedR
    REAL(wp), INTENT(IN) :: velL, velR
    REAL(wp), INTENT(OUT) :: momentum_path
    REAL(wp), INTENT(OUT) :: energy_path

    REAL(wp), PARAMETER :: gauss_abscissa(3) = [                           &
         0.5_wp * (1.0_wp - SQRT(3.0_wp / 5.0_wp)), 0.5_wp,                &
         0.5_wp * (1.0_wp + SQRT(3.0_wp / 5.0_wp)) ]
    REAL(wp), PARAMETER :: gauss_weight(3) = [                              &
         5.0_wp / 18.0_wp, 4.0_wp / 9.0_wp, 5.0_wp / 18.0_wp ]

    REAL(wp) :: delta_bed
    REAL(wp) :: delta_grav
    REAL(wp) :: h_path
    REAL(wp) :: gamma_path
    REAL(wp) :: grav_path
    REAL(wp) :: vel_path
    REAL(wp) :: momentum_integrand
    REAL(wp) :: s

    INTEGER :: i_quad

    delta_bed = bedR - bedL
    delta_grav = gravR - gravL

    momentum_path = 0.0_wp
    energy_path = 0.0_wp

    DO i_quad = 1, 3

       s = gauss_abscissa(i_quad)
       h_path = hL + s * (hR - hL)
       gamma_path = gammaL + s * (gammaR - gammaL)
       grav_path = gravL + s * delta_grav
       vel_path = velL + s * (velR - velL)

       momentum_integrand = -gamma_path * grav_path * h_path * delta_bed    &
            + 0.5_wp * gamma_path * h_path**2 * delta_grav

       momentum_path = momentum_path                                       &
            + gauss_weight(i_quad) * momentum_integrand
       energy_path = energy_path                                           &
            + gauss_weight(i_quad) * vel_path * momentum_integrand

    END DO

  END SUBROUTINE eval_hydrostatic_path_integral

END MODULE well_balanced_2d
