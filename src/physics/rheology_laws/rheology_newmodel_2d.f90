!> \brief D17 extension fixture: simple linear Newton-only drag.
!> \details Fixed coefficient k=0.1 s^-1. This demonstration law is intentionally
!> independent of model-specific namelist parameters and numerical assemblers.
MODULE rheology_newmodel_2d
  USE parameters_2d, ONLY : wp
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_rheology_newmodel_explicit
  PUBLIC :: eval_rheology_newmodel_semi_implicit
  PUBLIC :: eval_rheology_newmodel_implicit
  PUBLIC :: solve_rheology_newmodel_analytic
CONTAINS
  SUBROUTINE eval_rheology_newmodel_explicit( &
       h, T, u, v, mod_vel, mod_hor_vel, &
       mod_vel2, rhom, inv_rhom, red_grav, &
       grav_coeff, alphas, curvature_term, B_second_xx, B_second_xy, &
       B_second_yy, use_pore_pressure, exc_pore_pres, moving, use_stochastic, &
       stochastic_value, rate )
    REAL(wp), INTENT(IN) :: h, T, u, v
    REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
    REAL(wp), INTENT(IN) :: rhom, inv_rhom, red_grav
    REAL(wp), INTENT(IN) :: grav_coeff, alphas(:)
    REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
    REAL(wp), INTENT(IN) :: exc_pore_pres, stochastic_value
    LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
    REAL(wp), INTENT(OUT) :: rate(2)
    ! Keep the common four-treatment interface, but this fixture is Newton-only.
    rate=0.0_wp
    CALL fatal_error('rheology_newmodel: explicit treatment unsupported')
  END SUBROUTINE eval_rheology_newmodel_explicit

  SUBROUTINE eval_rheology_newmodel_semi_implicit(rate_x,rate_y)
    REAL(wp), INTENT(OUT) :: rate_x,rate_y
    rate_x=0.0_wp
    rate_y=0.0_wp
  END SUBROUTINE eval_rheology_newmodel_semi_implicit

  SUBROUTINE eval_rheology_newmodel_implicit(u,v, &
       rhom,rate_x,rate_y)
    COMPLEX(wp), INTENT(IN) :: u,v,rhom
    COMPLEX(wp), INTENT(OUT) :: rate_x,rate_y
    REAL(wp), PARAMETER :: drag_coefficient=0.1_wp
    rate_x=-drag_coefficient*rhom*u
    rate_y=-drag_coefficient*rhom*v
  END SUBROUTINE eval_rheology_newmodel_implicit

  SUBROUTINE solve_rheology_newmodel_analytic()
    CALL fatal_error('rheology_newmodel: analytic treatment unsupported')
  END SUBROUTINE solve_rheology_newmodel_analytic
END MODULE rheology_newmodel_2d
