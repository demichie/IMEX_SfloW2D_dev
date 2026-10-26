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
       flow_depth, temperature, velocity_x, velocity_y, speed, horizontal_speed, &
       speed_squared, mixture_density, inverse_density, reduced_gravity, &
       gravity_coefficient, solid_fractions, curvature_term, curvature_xx, curvature_xy, &
       curvature_yy, use_pore_pressure, excess_pore_pressure, moving, use_stochastic, &
       stochastic_value, rate )
    REAL(wp), INTENT(IN) :: flow_depth, temperature, velocity_x, velocity_y
    REAL(wp), INTENT(IN) :: speed, horizontal_speed, speed_squared
    REAL(wp), INTENT(IN) :: mixture_density, inverse_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, solid_fractions(:)
    REAL(wp), INTENT(IN) :: curvature_xx, curvature_xy, curvature_yy
    REAL(wp), INTENT(IN) :: excess_pore_pressure, stochastic_value
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

  SUBROUTINE eval_rheology_newmodel_implicit(velocity_x,velocity_y, &
       mixture_density,rate_x,rate_y)
    COMPLEX(wp), INTENT(IN) :: velocity_x,velocity_y,mixture_density
    COMPLEX(wp), INTENT(OUT) :: rate_x,rate_y
    REAL(wp), PARAMETER :: drag_coefficient=0.1_wp
    rate_x=-drag_coefficient*mixture_density*velocity_x
    rate_y=-drag_coefficient*mixture_density*velocity_y
  END SUBROUTINE eval_rheology_newmodel_implicit

  SUBROUTINE solve_rheology_newmodel_analytic()
    CALL fatal_error('rheology_newmodel: analytic treatment unsupported')
  END SUBROUTINE solve_rheology_newmodel_analytic
END MODULE rheology_newmodel_2d
