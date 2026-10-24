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
  SUBROUTINE eval_rheology_newmodel_explicit()
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
