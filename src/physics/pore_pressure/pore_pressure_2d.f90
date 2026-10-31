!> Dispatch pore-pressure models through the common four-treatment interface.
!> The model list contains canonical names, allowing future physical models
!> to be added without changing any caller. Darcy is the only model today.
MODULE pore_pressure_2d
  USE parameters_2d, ONLY : wp, n_solid
  USE diagnostics_2d, ONLY : fatal_error
  USE pore_pressure_darcy_2d, ONLY : &
       darcy_explicit => eval_pore_pressure_dissipation_explicit, &
       darcy_semi_implicit => eval_pore_pressure_dissipation_semi_implicit, &
       darcy_implicit => eval_pore_pressure_dissipation_implicit, &
       darcy_analytic => solve_pore_pressure_dissipation_analytic, &
       solve_frozen_darcy_mass_pressure, apply_frozen_darcy_withdrawal, &
       advance_pore_pressure_diffusion, eval_pore_pressure_gas_loss, &
       eval_pore_pressure_inlet_source, eval_pore_pressure_mass_exchange_source
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_pore_pressure_explicit, eval_pore_pressure_semi_implicit
  PUBLIC :: eval_pore_pressure_implicit, solve_pore_pressure_analytic
  ! Compatibility entry points while the full conservative source contract
  ! and the shared carrier-inventory limiters are migrated.
  PUBLIC :: eval_pore_pressure_dissipation_explicit
  PUBLIC :: eval_pore_pressure_dissipation_semi_implicit
  PUBLIC :: eval_pore_pressure_dissipation_implicit
  PUBLIC :: solve_pore_pressure_dissipation_analytic
  PUBLIC :: solve_frozen_darcy_mass_pressure, apply_frozen_darcy_withdrawal
  PUBLIC :: advance_pore_pressure_diffusion, eval_pore_pressure_gas_loss
  PUBLIC :: eval_pore_pressure_inlet_source, eval_pore_pressure_mass_exchange_source
CONTAINS

  SUBROUTINE eval_pore_pressure_explicit(models, flow_depth, temperature, &
       mixture_density, solid_fractions, excess_pressure, rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    REAL(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    REAL(wp), INTENT(IN) :: solid_fractions(n_solid), excess_pressure
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: model_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      model_rate=0.0_wp
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_explicit(flow_depth,temperature,mixture_density, &
             solid_fractions,excess_pressure,model_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown explicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_explicit

  SUBROUTINE eval_pore_pressure_semi_implicit(models,rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    REAL(wp), INTENT(OUT) :: rate
    REAL(wp) :: model_rate
    INTEGER :: i
    rate=0.0_wp
    DO i=1,SIZE(models)
      model_rate=0.0_wp
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_semi_implicit(model_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown semi-implicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_semi_implicit

  SUBROUTINE eval_pore_pressure_implicit(models, flow_depth, temperature, &
       mixture_density, solid_fractions, excess_pressure, rate)
    CHARACTER(*), INTENT(IN) :: models(:)
    COMPLEX(wp), INTENT(IN) :: flow_depth,temperature,mixture_density
    COMPLEX(wp), INTENT(IN) :: solid_fractions(n_solid),excess_pressure
    COMPLEX(wp), INTENT(OUT) :: rate
    COMPLEX(wp) :: model_rate
    INTEGER :: i
    rate=CMPLX(0.0_wp,0.0_wp,wp)
    DO i=1,SIZE(models)
      model_rate=CMPLX(0.0_wp,0.0_wp,wp)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_implicit(flow_depth,temperature,mixture_density, &
             solid_fractions,excess_pressure,model_rate)
      CASE DEFAULT
        CALL fatal_error('Unknown implicit pore-pressure model: '//TRIM(models(i)))
      END SELECT
      rate=rate+model_rate
    END DO
  END SUBROUTINE eval_pore_pressure_implicit

  SUBROUTINE solve_pore_pressure_analytic(models,excess_pressure,interval, &
       flow_depth,temperature,solid_fractions,pressure_after)
    CHARACTER(*), INTENT(IN) :: models(:)
    REAL(wp), INTENT(IN) :: excess_pressure,interval,flow_depth,temperature
    REAL(wp), INTENT(IN) :: solid_fractions(n_solid)
    REAL(wp), INTENT(OUT) :: pressure_after
    REAL(wp) :: current_pressure
    INTEGER :: i
    current_pressure=excess_pressure
    DO i=1,SIZE(models)
      SELECT CASE(TRIM(models(i)))
      CASE('pore_pressure_darcy')
        CALL darcy_analytic(current_pressure,interval,flow_depth, &
             temperature,solid_fractions,pressure_after)
      CASE DEFAULT
        CALL fatal_error('Unknown analytic pore-pressure model: '//TRIM(models(i)))
      END SELECT
      current_pressure=pressure_after
    END DO
    pressure_after=current_pressure
  END SUBROUTINE solve_pore_pressure_analytic

  ! Transitional adapters: existing tests and source assembly still call
  ! the old scalar-only interface until the full process rates are wired.
  SUBROUTINE eval_pore_pressure_dissipation_explicit(h,t,rho,alpha,p,rate)
    REAL(wp), INTENT(IN) :: h,t,rho,alpha(n_solid),p
    REAL(wp), INTENT(OUT) :: rate
    CALL eval_pore_pressure_explicit(['pore_pressure_darcy'],h,t,rho,alpha,p,rate)
  END SUBROUTINE eval_pore_pressure_dissipation_explicit

  SUBROUTINE eval_pore_pressure_dissipation_semi_implicit(rate)
    REAL(wp), INTENT(OUT) :: rate
    CALL eval_pore_pressure_semi_implicit(['pore_pressure_darcy'],rate)
  END SUBROUTINE eval_pore_pressure_dissipation_semi_implicit

  SUBROUTINE eval_pore_pressure_dissipation_implicit(h,t,rho,alpha,p,rate)
    COMPLEX(wp), INTENT(IN) :: h,t,rho,alpha(n_solid),p
    COMPLEX(wp), INTENT(OUT) :: rate
    CALL eval_pore_pressure_implicit(['pore_pressure_darcy'],h,t,rho,alpha,p,rate)
  END SUBROUTINE eval_pore_pressure_dissipation_implicit

  SUBROUTINE solve_pore_pressure_dissipation_analytic(p,dt,h,t,alpha,p_after)
    REAL(wp), INTENT(IN) :: p,dt,h,t,alpha(n_solid)
    REAL(wp), INTENT(OUT) :: p_after
    CALL solve_pore_pressure_analytic(['pore_pressure_darcy'],p,dt,h,t,alpha,p_after)
  END SUBROUTINE solve_pore_pressure_dissipation_analytic

END MODULE pore_pressure_2d
