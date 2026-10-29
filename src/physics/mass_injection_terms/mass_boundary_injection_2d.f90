!> \brief Physical treatment interfaces for boundary_injection.
MODULE mass_boundary_injection_2d
    USE parameters_2d, ONLY : idx_T, idx_h, idx_hu, idx_hv
  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var
  USE thermal_laws_2d, ONLY : eval_inlet_thermal_energy_source
  USE pore_pressure_laws_2d, ONLY : eval_pore_pressure_inlet_source
  USE geometry_2d, ONLY : pi_g, bottom_source_area_ratio
  USE mass_injection_coefficients_2d
  USE diagnostics_2d, ONLY : fatal_error
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_boundary_injection_explicit
  PUBLIC :: eval_boundary_injection_semi_implicit
  PUBLIC :: eval_boundary_injection_implicit
  PUBLIC :: solve_boundary_injection_analytic
CONTAINS
  SUBROUTINE eval_boundary_injection_explicit( time, vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , xs_source,    &
         xg_source, time_param

    USE geometry_2d, ONLY : pi_g

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff

    source_bdry = 0.0_wp
    source_bdry(idx_T) = T_source

    IF ( time .GE. time_param(4) ) THEN

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    t_coeff = 0.0_wp

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_wp ) ) )

       END IF

    END IF

    ! Preserve the common dry-state convention while the periodic source is
    ! inactive. In particular, no composition is attached to zero thickness.
    IF (t_coeff .LE. EPSILON(1.0_wp)) RETURN

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(idx_h) = t_coeff * h_source
    source_bdry(idx_hu) = t_coeff**1.5_wp * h_source * vel_source * vect_x
    source_bdry(idx_hv) = t_coeff**1.5_wp * h_source * vel_source * vect_y
    source_bdry(idx_T) = T_source

    source_bdry(idx_solid_first:idx_solid_last) = xs_source(1:n_solid)
    source_bdry(idx_add_gas_first:idx_add_gas_last) = xg_source(1:n_add_gas)

    source_bdry(idx_u) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(idx_v) = t_coeff**0.5_wp * vel_source * vect_y

    RETURN

  END SUBROUTINE eval_boundary_injection_explicit

  SUBROUTINE eval_boundary_injection_semi_implicit(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: semi_implicit treatment is unsupported')
  END SUBROUTINE eval_boundary_injection_semi_implicit

  SUBROUTINE eval_boundary_injection_implicit(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: implicit treatment is unsupported')
  END SUBROUTINE eval_boundary_injection_implicit

  SUBROUTINE solve_boundary_injection_analytic(time, vect_x, vect_y, source_bdry)
    REAL(wp), INTENT(IN) :: time, vect_x, vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)
    CALL fatal_error('boundary_injection: analytic treatment is unsupported')
  END SUBROUTINE solve_boundary_injection_analytic
END MODULE mass_boundary_injection_2d
