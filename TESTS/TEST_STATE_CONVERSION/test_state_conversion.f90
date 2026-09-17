PROGRAM test_state_conversion

  USE parameters_2d
  USE constitutive_2d

  IMPLICIT NONE

  CALL initialize_test_properties

  CALL run_wet_case('gas-solid alpha energy', .FALSE., .TRUE., .TRUE.)
  CALL run_wet_case('gas-solid h-alpha temperature', .FALSE., .FALSE., .FALSE.)
  CALL run_wet_case('gas-liquid-solid alpha energy', .TRUE., .TRUE., .TRUE.)
  CALL run_wet_case('gas-liquid-solid h-alpha temperature', .TRUE., .FALSE., .FALSE.)
  CALL run_zero_carrier_case
  CALL run_dry_cases
  CALL run_complex_step_check

  CALL finalize_test_properties

  WRITE (*, *) 'PASS: state conversion verified'

CONTAINS

  SUBROUTINE initialize_test_properties

    n_solid = 2
    n_add_gas = 1
    n_stoch_vars = 1
    n_pore_vars = 1

    ALLOCATE (rho_s(n_solid), inv_rho_s(n_solid))
    ALLOCATE (diam_s(n_solid), sphericity_s(n_solid), sp_heat_s(n_solid))
    ALLOCATE (sp_heat_g(n_add_gas), sp_gas_const_g(n_add_gas))

    rho_s = [2500.0_wp, 1800.0_wp]
    inv_rho_s = 1.0_wp/rho_s
    diam_s = [1.0E-4_wp, 2.0E-4_wp]
    sphericity_s = 1.0_wp
    sp_heat_s = [900.0_wp, 1100.0_wp]

    sp_heat_a = 1005.0_wp
    sp_gas_const_a = 287.05_wp
    sp_heat_g = [1850.0_wp]
    sp_gas_const_g = [461.5_wp]
    sp_heat_l = 4180.0_wp
    rho_l = 1000.0_wp
    inv_rho_l = 1.0_wp/rho_l

    pres = 101325.0_wp
    inv_pres = 1.0_wp/pres
    T_ambient = 300.0_wp
    rho_a_amb = pres/(sp_gas_const_a*T_ambient)
    grav = 9.81_wp

    eps_sing = 1.0E-8_wp
    eps_sing4 = eps_sing**4
    sutherland_flag = .FALSE.
    rheology_flag = .FALSE.
    slope_correction_flag = .FALSE.
    curvature_term_flag = .FALSE.
    stochastic_flag = .TRUE.
    stoch_transport_flag = .TRUE.
    pore_pressure_flag = .TRUE.
    gas_flag = .TRUE.

  END SUBROUTINE initialize_test_properties

  SUBROUTINE configure_layout(has_liquid, has_energy, stores_alpha)

    LOGICAL, INTENT(IN) :: has_liquid, has_energy, stores_alpha

    liquid_flag = has_liquid
    energy_flag = has_energy
    alpha_flag = stores_alpha

    idx_h = 1
    idx_hu = 2
    idx_hv = 3
    idx_T = 4
    idx_alfas_first = 5
    idx_alfas_last = 4 + n_solid
    idx_addGas_first = idx_alfas_last + 1
    idx_addGas_last = idx_addGas_first + n_add_gas - 1
    idx_stoch = idx_addGas_last + 1
    idx_pore = idx_stoch + 1

    n_vars = 4 + n_solid + n_add_gas + n_stoch_vars + n_pore_vars
    IF (liquid_flag) n_vars = n_vars + 1
    n_eqns = n_vars
    idx_u = n_vars + 1
    idx_v = n_vars + 2

  END SUBROUTINE configure_layout

  SUBROUTINE make_wet_state(qp)

    REAL(wp), INTENT(OUT) :: qp(n_vars + 2)

    REAL(wp), PARAMETER :: h = 2.0_wp
    REAL(wp), PARAMETER :: u = 1.5_wp
    REAL(wp), PARAMETER :: v = -0.25_wp
    REAL(wp) :: alphas(n_solid), alphag(n_add_gas)

    alphas = [0.10_wp, 0.05_wp]
    alphag = [0.08_wp]

    qp = 0.0_wp
    qp(1) = h
    qp(2) = h*u
    qp(3) = h*v
    qp(4) = 450.0_wp

    IF (alpha_flag) THEN
      qp(idx_alfas_first:idx_alfas_last) = alphas
      qp(idx_addGas_first:idx_addGas_last) = alphag
      IF (liquid_flag) qp(n_vars) = 0.20_wp
    ELSE
      qp(idx_alfas_first:idx_alfas_last) = h*alphas
      qp(idx_addGas_first:idx_addGas_last) = h*alphag
      IF (liquid_flag) qp(n_vars) = h*0.20_wp
    END IF

    qp(idx_stoch) = 0.37_wp
    qp(idx_pore) = 1250.0_wp
    qp(idx_u) = u
    qp(idx_v) = v

  END SUBROUTINE make_wet_state

  SUBROUTINE run_wet_case(label, has_liquid, has_energy, stores_alpha)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: has_liquid, has_energy, stores_alpha

    REAL(wp), ALLOCATABLE :: qp0(:), qp1(:), q0(:), q1(:)
    REAL(wp) :: p_dyn

    CALL configure_layout(has_liquid, has_energy, stores_alpha)
    ALLOCATE (qp0(n_vars + 2), qp1(n_vars + 2), q0(n_vars), q1(n_vars))

    CALL make_wet_state(qp0)
    CALL qp_to_qc(qp0, q0)
    CALL qc_to_qp(q0, qp1, p_dyn)
    CALL qp_to_qc(qp1, q1)

    CALL assert_close_vector(TRIM(label)//' q-qp-q', q1, q0, 2.0E-12_wp)
    CALL assert_close_vector(TRIM(label)//' qp', qp1, qp0, 2.0E-12_wp)
    CALL assert_true(TRIM(label)//' dynamic pressure finite', &
                     ieee_is_finite(p_dyn) .AND. p_dyn .GE. 0.0_wp)
    CALL check_real_complex(TRIM(label), q0, 2.0E-12_wp)

    DEALLOCATE (qp0, qp1, q0, q1)

  END SUBROUTINE run_wet_case

  SUBROUTINE run_zero_carrier_case

    REAL(wp), ALLOCATABLE :: q(:), real_outputs(:)
    COMPLEX(wp), ALLOCATABLE :: cq(:), complex_outputs(:)
    REAL(wp), PARAMETER :: test_temperature = 400.0_wp

    CALL configure_layout(.FALSE., .TRUE., .TRUE.)
    ALLOCATE (q(n_vars), real_outputs(5), cq(n_vars), complex_outputs(5))

    q = 0.0_wp
    q(1) = rho_s(1)
    q(2) = 0.2_wp*q(1)
    q(3) = -0.1_wp*q(1)
    q(4) = q(1)*sp_heat_s(1)*test_temperature +                         &
           0.5_wp*(q(2)**2 + q(3)**2)/q(1)
    q(idx_alfas_first) = q(1)

    CALL real_conversion_outputs(q, real_outputs)
    cq = CMPLX(q, 0.0_wp, wp)
    CALL complex_conversion_outputs(cq, complex_outputs)

    CALL assert_true('zero carrier real conversion finite',                    &
                     ALL(ieee_is_finite(real_outputs)))
    CALL assert_true('zero carrier complex conversion finite',                 &
                     ALL(ieee_is_finite(REAL(complex_outputs))) .AND.          &
                     ALL(ieee_is_finite(AIMAG(complex_outputs))))
    CALL assert_close_vector('zero carrier real-complex agreement',            &
                             REAL(complex_outputs), real_outputs, 2.0E-12_wp)

    DEALLOCATE (q, real_outputs, cq, complex_outputs)

  END SUBROUTINE run_zero_carrier_case

  SUBROUTINE run_dry_cases

    REAL(wp), ALLOCATABLE :: q0(:), q1(:), qp(:)
    REAL(wp) :: p_dyn

    CALL configure_layout(.FALSE., .TRUE., .TRUE.)
    ALLOCATE (q0(n_vars), q1(n_vars), qp(n_vars + 2))

    q0 = 0.0_wp
    CALL qc_to_qp(q0, qp, p_dyn)
    CALL qp_to_qc(qp, q1)
    CALL assert_close_vector('dry q-qp-q', q1, q0, 0.0_wp)
    CALL assert_true('dry physical defaults', qp(1) .EQ. 0.0_wp .AND. &
                     qp(idx_u) .EQ. 0.0_wp .AND. qp(idx_v) .EQ. 0.0_wp .AND. &
                     qp(4) .EQ. T_ambient)

    q0 = 0.0_wp
    q0(1) = 10.0_wp*EPSILON(1.0_wp)
    q0(2) = 0.25_wp*q0(1)
    q0(3) = -0.10_wp*q0(1)
    q0(4) = q0(1)*sp_heat_a*T_ambient
    CALL qc_to_qp(q0, qp, p_dyn)
    CALL assert_true('near-dry conversion finite', ALL(ieee_is_finite(qp)) .AND. &
                     ieee_is_finite(p_dyn) .AND. qp(1) .GE. 0.0_wp)

    DEALLOCATE (q0, q1, qp)

  END SUBROUTINE run_dry_cases

  SUBROUTINE check_real_complex(label, q, tolerance)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: q(n_vars), tolerance

    REAL(wp) :: h, u, v, T, rho_m, alphal, red_grav, p_dyn, Zs, pore
    REAL(wp) :: alphas(n_solid), alphag(n_add_gas)
    COMPLEX(wp) :: cq(n_vars), ch, cu, cv, cT, crho_m, cinv_rhom, cZs, cpore
    COMPLEX(wp) :: calphas(n_solid), calphag(n_add_gas)
    cq = CMPLX(q, 0.0_wp, wp)
    CALL r_phys_var(q, h, u, v, alphas, rho_m, T, alphal, alphag, red_grav, &
                    p_dyn, Zs, pore)
    CALL c_phys_var(cq, ch, cu, cv, cT, crho_m, calphas, calphag, cinv_rhom, &
                    cZs, cpore)

    CALL assert_close_scalar(TRIM(label)//' real-complex u', REAL(cu), u, tolerance)
    CALL assert_close_scalar(TRIM(label)//' real-complex v', REAL(cv), v, tolerance)
    CALL assert_close_scalar(TRIM(label)//' real-complex stochastic', REAL(cZs), &
                             Zs, tolerance)
    CALL assert_close_scalar(TRIM(label)//' real-complex pore', REAL(cpore), pore, &
                             tolerance)

    CALL assert_close_scalar(TRIM(label)//' real-complex h', REAL(ch), h, tolerance)
    CALL assert_close_scalar(TRIM(label)//' real-complex T', REAL(cT), T, tolerance)
    CALL assert_close_scalar(TRIM(label)//' real-complex rho', REAL(crho_m), &
                             rho_m, tolerance)
    CALL assert_close_vector(TRIM(label)//' real-complex solids', REAL(calphas), &
                             alphas, tolerance)
    CALL assert_close_vector(TRIM(label)//' real-complex gases', REAL(calphag), &
                             alphag, tolerance)
    CALL assert_true(TRIM(label)//' zero imaginary baseline', &
                     ABS(AIMAG(ch)) + ABS(AIMAG(cu)) + ABS(AIMAG(cv)) + &
                     ABS(AIMAG(cT)) + ABS(AIMAG(crho_m)) .EQ. 0.0_wp)

  END SUBROUTINE check_real_complex

  SUBROUTINE run_complex_step_check

    INTEGER :: component, output_idx
    INTEGER :: components(6)
    REAL(wp), PARAMETER :: complex_step = 1.0E-30_wp
    REAL(wp) :: finite_step
    REAL(wp), ALLOCATABLE :: qp(:), q(:), q_plus(:), q_minus(:)
    REAL(wp) :: f_plus(5), f_minus(5), fd_derivative(5), cs_derivative(5)
    REAL(wp) :: relative_difference
    COMPLEX(wp), ALLOCATABLE :: cq(:)
    COMPLEX(wp) :: cf(5)
    CHARACTER(LEN=80) :: derivative_label

    CALL configure_layout(.FALSE., .TRUE., .TRUE.)
    ALLOCATE (qp(n_vars + 2), q(n_vars), q_plus(n_vars), q_minus(n_vars), cq(n_vars))
    CALL make_wet_state(qp)
    CALL qp_to_qc(qp, q)

    components = [1, 2, 4, idx_alfas_first, idx_addGas_first, idx_pore]
    DO component = 1, SIZE(components)
      finite_step = 1.0E-6_wp*MAX(1.0_wp, ABS(q(components(component))))
      q_plus = q
      q_minus = q
      q_plus(components(component)) = q_plus(components(component)) + finite_step
      q_minus(components(component)) = q_minus(components(component)) - finite_step
      CALL real_conversion_outputs(q_plus, f_plus)
      CALL real_conversion_outputs(q_minus, f_minus)
      fd_derivative = (f_plus - f_minus)/(2.0_wp*finite_step)

      cq = CMPLX(q, 0.0_wp, wp)
      cq(components(component)) = cq(components(component)) + &
                                  CMPLX(0.0_wp, complex_step, wp)
      CALL complex_conversion_outputs(cq, cf)
      cs_derivative = AIMAG(cf)/complex_step

      DO output_idx = 1, SIZE(fd_derivative)
        WRITE (derivative_label, '(A,I0,A,I0)') 'complex-step component ', &
          components(component), ' output ', output_idx
        relative_difference = ABS(cs_derivative(output_idx) - &
                                  fd_derivative(output_idx)) / &
                              MAX(1.0_wp, ABS(fd_derivative(output_idx)))
        CALL assert_true(TRIM(derivative_label), &
                         relative_difference .LE. 2.0E-6_wp)
      END DO
    END DO

    DEALLOCATE (qp, q, q_plus, q_minus, cq)

  END SUBROUTINE run_complex_step_check

  SUBROUTINE real_conversion_outputs(q, outputs)

    REAL(wp), INTENT(IN) :: q(n_vars)
    REAL(wp), INTENT(OUT) :: outputs(5)
    REAL(wp) :: h, u, v, T, rho_m, alphal, red_grav, p_dyn, Zs, pore
    REAL(wp) :: alphas(n_solid), alphag(n_add_gas)

    CALL r_phys_var(q, h, u, v, alphas, rho_m, T, alphal, alphag, red_grav, &
                    p_dyn, Zs, pore)
    outputs = [h, u, v, T, rho_m]

  END SUBROUTINE real_conversion_outputs

  SUBROUTINE complex_conversion_outputs(q, outputs)

    COMPLEX(wp), INTENT(IN) :: q(n_vars)
    COMPLEX(wp), INTENT(OUT) :: outputs(5)
    COMPLEX(wp) :: h, u, v, T, rho_m, inv_rhom, Zs, pore
    COMPLEX(wp) :: alphas(n_solid), alphag(n_add_gas)

    CALL c_phys_var(q, h, u, v, T, rho_m, alphas, alphag, inv_rhom, Zs, pore)
    outputs = [h, u, v, T, rho_m]

  END SUBROUTINE complex_conversion_outputs

  SUBROUTINE assert_close_vector(label, actual, expected, tolerance)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual(:), expected(:), tolerance
    INTEGER :: max_index(1)
    REAL(wp) :: scale, max_error

    scale = MAX(1.0_wp, MAXVAL(ABS(expected)))
    max_error = MAXVAL(ABS(actual - expected))
    IF (SIZE(actual) .NE. SIZE(expected) .OR. max_error .GT. tolerance*scale) THEN
      max_index = MAXLOC(ABS(actual - expected))
      WRITE (*, *) '  max error/index/actual/expected:', max_error, max_index(1), &
                   actual(max_index(1)), expected(max_index(1))
      CALL assert_true(label, .FALSE.)
    END IF

  END SUBROUTINE assert_close_vector

  SUBROUTINE assert_close_scalar(label, actual, expected, tolerance)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual, expected, tolerance

    IF (ABS(actual - expected) .GT. tolerance*MAX(1.0_wp, ABS(expected))) THEN
      WRITE (*, *) '  error/actual/expected:', ABS(actual - expected), actual, expected
      CALL assert_true(label, .FALSE.)
    END IF

  END SUBROUTINE assert_close_scalar

  SUBROUTINE assert_true(label, condition)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF (.NOT. condition) THEN
      WRITE (*, *) 'FAIL: ', label
      ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

  ELEMENTAL LOGICAL FUNCTION ieee_is_finite(value)

    REAL(wp), INTENT(IN) :: value

    ieee_is_finite = (value .EQ. value) .AND. (ABS(value) .LE. HUGE(value))

  END FUNCTION ieee_is_finite

  SUBROUTINE finalize_test_properties

    DEALLOCATE (rho_s, inv_rho_s)
    DEALLOCATE (diam_s, sphericity_s, sp_heat_s)
    DEALLOCATE (sp_heat_g, sp_gas_const_g)

  END SUBROUTINE finalize_test_properties

END PROGRAM test_state_conversion
