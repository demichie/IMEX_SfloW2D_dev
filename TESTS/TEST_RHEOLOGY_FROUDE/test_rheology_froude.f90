!> \brief Regression for the domain of the Froude-dependent basal friction law.
PROGRAM test_rheology_froude

  USE parameters_2d
  USE constitutive_parameters_2d
  USE state_conversion_2d, ONLY : mixt_var, qp_to_qc
  USE equation_terms_2d, ONLY : eval_nh_semi_impl_terms, eval_implicit_terms
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  REAL(wp) :: qp(8), source(6), expected(6), q(6), implicit_source(6)
  COMPLEX(wp) :: complex_q(6), complex_source(6)
  REAL(wp) :: Ri, rho_m, rho_c, g_reduced, cp_c, cp_mix
  REAL(wp) :: speed, Fr, friction, slope_x, slope_y, G, z
  INTEGER :: mode

  CALL initialize_properties

  ! Set the ambient density from the same thermodynamic closure, guaranteeing
  ! exact neutral buoyancy independently of compiler rounding of ideal-gas rho.
  CALL make_state(1.0_wp, T_ambient, 0.0_wp, 3.0_wp, -4.0_wp)
  CALL mixt_var(qp, Ri, rho_m, rho_c, g_reduced, cp_c, cp_mix)
  rho_a_amb = rho_m

  DO mode = 0, 1
     stochastic_flag = mode == 1
     CALL make_state(1.0_wp, T_ambient, 0.0_wp, 3.0_wp, -4.0_wp)
     CALL assert_zero_source('moving neutral gas', 0.0_wp, 0.0_wp)

     ! Match the failing restart's tiny, but not machine-dry, ambient-gas tail.
     CALL make_state(3.3491779338939385E-12_wp, T_ambient, 0.0_wp,          &
          -3.8799233458989588E-10_wp, -3.7378768072356794E-6_wp)
     CALL assert_zero_source('neutral transported tail', 0.0_wp, 0.0_wp)

     CALL make_state(1.0_wp, 2.0_wp*T_ambient, 0.0_wp, 3.0_wp, -4.0_wp)
     CALL mixt_var(qp, Ri, rho_m, rho_c, g_reduced, cp_c, cp_mix)
     IF (g_reduced >= 0.0_wp) ERROR STOP 'buoyant fixture is not buoyant'
     CALL assert_zero_source('moving buoyant gas', 0.35_wp, -0.22_wp)

     CALL make_state(0.0_wp, T_ambient, 0.5_wp, 3.0_wp, -4.0_wp)
     CALL assert_zero_source('dry cell', 0.0_wp, 0.0_wp)

     CALL make_state(1.0_wp, T_ambient, 0.5_wp, 0.0_wp, 0.0_wp)
     CALL assert_zero_source('dense cell at rest', 0.0_wp, 0.0_wp)
  END DO

  ! In the positive-gravity domain preserve the existing formula, with and
  ! without slope correction and with deterministic/shifted/clipped Froude.
  CALL make_state(2.0_wp, T_ambient, 0.5_wp, 3.0_wp, -4.0_wp)
  CALL mixt_var(qp, Ri, rho_m, rho_c, g_reduced, cp_c, cp_mix)
  IF (g_reduced <= 0.0_wp) ERROR STOP 'dense fixture has no buoyancy load'
  DO mode = 0, 5
     slope_correction_flag = MOD(mode, 2) == 1
     stochastic_flag = mode >= 2
     slope_x = 0.35_wp
     slope_y = -0.22_wp
     G = 1.0_wp
     speed = 5.0_wp
     IF (slope_correction_flag) THEN
        G = 1.0_wp/(1.0_wp + slope_x**2 + slope_y**2)
        speed = SQRT(25.0_wp + (3.0_wp*slope_x - 4.0_wp*slope_y)**2)
     END IF
     z = 0.2_wp
     IF (mode >= 4) z = -100.0_wp
     Fr = speed/SQRT(g_reduced*qp(1))
     IF (stochastic_flag) Fr = MAX(0.0_wp, Fr + z)
     friction = rho_m*qp(1)*G*g_reduced                              &
          * (mu_inf + (mu_0 - mu_inf)*EXP(-Fr/Fr_0))
     expected = 0.0_wp
     expected(2:3) = -friction*[3.0_wp, -4.0_wp]/5.0_wp
     CALL eval_nh_semi_impl_terms(slope_x, slope_y, 0.0_wp, 0.0_wp,    &
          0.0_wp, G, qp, source, z)
     CALL assert_finite(source)
     IF (MAXVAL(ABS(source-expected)) >                             &
          128.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp, MAXVAL(ABS(expected)))) &
          ERROR STOP 'positive-gravity friction formula changed'
     IF (DOT_PRODUCT(source(2:3), qp(idx_u:idx_v)) >= 0.0_wp)        &
          ERROR STOP 'basal friction must dissipate momentum'
  END DO

  DEALLOCATE(rho_s, inv_rho_s, sp_heat_s, sp_heat_g, sp_gas_const_g)
  WRITE(*,*) 'PASS: Froude friction domain and dense-state formula verified'

CONTAINS

  !> \brief Configure one solid and one transported stochastic scalar.
  SUBROUTINE initialize_properties
    n_solid = 1
    n_add_gas = 0
    n_stoch_vars = 1
    n_pore_vars = 0
    n_vars = 6
    n_eqns = 6
    idx_h = 1
    idx_hu = 2
    idx_hv = 3
    idx_T = 4
    idx_solid_first = 5
    idx_solid_last = 5
    idx_add_gas_first = 6
    idx_add_gas_last = 5
    idx_stoch = 6
    idx_u = 7
    idx_v = 8
    ALLOCATE(rho_s(1), inv_rho_s(1), sp_heat_s(1))
    ALLOCATE(sp_heat_g(0), sp_gas_const_g(0))
    rho_s = 2000.0_wp
    inv_rho_s = 1.0_wp/rho_s
    sp_heat_s = 1617.0_wp
    sp_heat_a = 998.0_wp
    sp_gas_const_a = 287.051_wp
    pres = 101300.0_wp
    inv_pres = 1.0_wp/pres
    T_ambient = 300.0_wp
    rho_a_amb = pres/(sp_gas_const_a*T_ambient)
    grav = 9.81_wp
    eps_sing = 1.0E-8_wp
    eps_sing4 = eps_sing**4
    gas_flag = .TRUE.
    liquid_flag = .FALSE.
    rheology_flag = .TRUE.
    rheology_model = 9
    slope_correction_flag = .FALSE.
    curvature_term_flag = .FALSE.
    pore_pressure_flag = .FALSE.
    stoch_transport_flag = .TRUE.
    mu_0 = 0.1_wp
    mu_inf = 0.5_wp
    Fr_0 = 1.0_wp
  END SUBROUTINE initialize_properties

  !> \brief Populate the canonical primitive state and cached velocities.
  !> \param[in] height Flow thickness [m].
  !> \param[in] temperature Gas temperature [K].
  !> \param[in] solid_fraction Solid mass fraction.
  !> \param[in] u,v Horizontal velocity components [m s^-1].
  SUBROUTINE make_state(height, temperature, solid_fraction, u, v)
    REAL(wp), INTENT(IN) :: height, temperature, solid_fraction, u, v
    qp = 0.0_wp
    qp(1) = height
    qp(2) = height*u
    qp(3) = height*v
    qp(4) = temperature
    qp(5) = solid_fraction
    qp(idx_u) = u
    qp(idx_v) = v
  END SUBROUTINE make_state

  !> \brief Assert finite zero semi-implicit friction and both implicit API paths.
  !> \param[in] label Name printed on a failed assertion.
  !> \param[in] bx,by Bed slopes used by the friction evaluation.
  SUBROUTINE assert_zero_source(label, bx, by)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: bx, by
    CALL eval_nh_semi_impl_terms(bx, by, 0.0_wp, 0.0_wp, 0.0_wp,       &
         1.0_wp, qp, source, 0.2_wp)
    CALL assert_finite(source)
    IF (ANY(source /= 0.0_wp)) THEN
       WRITE(*,*) 'FAIL: ', label, source
       ERROR STOP 1
    END IF

    ! Model 9 is purely semi-implicit: the REAL and COMPLEX implicit paths
    ! must remain finite and zero, including a complex-step velocity probe.
    CALL qp_to_qc(qp, q)
    CALL eval_implicit_terms(bx, by, 0.2_wp, r_qj=q,               &
         r_nh_term_impl=implicit_source)
    CALL assert_finite(implicit_source)
    IF (ANY(implicit_source /= 0.0_wp)) ERROR STOP 'real implicit model 9 source'
    complex_q = CMPLX(q, 0.0_wp, wp)
    complex_q(2) = complex_q(2) + CMPLX(0.0_wp, 1.0E-20_wp, wp)
    CALL eval_implicit_terms(bx, by, 0.2_wp, c_qj=complex_q,       &
         c_nh_term_impl=complex_source)
    CALL assert_finite(REAL(complex_source,wp))
    CALL assert_finite(AIMAG(complex_source))
    IF (ANY(complex_source /= CMPLX(0.0_wp,0.0_wp,wp)))           &
         ERROR STOP 'complex implicit model 9 source'
  END SUBROUTINE assert_zero_source

  !> \brief Reject NaN/Inf, including optimized builds where comparisons can hide them.
  !> \param[in] values Source vector to inspect.
  SUBROUTINE assert_finite(values)
    REAL(wp), INTENT(IN) :: values(:)
    IF (.NOT. ALL(ieee_is_finite(values))) ERROR STOP 'nonfinite friction source'
  END SUBROUTINE assert_finite

END PROGRAM test_rheology_froude
