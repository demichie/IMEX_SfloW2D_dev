PROGRAM test_stochastic_evolution

  USE constitutive_parameters_2d, ONLY : T_ambient
  USE geometry_2d, ONLY : cell_size, comp_cells_x, comp_cells_y
  USE model_layout_2d, ONLY : model_layout_type
  USE parameters_2d, ONLY : wp, idx_stoch, length_spatial_corr,             &
       n_add_gas, n_dyn_pres_levels, n_eqns, n_pore_vars, n_solid,         &
       n_stoch_vars, n_thickness_levels, n_vars, output_stoch_vars_flag,   &
       rheology_model, stochastic_flag, stoch_transport_flag, idx_u, idx_v
  USE state_2d, ONLY : state_type
  USE stochastic_module, ONLY : stochastic_workspace_type, sym_noise,     &
       noise_pow_val, std_max, std_min, std_slope_factor, tau_stochastic, &
       convolve_2d, expFormNoise, getSigmaNoise

  IMPLICIT NONE

  TYPE(state_type) :: state
  TYPE(stochastic_workspace_type) :: stochastic
  TYPE(model_layout_type) :: model_layout
  REAL(wp), PARAMETER :: tolerance = 64.0_wp * EPSILON(1.0_wp)
  REAL(wp) :: expected_sigma
  REAL(wp) :: rectangular_field(3,2), identity_kernel(1,1)
  REAL(wp) :: convolution_result(3,2)

  comp_cells_x = 3
  comp_cells_y = 2
  n_solid = 0
  n_add_gas = 0
  n_stoch_vars = 1
  n_pore_vars = 0
  n_eqns = 5
  n_vars = 5
  idx_stoch = 5
  idx_u = n_vars + 1
  idx_v = n_vars + 2
  n_thickness_levels = 1
  n_dyn_pres_levels = 1
  T_ambient = 273.15_wp

  stochastic_flag = .TRUE.
  stoch_transport_flag = .TRUE.
  output_stoch_vars_flag = .FALSE.
  rheology_model = 0
  length_spatial_corr = 0.0_wp
  tau_stochastic = 2.0_wp
  std_max = 0.0_wp
  sym_noise = 1.0_wp
  noise_pow_val = 2.0_wp

  CALL model_layout%initialize(1, n_vars, n_eqns, n_solid, n_add_gas,      &
       n_stoch_vars, n_pore_vars, .FALSE.)
  CALL state%initialize(model_layout)
  CALL stochastic%initialize

  state%q = 0.0_wp
  state%qp = 0.0_wp
  state%q(1,1,1) = 2.0_wp
  state%q(idx_stoch,1,1) = 6.0_wp
  stochastic%Z = 4.0_wp

  CALL stochastic%prepare_timestep(state, 0.5_wp)

  CALL assert_close('wet cell starts from transported hZ/h',               &
       stochastic%Z(1,1), 2.25_wp)
  CALL assert_close('dry cell retains and evolves local OU state',         &
       stochastic%Z(2,1), 3.0_wp)
  CALL assert_close('transport variable receives updated raw OU state',    &
       state%q(idx_stoch,1,1), 4.5_wp)
  CALL assert_close('dry transport variable remains zero',                 &
       state%q(idx_stoch,2,1), 0.0_wp)
  CALL assert_close('nonlinear map is applied only to friction field',     &
       stochastic%effective_Z(1,1), 2.25_wp**2)

  CALL stochastic%prepare_timestep(state, 0.5_wp)
  CALL assert_close('effective field is not fed back into OU state',       &
       stochastic%Z(1,1), 1.6875_wp)

  std_min = 2.0_wp
  std_max = 10.0_wp
  std_slope_factor = 25.0_wp
  CALL assert_close('sigma at zero velocity', expFormNoise(0.0_wp),        &
       std_min)
  expected_sigma = std_max + (std_min - std_max) * EXP(-1.0_wp)
  CALL assert_close('sigma at activation velocity', expFormNoise(25.0_wp), &
       expected_sigma)

  state%q(1,1,1) = 1.0_wp
  state%qp(idx_u,1,1) = 3.0_wp
  state%qp(idx_v,1,1) = 4.0_wp
  rheology_model = 9
  CALL assert_close('sigma uses squared physical speed',                  &
       getSigmaNoise(state,1,1), expected_sigma)

  state%q(1,2,1) = 0.0_wp
  CALL assert_close('dry cells use sigma_max', getSigmaNoise(state,2,1),   &
       std_max)

  cell_size = 1.0_wp
  length_spatial_corr = 1.0_wp
  CALL stochastic%generate_kernel
  CALL assert_true('kernel size is odd',                                   &
       MOD(SIZE(stochastic%conv_kernel,1), 2) == 1 .AND.                    &
       MOD(SIZE(stochastic%conv_kernel,2), 2) == 1)
  CALL assert_close('kernel is normalized', SUM(stochastic%conv_kernel),   &
       1.0_wp)
  CALL assert_close('kernel bandwidth is half the correlation length',     &
       stochastic%conv_kernel(4,3) / stochastic%conv_kernel(3,3),          &
       EXP(-2.0_wp))
  CALL assert_true('kernel is symmetric',                                  &
       MAXVAL(ABS(stochastic%conv_kernel -                                 &
       stochastic%conv_kernel(SIZE(stochastic%conv_kernel,1):1:-1,        &
       SIZE(stochastic%conv_kernel,2):1:-1))) <= tolerance)

  rectangular_field = RESHAPE([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp,   &
       6.0_wp], SHAPE(rectangular_field))
  identity_kernel = 1.0_wp
  CALL convolve_2d(rectangular_field, identity_kernel, convolution_result)
  CALL assert_true('convolution preserves rectangular (x,y) layout',      &
       MAXVAL(ABS(convolution_result - rectangular_field)) <= tolerance)

  CALL stochastic%finalize
  CALL state%finalize
  CALL model_layout%finalize

  WRITE(*,*) 'PASS: stochastic OU evolution and transport coupling verified'

CONTAINS

  SUBROUTINE assert_close(label, actual, expected)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual
    REAL(wp), INTENT(IN) :: expected

    IF (ABS(actual - expected) > tolerance * MAX(1.0_wp, ABS(expected))) THEN
       WRITE(*,*) 'FAIL: ', label, actual, expected
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_close

  SUBROUTINE assert_true(label, condition)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF (.NOT. condition) THEN
       WRITE(*,*) 'FAIL: ', label
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

END PROGRAM test_stochastic_evolution
