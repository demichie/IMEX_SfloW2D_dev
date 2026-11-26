PROGRAM test_stochastic_evolution

  USE constitutive_parameters_2d, ONLY : T_ambient
  USE geometry_2d, ONLY : cell_size, comp_cells_x, comp_cells_y
  USE model_layout_2d, ONLY : model_layout_type
  USE parameters_2d, ONLY : wp, idx_stoch,             &
       n_add_gas, n_dyn_pres_levels, n_eqns, n_pore_vars, n_solid,         &
       n_stoch_vars, n_thickness_levels, n_vars, configure_canonical_state_indices, &
       idx_u, idx_v
  USE rheology_parameters_2d, ONLY : rheology_l1,rheology_l2
   USE state_fields, ONLY : state_type
  USE stochastic_module, ONLY : stochastic_workspace_type,stochastic_config_type, &
       stochastic_l1,stochastic_l2,convolve_2d,expFormNoise,getSigmaNoise

  IMPLICIT NONE

  TYPE(state_type) :: state,dual_state
  TYPE(stochastic_workspace_type) :: stochastic,upper,upper_reference
  TYPE(model_layout_type) :: model_layout
  REAL(wp), PARAMETER :: tolerance = 64.0_wp * EPSILON(1.0_wp)
  REAL(wp) :: expected_sigma
  REAL(wp) :: rectangular_field(3,2), identity_kernel(1,1)
  REAL(wp) :: convolution_result(3,2)
  REAL(wp) :: lower_before(3,2),upper_noise(3,2)
  INTEGER :: seed_size
  INTEGER, ALLOCATABLE :: seed_before(:),seed_after(:)

  comp_cells_x = 3
  comp_cells_y = 2
  n_solid = 0
  n_add_gas = 0
  n_stoch_vars = 1
  n_pore_vars = 0
  n_eqns = 5
  n_vars = 5
  CALL configure_canonical_state_indices()
  n_thickness_levels = 1
  n_dyn_pres_levels = 1
  T_ambient = 273.15_wp

  rheology_l1%stochastic_rheology_flag = .TRUE.
  stochastic_l1%stochastic_transport_flag = .TRUE.
  rheology_l1%rheology_model='none'
  stochastic_l1%length_spatial_corr = 0.0_wp
  stochastic_l1%tau_stochastic = 2.0_wp
  stochastic_l1%std_max = 0.0_wp
  stochastic_l1%sym_noise = 1.0_wp
  stochastic_l1%noise_pow_val = 2.0_wp

  CALL model_layout%initialize(1, n_vars, n_eqns, n_solid, n_add_gas,      &
       n_stoch_vars, n_pore_vars)
  CALL state%initialize(model_layout)
  CALL stochastic%initialize(stochastic_l1,rheology_l1)

  state%q(:,1,:,:) = 0.0_wp
  state%qp(:,1,:,:) = 0.0_wp
  state%q(1,1,1,1) = 2.0_wp
  state%q(idx_stoch,1,1,1) = 6.0_wp
  stochastic%Z = 4.0_wp

  CALL stochastic%prepare_timestep(state, 0.5_wp, .TRUE.)

  CALL assert_close('wet cell starts from transported hZ/h',               &
       stochastic%Z(1,1), 2.25_wp)
  CALL assert_close('dry cell retains and evolves local OU state',         &
       stochastic%Z(2,1), 3.0_wp)
  CALL assert_close('transport variable receives updated raw OU state',    &
       state%q(idx_stoch,1,1,1), 4.5_wp)
  CALL assert_close('dry transport variable remains zero',                 &
       state%q(idx_stoch,1,2,1), 0.0_wp)
  CALL assert_close('nonlinear map is applied only to friction field',     &
       stochastic%effective_Z(1,1), 2.25_wp**2)

  CALL stochastic%prepare_timestep(state, 0.5_wp, .TRUE.)
  CALL assert_close('effective field is not fed back into OU state',       &
       stochastic%Z(1,1), 1.6875_wp)

  ! The common optional slot exists even for a layer without transport.
  state%q(idx_stoch,1,1,1)=8.0_wp
  state%qp(idx_stoch,1,1,1)=4.0_wp
  stochastic_l1%stochastic_transport_flag=.TRUE.
  CALL stochastic%prepare_timestep(state,0.5_wp,.FALSE.)
  CALL assert_close('disabled layer clears conservative stochastic slot', &
       state%q(idx_stoch,1,1,1),0.0_wp)
  CALL assert_close('disabled layer clears primitive stochastic slot', &
       state%qp(idx_stoch,1,1,1),0.0_wp)

  ! Enabled transport follows the argument, not the legacy global flag.
  stochastic_l1%stochastic_transport_flag=.FALSE.
  state%q(idx_stoch,1,1,1)=6.0_wp
  stochastic%Z=4.0_wp
  CALL stochastic%prepare_timestep(state,0.5_wp,.TRUE.)
  CALL assert_close('selected transport ignores legacy global flag', &
       state%q(idx_stoch,1,1,1),4.5_wp)
  stochastic_l1%stochastic_transport_flag=.TRUE.

  stochastic_l1%std_min = 2.0_wp
  stochastic_l1%std_max = 10.0_wp
  stochastic_l1%std_slope_factor = 25.0_wp
  CALL assert_close('sigma at zero velocity', expFormNoise(0.0_wp,stochastic_l1),        &
       stochastic_l1%std_min)
  expected_sigma = stochastic_l1%std_max + (stochastic_l1%std_min - stochastic_l1%std_max) * EXP(-1.0_wp)
  CALL assert_close('sigma at activation velocity', expFormNoise(25.0_wp,stochastic_l1), &
       expected_sigma)

  state%q(1,1,1,1) = 1.0_wp
  state%qp(idx_u,1,1,1) = 3.0_wp
  state%qp(idx_v,1,1,1) = 4.0_wp
  rheology_l1%rheology_model='froude'
  CALL assert_close('sigma uses squared physical speed',                  &
       getSigmaNoise(state,1,1,stochastic_l1,rheology_l1), expected_sigma)

  state%q(1,1,2,1) = 0.0_wp
  CALL assert_close('dry cells use sigma_max', getSigmaNoise(state,2,1,stochastic_l1,rheology_l1),   &
       stochastic_l1%std_max)

  cell_size = 1.0_wp
  stochastic_l1%length_spatial_corr = 1.0_wp
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

  ! Independent layer coefficients are not aliased by the live L1 workspace.
  stochastic_l2=stochastic_config_type()
  rheology_l2%rheology_model='none'
  stochastic_l2%tau_stochastic=4.0_wp
  stochastic_l2%length_spatial_corr=3.0_wp
  CALL assert_close('L2 OU time is independent of L1',stochastic_l1%tau_stochastic,2.0_wp)
  CALL assert_close('L2 OU spatial scale is independent of L1', &
       stochastic_l1%length_spatial_corr,1.0_wp)
  CALL assert_close('independent L2 OU time',stochastic_l2%tau_stochastic,4.0_wp)

  ! Transport and stochastic rheology are independent properties of each layer.
  stochastic_l1%stochastic_transport_flag=.TRUE.
  stochastic_l2%stochastic_transport_flag=.FALSE.
  rheology_l1%stochastic_rheology_flag=.FALSE.
  CALL assert_true('L1 stochastic transport without stochastic rheology', &
       stochastic_l1%stochastic_transport_flag .AND. .NOT.rheology_l1%stochastic_rheology_flag)
  stochastic_l2%stochastic_transport_flag=.TRUE.
  CALL assert_true('L2 transport does not change L1 transport',stochastic_l1%stochastic_transport_flag)
  stochastic_l1%stochastic_transport_flag=.FALSE.
  rheology_l2%stochastic_rheology_flag=.TRUE.
  CALL assert_true('transport enabled only for L2', &
       .NOT.stochastic_l1%stochastic_transport_flag .AND. stochastic_l2%stochastic_transport_flag)
  CALL assert_true('stochastic rheology enabled only for L2', &
       .NOT.rheology_l1%stochastic_rheology_flag .AND. rheology_l2%stochastic_rheology_flag)
  CALL assert_true('L2 and L1 OU configurations remain independent', &
       stochastic_l2%length_spatial_corr/=stochastic_l1%length_spatial_corr)
  CALL assert_true('workspace references canonical L1 OU config',ASSOCIATED(stochastic%config,stochastic_l1))

  ! Exercise independent L2 OU evolution and conservative transport without
  ! passing through the production-only model-layout gate.
  ALLOCATE(dual_state%q(n_vars,2,comp_cells_x,comp_cells_y))
  ALLOCATE(dual_state%qp(n_vars,2,comp_cells_x,comp_cells_y))
  dual_state%q=0.0_wp
  dual_state%qp=0.0_wp
  dual_state%q(1,2,1,1)=2.0_wp
  dual_state%q(idx_stoch,2,1,1)=6.0_wp
  dual_state%q(idx_stoch,1,1,1)=11.0_wp
  stochastic_l2%length_spatial_corr=0.0_wp
  stochastic_l2%std_max=0.0_wp
  stochastic_l2%stochastic_seed=246813
  CALL random_seed(SIZE=seed_size)
  ALLOCATE(seed_before(seed_size),seed_after(seed_size))
  CALL random_seed(GET=seed_before)
  CALL upper%initialize(stochastic_l2,rheology_l2,layer=2)
  CALL random_seed(GET=seed_after)
  CALL assert_true('L2 initialization preserves the global L1 random stream',ALL(seed_before==seed_after))
  lower_before=dual_state%q(idx_stoch,1,:,:)
  CALL upper%prepare_timestep(dual_state,0.5_wp,.TRUE.)
  CALL assert_close('L2 OU uses its own correlation time',upper%Z(1,1),2.625_wp)
  CALL assert_close('L2 conservative transport uses L2 mass',dual_state%q(idx_stoch,2,1,1),5.25_wp)
  CALL assert_true('L2 OU never changes L1 transported state', &
       ALL(dual_state%q(idx_stoch,1,:,:)==lower_before))
  CALL assert_close('L2 sigma uses L2 mass and velocity', &
       getSigmaNoise(dual_state,1,1,stochastic_l2,rheology_l2,layer=2),0.0_wp)

  ! Equal seeds define equal independent streams even if an L1 update occurs
  ! between the two L2 updates. Neither stream depends on OpenMP scheduling.
  stochastic_l2%std_max=1.0_wp
  CALL upper%initialize(stochastic_l2,rheology_l2,layer=2)
  CALL upper_reference%initialize(stochastic_l2,rheology_l2,layer=2)
  CALL random_seed(GET=seed_before)
  CALL upper%update(dual_state,0.5_wp)
  CALL random_seed(GET=seed_after)
  CALL assert_true('L2 update preserves the global L1 random stream',ALL(seed_before==seed_after))
  upper_noise=upper%Z
  stochastic_l1%std_max=1.0_wp
  stochastic_l1%length_spatial_corr=0.0_wp
  stochastic_l1%stochastic_seed=111
  CALL stochastic%initialize(stochastic_l1,rheology_l1,layer=1)
  CALL stochastic%update(state,0.5_wp)
  CALL upper_reference%update(dual_state,0.5_wp)
  CALL assert_true('Independent L2 random stream survives intervening L1 update', &
       ALL(upper_reference%Z==upper_noise))
  DEALLOCATE(seed_before,seed_after)
  CALL upper_reference%finalize
  CALL upper%finalize
  DEALLOCATE(dual_state%q,dual_state%qp)

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
