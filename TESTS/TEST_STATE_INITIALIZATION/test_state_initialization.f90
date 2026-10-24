PROGRAM test_state_initialization

  USE constitutive_parameters_2d, ONLY : T_ambient
  USE domain_2d, ONLY : domain_type
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, comp_interfaces_x,     &
       comp_interfaces_y
  USE model_layout_2d, ONLY : model_layout_type
  USE parameters_2d, ONLY : wp, n_thickness_levels,                         &
       n_dyn_pres_levels
  USE state_2d, ONLY : state_type, correct_friction_momentum

  IMPLICIT NONE

  TYPE(state_type) :: state
  TYPE(domain_type) :: domain
  TYPE(model_layout_type) :: model_layout
  REAL(wp) :: before_momentum(2), after_momentum(2)

  comp_cells_x = 5
  comp_cells_y = 4
  comp_interfaces_x = comp_cells_x + 1
  comp_interfaces_y = comp_cells_y + 1

  n_thickness_levels = 2
  n_dyn_pres_levels = 3
  T_ambient = 273.15_wp

  CALL model_layout%initialize(1, 4, 4, 0, 0, 0, 0, .FALSE.)
  CALL state%initialize(model_layout)
  CALL domain%initialize

  CALL assert_true('conservative state', ALL(state%q .EQ. 0.0_wp))
  CALL assert_true('physical state',                                      &
       ALL(state%qp(1:3,:,:) .EQ. 0.0_wp) .AND.                           &
       ALL(state%qp(4,:,:) .EQ. T_ambient) .AND.                          &
       ALL(state%qp(5:6,:,:) .EQ. 0.0_wp))
  CALL assert_true('positive-thickness masks',                             &
       .NOT. ANY(state%hpos) .AND. .NOT. ANY(state%hpos_old))
  CALL assert_true('maximum diagnostics',                                 &
       ALL(state%hmax .EQ. 0.0_wp) .AND.                                  &
       ALL(state%pdynmax .EQ. 0.0_wp) .AND.                               &
       ALL(state%mod_vel_max .EQ. 0.0_wp))
  CALL assert_true('vulnerability diagnostics',                           &
       .NOT. ANY(state%vuln_table))

  CALL assert_true('domain counters', domain%solve_cells .EQ. 0 .AND.     &
       domain%solve_interfaces_x .EQ. 0 .AND.                             &
       domain%solve_interfaces_y .EQ. 0)
  CALL assert_true('release-time mask',                                   &
       ALL(domain%solve_mask_time .EQ. 0.0_wp))
  CALL assert_true('interior cell mask',                                  &
       .NOT. ANY(domain%solve_mask(2:comp_cells_x-1,                       &
       2:comp_cells_y-1)))
  CALL assert_true('boundary cell mask',                                  &
       ALL(domain%solve_mask(1,:)) .AND.                                  &
       ALL(domain%solve_mask(comp_cells_x,:)) .AND.                       &
       ALL(domain%solve_mask(:,1)) .AND.                                  &
       ALL(domain%solve_mask(:,comp_cells_y)))
  CALL assert_true('working masks', .NOT. ANY(domain%solve_mask_temp)     &
       .AND. .NOT. ANY(domain%solve_mask_x)                               &
       .AND. .NOT. ANY(domain%solve_mask_y))
  CALL assert_true('compact index storage', ALL(domain%j_cent .EQ. 0)     &
       .AND. ALL(domain%k_cent .EQ. 0)                                    &
       .AND. ALL(domain%j_stag_x .EQ. 0)                                  &
       .AND. ALL(domain%k_stag_x .EQ. 0)                                  &
       .AND. ALL(domain%j_stag_y .EQ. 0)                                  &
       .AND. ALL(domain%k_stag_y .EQ. 0))

  ! Isolated SI projection: preserve incoming direction and disallow growth.
  before_momentum=[3.0_wp,4.0_wp]
  after_momentum=[2.0_wp,1.0_wp]
  CALL correct_friction_momentum(before_momentum,after_momentum)
  CALL assert_true('post-SI projection onto prior direction', &
       MAXVAL(ABS(after_momentum-[1.2_wp,1.6_wp]))<1.0e-12_wp)
  ! Post-Newton projection uses the *post-SI* momentum, not the stage input.
  before_momentum=after_momentum
  after_momentum=[10.0_wp,10.0_wp]
  CALL correct_friction_momentum(before_momentum,after_momentum)
  CALL assert_true('post-Newton correction clamps magnitude', &
       MAXVAL(ABS(after_momentum-before_momentum))<1.0e-12_wp)
  after_momentum=-before_momentum
  CALL correct_friction_momentum(before_momentum,after_momentum)
  CALL assert_true('momentum reversal rejected',ALL(after_momentum==0.0_wp))
  before_momentum=0.0_wp
  after_momentum=[4.0_wp,-5.0_wp]
  CALL correct_friction_momentum(before_momentum,after_momentum)
  CALL assert_true('zero reference leaves zero momentum',ALL(after_momentum==0.0_wp))
  CALL domain%finalize
  CALL state%finalize
  CALL model_layout%finalize

  WRITE(*,*) 'PASS: state and domain initialization verified'

CONTAINS

  SUBROUTINE assert_true(label, condition)

    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition

    IF (.NOT. condition) THEN
       WRITE(*,*) 'FAIL: ', label
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_true

END PROGRAM test_state_initialization
