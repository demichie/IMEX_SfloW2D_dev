PROGRAM test_curvature_source

  USE parameters_2d
  USE constitutive_parameters_2d, ONLY : rho_l, inv_rho_l, rho_a_amb, grav, &
       sp_heat_l, inv_rho_s, sp_heat_s, sp_heat_g, sp_gas_const_g
  USE geometry_2d, ONLY : bottom_source_area_ratio, pi_g
  USE rheology_2d, ONLY : curvature_acceleration
  USE source_terms_2d, ONLY : eval_geometry_curvature_explicit
  USE source_terms_2d, ONLY : eval_source_explicit
  USE mass_injection_coefficients_2d, ONLY : set_source_time_coefficients
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT, SOURCE_CADENCE_STAGE
  USE source_registry_2d, ONLY : source_group_type

  IMPLICIT NONE

  REAL(wp), PARAMETER :: Bx = 0.35_wp
  REAL(wp), PARAMETER :: By = -0.22_wp
  REAL(wp), PARAMETER :: Bxx = 0.08_wp
  REAL(wp), PARAMETER :: Bxy = -0.035_wp
  REAL(wp), PARAMETER :: Byy = 0.06_wp
  REAL(wp), PARAMETER :: G = 0.42_wp
  REAL(wp), PARAMETER :: mixture_mass = 1800.0_wp
  REAL(wp), PARAMETER :: u = 2.4_wp
  REAL(wp), PARAMETER :: v = -1.7_wp

  REAL(wp) :: acceleration, expected_acceleration
  REAL(wp) :: source_x, source_y, expected_x, expected_y
  REAL(wp) :: tolerance, ramp_mean
  REAL(wp) :: physical(7), original(7), explicit_source(5)
  TYPE(source_group_type) :: curvature_group, circular_injection_group
  TYPE(source_group_type) :: fissural_injection_group, radial_injection_group
  INTEGER :: sample

  curvature_group%method=SOURCE_METHOD_EXPLICIT
  curvature_group%cadence=SOURCE_CADENCE_STAGE
  curvature_group%source_names=[CHARACTER(LEN=64)::'geometry_curvature']
  circular_injection_group%method=SOURCE_METHOD_EXPLICIT
  circular_injection_group%cadence=SOURCE_CADENCE_STAGE
  circular_injection_group%source_names=[CHARACTER(LEN=64)::'bottom_circular_injection']
  fissural_injection_group%method=SOURCE_METHOD_EXPLICIT
  fissural_injection_group%cadence=SOURCE_CADENCE_STAGE
  fissural_injection_group%source_names=[CHARACTER(LEN=64)::'bottom_fissural_injection']
  radial_injection_group%method=SOURCE_METHOD_EXPLICIT
  radial_injection_group%cadence=SOURCE_CADENCE_STAGE
  radial_injection_group%source_names=[CHARACTER(LEN=64)::'radial_injection']

  expected_acceleration = Bxx*u*u + 2.0_wp*Bxy*u*v + Byy*v*v
  CALL assert_small( 'curvature acceleration',                             &
       ABS(curvature_acceleration(Bxx,Bxy,Byy,u,v)-expected_acceleration),  &
       64.0_wp*EPSILON(1.0_wp)*ABS(expected_acceleration) )

  acceleration = curvature_acceleration(Bxx,Bxy,Byy,u,v)
  CALL eval_geometry_curvature_explicit( Bx, By, Bxx, Bxy, Byy, G,           &
       mixture_mass, u, v, source_x, source_y )

  expected_x = -G*mixture_mass*expected_acceleration*Bx
  expected_y = -G*mixture_mass*expected_acceleration*By
  tolerance = 128.0_wp*EPSILON(1.0_wp)                                    &
       * MAX(1.0_wp,ABS(expected_x),ABS(expected_y))

  CALL assert_small('returned acceleration',                              &
       ABS(acceleration-expected_acceleration),tolerance)
  CALL assert_small('x curvature source',ABS(source_x-expected_x),tolerance)
  CALL assert_small('y curvature source',ABS(source_y-expected_y),tolerance)

  ! The mixed-Hessian contribution must not be dropped in two dimensions.
  CALL assert_small( 'mixed Hessian contribution',                         &
       ABS(expected_acceleration-(Bxx*u*u+Byy*v*v)-2.0_wp*Bxy*u*v),         &
       64.0_wp*EPSILON(1.0_wp)*ABS(expected_acceleration) )

  ! Exercise the actual local-source caller, not a duplicate test-only guard.
  n_vars=5; n_eqns=5; n_solid=1; n_add_gas=0; n_fissures=0
  n_stoch_vars=0; n_pore_vars=0
  CALL configure_canonical_state_indices()
  liquid_flag=.TRUE.; gas_flag=.FALSE.; curvature_term_flag=.TRUE.
  bottom_radial_source_flag=.FALSE.; bottom_fissural_source_flag=.FALSE.
  radial_source_flag=.FALSE.; pore_pressure_flag=.FALSE.
  stoch_transport_flag=.FALSE.; n_intervals=0
  rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; rho_a_amb=1.2_wp; grav=9.81_wp
  sp_heat_l=4180.0_wp
  ALLOCATE(inv_rho_s(1),sp_heat_s(1),sp_heat_g(0),sp_gas_const_g(0))
  inv_rho_s=1.0_wp/2500.0_wp; sp_heat_s=1100.0_wp
  physical=0.0_wp; physical(4)=300.0_wp
  DO sample=1,2
     physical(1)=dry_thickness_tolerance*REAL(sample,wp)/2.0_wp
     physical(6)=1.0E150_wp; physical(7)=-1.0E150_wp
     original=physical
     CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
          0.0_wp,0.0_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
          curvature_group)
     CALL assert_small('dry curvature source',MAXVAL(ABS(explicit_source)),0.0_wp)
     CALL assert_small('dry physical state retained',MAXVAL(ABS(physical-original)),0.0_wp)
  END DO
  physical(1)=2.0_wp*dry_thickness_tolerance; physical(6)=u; physical(7)=v
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.0_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       curvature_group)
  CALL eval_geometry_curvature_explicit(Bx,By,Bxx,Bxy,Byy,G,          &
       rho_l*physical(1),u,v,expected_x,expected_y)
  tolerance=128.0_wp*EPSILON(1.0_wp)*MAX(ABS(expected_x),ABS(expected_y))
  CALL assert_small('resolved x source retained',ABS(explicit_source(2)-expected_x),tolerance)
  CALL assert_small('resolved y source retained',ABS(explicit_source(3)-expected_y),tolerance)
  IF (explicit_source(2)==0.0_wp .OR. explicit_source(3)==0.0_wp) ERROR STOP 1

  ! The curvature guard must not disable injection into an unresolved cell.
  xs_source=0.0_wp; xg_source=0.0_wp
  bottom_radial_source_flag=.TRUE.; vel_source=0.01_wp; T_source=300.0_wp
  time_param=[10.0_wp,5.0_wp,0.0_wp,0.0_wp]
  physical(1)=dry_thickness_tolerance/2.0_wp
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       curvature_group)
  CALL assert_small('unselected circular injection excluded', &
       MAXVAL(ABS(explicit_source)),0.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       circular_injection_group)
  expected_x=0.25_wp*vel_source*rho_l
  CALL assert_small('dry-cell mass injection retained',ABS(explicit_source(1)-expected_x), &
       64.0_wp*EPSILON(1.0_wp)*expected_x)
  expected_y=expected_x*sp_heat_l*T_source
  CALL assert_small('dry-cell thermal injection retained',ABS(explicit_source(4)-expected_y), &
       64.0_wp*EPSILON(1.0_wp)*expected_y)
  CALL assert_small('dry curvature with injection',MAXVAL(ABS(explicit_source(2:3))),0.0_wp)

  ! Step averaging must integrate activation across the interval, not sample
  ! only its active start or endpoint.
  bottom_source_area_ratio=1.0_wp
  pi_g=ACOS(-1.0_wp)
  time_param=[10.0_wp,5.0_wp,0.0_wp,0.0_wp]
  CALL set_source_time_coefficients(4.0_wp,2.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       4.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       circular_injection_group)
  expected_x=0.5_wp*0.25_wp*vel_source*rho_l
  CALL assert_small('partially active step-averaged mass source', &
       ABS(explicit_source(1)-expected_x),64.0_wp*EPSILON(1.0_wp)*expected_x)

  time_param=[10.0_wp,5.0_wp,1.0_wp,0.0_wp]
  CALL set_source_time_coefficients(0.0_wp,2.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       circular_injection_group)
  expected_x=0.75_wp*0.25_wp*vel_source*rho_l
  CALL assert_small('linearly ramped step-averaged mass source', &
       ABS(explicit_source(1)-expected_x),64.0_wp*EPSILON(1.0_wp)*expected_x)

  time_param=[10.0_wp,5.0_wp,0.0_wp,0.0_wp]
  CALL set_source_time_coefficients(9.0_wp,2.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       9.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       circular_injection_group)
  expected_x=0.5_wp*0.25_wp*vel_source*rho_l
  CALL assert_small('period-wrapped step-averaged mass source', &
       ABS(explicit_source(1)-expected_x),64.0_wp*EPSILON(1.0_wp)*expected_x)

  CALL set_source_time_coefficients(5.0_wp,2.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       5.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       circular_injection_group)
  CALL assert_small('inactive step-averaged mass source', &
       MAXVAL(ABS(explicit_source)),0.0_wp)

  bottom_radial_source_flag=.FALSE.; bottom_fissural_source_flag=.TRUE.
  n_fissures=1
  ALLOCATE(xs_fissures(1,1),xg_fissures(0,1))
  xs_fissures=0.0_wp; xg_fissures=0.0_wp
  linear_vel_fissures(1)=0.02_wp; T_fissures(1)=T_source
  time_param_fissures(:,1)=[10.0_wp,5.0_wp,0.0_wp,0.0_wp]
  CALL set_source_time_coefficients(4.0_wp,2.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       4.0_wp,0.0_wp,[0.25_wp],0.0_wp,0.0_wp,0.0_wp,1.0_wp, &
       fissural_injection_group)
  expected_x=0.5_wp*0.25_wp*linear_vel_fissures(1)*rho_l
  CALL assert_small('partially active fissural source', &
       ABS(explicit_source(1)-expected_x),64.0_wp*EPSILON(1.0_wp)*expected_x)

  bottom_fissural_source_flag=.FALSE.; radial_source_flag=.TRUE.
  h_source=2.0_wp; vel_source=0.01_wp
  time_param=[10.0_wp,4.0_wp,2.0_wp,3.0_wp]
  CALL set_source_time_coefficients(0.0_wp,1.0_wp)
  CALL eval_source_explicit(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.0_wp,[REAL(wp)::],2.0_wp,1.0_wp,0.0_wp,1.0_wp, &
       radial_injection_group)
  ramp_mean=0.5_wp-1.0_wp/ACOS(-1.0_wp)
  expected_x=ramp_mean*h_source*vel_source*2.0_wp*rho_l
  CALL assert_small('cosine-ramped radial source', &
       ABS(explicit_source(1)-expected_x),64.0_wp*EPSILON(1.0_wp)*expected_x)

  WRITE(*,*) 'PASS: curvature, dry-state guard, and step-averaged source windows verified'

CONTAINS

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF ( value .GT. limit_value ) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,' > ',limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_curvature_source
