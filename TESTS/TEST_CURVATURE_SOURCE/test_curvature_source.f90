PROGRAM test_curvature_source

  USE parameters_2d
  USE constitutive_parameters_2d, ONLY : rho_l, inv_rho_l, rho_a_amb, grav, &
       sp_heat_l, inv_rho_s, sp_heat_s, sp_heat_g, sp_gas_const_g
  USE equation_terms_2d, ONLY : curvature_acceleration
  USE equation_terms_2d, ONLY : eval_curvature_momentum_source
  USE equation_terms_2d, ONLY : eval_expl_terms

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
  REAL(wp) :: tolerance
  REAL(wp) :: physical(7), original(7), explicit_source(5)
  INTEGER :: sample

  expected_acceleration = Bxx*u*u + 2.0_wp*Bxy*u*v + Byy*v*v
  CALL assert_small( 'curvature acceleration',                             &
       ABS(curvature_acceleration(Bxx,Bxy,Byy,u,v)-expected_acceleration),  &
       64.0_wp*EPSILON(1.0_wp)*ABS(expected_acceleration) )

  acceleration = curvature_acceleration(Bxx,Bxy,Byy,u,v)
  CALL eval_curvature_momentum_source( Bx, By, Bxx, Bxy, Byy, G,           &
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
  idx_solid_first=5; idx_solid_last=5
  idx_add_gas_first=6; idx_add_gas_last=5; idx_u=6; idx_v=7
  idx_solidEqn_first=5; idx_solidEqn_last=5
  idx_addGasEqn_first=6; idx_addGasEqn_last=5
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
     CALL eval_expl_terms(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
          0.0_wp,0.0_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp)
     CALL assert_small('dry curvature source',MAXVAL(ABS(explicit_source)),0.0_wp)
     CALL assert_small('dry physical state retained',MAXVAL(ABS(physical-original)),0.0_wp)
  END DO
  physical(1)=2.0_wp*dry_thickness_tolerance; physical(6)=u; physical(7)=v
  CALL eval_expl_terms(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.0_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp)
  CALL eval_curvature_momentum_source(Bx,By,Bxx,Bxy,Byy,G,          &
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
  CALL eval_expl_terms(Bx,By,Bxx,Bxy,Byy,G,physical,explicit_source, &
       0.0_wp,0.25_wp,[REAL(wp)::],0.0_wp,0.0_wp,0.0_wp,1.0_wp)
  expected_x=0.25_wp*vel_source*rho_l
  CALL assert_small('dry-cell mass injection retained',ABS(explicit_source(1)-expected_x), &
       64.0_wp*EPSILON(1.0_wp)*expected_x)
  expected_y=expected_x*sp_heat_l*T_source
  CALL assert_small('dry-cell thermal injection retained',ABS(explicit_source(4)-expected_y), &
       64.0_wp*EPSILON(1.0_wp)*expected_y)
  CALL assert_small('dry curvature with injection',MAXVAL(ABS(explicit_source(2:3))),0.0_wp)

  WRITE(*,*) 'PASS: slope-weighted curvature and shared dry-depth guard verified'

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
