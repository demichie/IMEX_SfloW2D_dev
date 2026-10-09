!> \brief Exercise the production spatial operator, active worksets and CFL API without time stepping.
PROGRAM test_spatial_operator
  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type
  USE equation_terms_2d, ONLY : eval_expl_terms
  USE omp_lib
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  IMPLICIT NONE

  TYPE(domain_type) :: domain
  TYPE(spatial_operator_type) :: spatial
  TYPE(reconstruction_workspace_type) :: reconstruction
  TYPE(hyperbolic_workspace_type) :: hyperbolic
  REAL(wp), ALLOCATABLE :: q(:,:,:), q0(:,:,:), qp(:,:,:), term(:,:,:), direct(:,:,:), cell_sources(:,:,:)
  REAL(wp) :: dt, repeated_dt, saved_cfl, capped_dt, no_fissures(0)
  CHARACTER(LEN=32) :: argument
  INTEGER :: nx,ny,slope,curvature,limiter_id,input_unit,output_unit,j,k,threads,actual_threads

  CALL get_command_argument(1,argument)
  READ(argument,*) threads
  actual_threads=0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual_threads=omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual_threads/=threads) ERROR STOP 'incorrect actual OpenMP team'
  WRITE(*,*) 'Actual OpenMP team:',actual_threads

  OPEN(NEWUNIT=input_unit,FILE='fixture.inp',STATUS='OLD',ACTION='READ')
  READ(input_unit,*) nx,ny,slope,curvature,limiter_id
  READ(input_unit,*) dx,dy
  CALL initialize_fixture
  READ(input_unit,*) B_vertex
  READ(input_unit,*) q
  CLOSE(input_unit)
  q0=q
  CALL refresh_topography_geometry
  CALL domain%check_solve(q,0.0_wp,.FALSE.)
  IF (domain%solve_cells<=0 .OR. domain%reconstruction_cells<domain%solve_cells) &
       ERROR STOP 'invalid production worksets'

  ! Only the read halo is deliberately poisoned. The timestep API must refresh
  ! wet states and dry halo entries; unrelated inactive caches are not inputs.
  qp=0.0_wp; qp(4,:,:)=T_ambient
  DO k=1,ny
     DO j=1,nx
        IF (domain%reconstruction_mask(j,k)) qp(:,j,k)=-999.0_wp
     END DO
  END DO
  CALL spatial%compute_timestep(q,qp,0.0_wp,dt,domain)
  IF (.NOT.ieee_is_finite(dt) .OR. dt<=0.0_wp) ERROR STOP 'invalid production timestep'
  IF (ANY(qp(1,:,:) < 0.0_wp)) ERROR STOP 'physical cache was not refreshed'
  CALL same_bits('timestep input is read-only',RESHAPE(q,[SIZE(q)]),RESHAPE(q0,[SIZE(q0)]))
  CALL spatial%evaluate(qp,term,0.0_wp,domain)
  CALL hyperbolic%evaluate_terms(reconstruction,qp,direct,0.0_wp,          &
       domain%solve_cells,domain%j_cent,domain%k_cent,                    &
       domain%solve_interfaces_x,domain%j_stag_x,domain%k_stag_x,          &
       domain%solve_interfaces_y,domain%j_stag_y,domain%k_stag_y)

  ! The API defines the returned term only on solve_cells. Initialize the
  ! remaining diagnostic payload explicitly, rather than reading undefined
  ! inactive entries from an INTENT(OUT) array.
  DO k=1,ny
     DO j=1,nx
        IF (.NOT.domain%solve_mask(j,k)) THEN
           term(:,j,k)=0.0_wp; direct(:,j,k)=0.0_wp
        END IF
     END DO
  END DO
  CALL same_bits('public operator delegates complete spatial term',     &
       RESHAPE(term,[SIZE(term)]),RESHAPE(direct,[SIZE(direct)]))
  IF (.NOT.ALL(ieee_is_finite(term))) ERROR STOP 'nonfinite spatial term'
  cell_sources=0.0_wp
  DO k=1,ny
     DO j=1,nx
        IF (.NOT.domain%solve_mask(j,k)) CYCLE
        CALL eval_expl_terms(B_prime_x_geom(j,k),B_prime_y_geom(j,k),    &
             B_second_xx_geom(j,k),B_second_xy_geom(j,k),              &
             B_second_yy_geom(j,k),grav_coeff(j,k),qp(:,j,k),           &
             cell_sources(:,j,k),0.0_wp,0.0_wp,no_fissures,            &
             0.0_wp,0.0_wp,0.0_wp,dx*dy)
     END DO
  END DO
  IF (.NOT.ALL(ieee_is_finite(cell_sources))) ERROR STOP 'nonfinite local source'
  CALL spatial%compute_timestep(q,qp,0.0_wp,repeated_dt,domain)
  CALL same_bits('CFL refresh/evaluation is repeatable',[dt],[repeated_dt])

  ! Test the existing max_dt cap and prescribed-step sentinel independently
  ! of the reference's explicit SSPRK2 timestep combination policy.
  max_dt=0.5_wp*dt
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('max_dt cap',[capped_dt],[max_dt])
  saved_cfl=cfl; cfl=-1.0_wp
  CALL spatial%compute_timestep(q,qp,0.0_wp,capped_dt,domain)
  CALL same_bits('prescribed timestep sentinel',[capped_dt],[max_dt])
  cfl=saved_cfl

  CALL normalize_inactive_face_payloads
  OPEN(NEWUNIT=output_unit,FILE='result.bin',ACCESS='STREAM',           &
       FORM='UNFORMATTED',STATUS='REPLACE')
  ! Versioned test-only stream order, decoded by compare_reference.py.
  ! No production output or private spatial-operator interface is changed.
  WRITE(output_unit) dt
  WRITE(output_unit) q,qp,term,cell_sources
  WRITE(output_unit) MERGE(1.0_wp,0.0_wp,domain%solve_mask)
  WRITE(output_unit) B_vertex,B_cent,B_face_x,B_face_y
  WRITE(output_unit) B_prime_x_geom,B_prime_y_geom,B_second_xx_geom, &
       B_second_xy_geom,B_second_yy_geom,grav_coeff,                &
       grav_coeff_stag_x,grav_coeff_stag_y
  WRITE(output_unit) hyperbolic%a_interface_xNeg,hyperbolic%a_interface_xPos, &
       hyperbolic%a_interface_yNeg,hyperbolic%a_interface_yPos
  WRITE(output_unit) reconstruction%q_interfaceL,reconstruction%q_interfaceR, &
       reconstruction%q_interfaceB,reconstruction%q_interfaceT
  WRITE(output_unit) reconstruction%qp_interfaceL,reconstruction%qp_interfaceR, &
       reconstruction%qp_interfaceB,reconstruction%qp_interfaceT
  WRITE(output_unit) reconstruction%eta_interfaceL,reconstruction%eta_interfaceR, &
       reconstruction%eta_interfaceB,reconstruction%eta_interfaceT
  CLOSE(output_unit)
  WRITE(*,*) 'Production worksets:',domain%solve_cells,domain%reconstruction_cells
  CALL spatial%finalize
  CALL hyperbolic%finalize
  CALL reconstruction%finalize
  CALL domain%finalize
  WRITE(*,*) 'PASS: production spatial API, active caches, repeatable CFL and caps'

CONTAINS

  !> \brief Configure a one-solid liquid model with zero solid fraction and no physical source except curvature.
  SUBROUTINE initialize_fixture
    INTEGER :: i
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    cell_size=dx; dx2=0.5_wp*dx; dy2=0.5_wp*dy
    one_by_dx=1.0_wp/dx; one_by_dy=1.0_wp/dy
    n_vars=5; n_eqns=5; n_solid=1; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    n_RK=2; idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
    idx_solid_first=5; idx_solid_last=5; idx_solidEqn_first=5; idx_solidEqn_last=5
    idx_add_gas_first=6; idx_add_gas_last=5
    idx_addGasEqn_first=6; idx_addGasEqn_last=5
    idx_u=6; idx_v=7; idx_stoch=0; idx_pore=0
    liquid_flag=.TRUE.; gas_flag=.FALSE.; rheology_flag=.FALSE.; rheology_model=0
    slope_correction_flag=slope/=0; curvature_term_flag=curvature/=0
    stochastic_flag=.FALSE.; stoch_transport_flag=.FALSE.; pore_pressure_flag=.FALSE.
    gas_loss_flag=.FALSE.; liquid_vaporization_flag=.FALSE.; entrainment_flag=.FALSE.
    radial_source_flag=.FALSE.; lateral_source_flag=.FALSE.
    bottom_radial_source_flag=.FALSE.; bottom_fissural_source_flag=.FALSE.
    sutherland_flag=.FALSE.; topo_change_flag=.FALSE.; verbose_level=-1
    n_fissures=0; n_intervals=0
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; rho_a_amb=1.2_wp; grav=9.81_wp
    sp_heat_l=4180.0_wp; sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    T_ambient=300.0_wp; pres=101300.0_wp; inv_pres=1.0_wp/pres
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
    limiter=limiter_id; theta=1.3_wp; reconstr_coeff=1.0_wp
    cfl=0.24_wp; max_dt=1.0E6_wp; maximum_solid_packing=0.6_wp
    ALLOCATE(rho_s(1),inv_rho_s(1),sp_heat_s(1),sp_heat_g(0),sp_gas_const_g(0))
    rho_s=2500.0_wp; inv_rho_s=1.0_wp/rho_s; sp_heat_s=1100.0_wp
    ALLOCATE(bcW(n_vars),bcE(n_vars),bcS(n_vars),bcN(n_vars))
    bcW%flag=1; bcW%value=0.0_wp; bcE%flag=1; bcE%value=0.0_wp
    bcS%flag=1; bcS%value=0.0_wp; bcN%flag=1; bcN%value=0.0_wp
    ALLOCATE(q(n_vars,nx,ny),q0(n_vars,nx,ny),qp(n_vars+2,nx,ny))
    ALLOCATE(term(n_eqns,nx,ny),direct(n_eqns,nx,ny),cell_sources(n_eqns,nx,ny))
    ALLOCATE(B_vertex(nx+1,ny+1),B_cent(nx,ny),B_face_x(nx+1,ny),B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny),B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny),B_second_xy_geom(nx,ny),B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny),grav_coeff_stag_x(nx+1,ny),grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(source_cell(nx,ny),sourceW(nx,ny),sourceE(nx,ny),sourceS(nx,ny),sourceN(nx,ny))
    source_cell=0; sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
    ALLOCATE(x_comp(nx),y_comp(ny),x_stag(nx+1),y_stag(ny+1))
    DO i=1,nx+1
       x_stag(i)=REAL(i-1,wp)*dx
       IF (i<=nx) x_comp(i)=(REAL(i,wp)-0.5_wp)*dx
    END DO
    DO i=1,ny+1
       y_stag(i)=REAL(i-1,wp)*dy
       IF (i<=ny) y_comp(i)=(REAL(i,wp)-0.5_wp)*dy
    END DO
    CALL domain%initialize
    CALL spatial%initialize
    CALL reconstruction%initialize
    CALL hyperbolic%initialize
  END SUBROUTINE initialize_fixture

  !> \brief Zero diagnostic entries outside the production active-face contract before serialization.
  SUBROUTINE normalize_inactive_face_payloads
    DO k=1,ny
       DO j=1,nx+1
          IF (domain%solve_mask_x(j,k)) CYCLE
          reconstruction%q_interfaceL(:,j,k)=0.0_wp; reconstruction%q_interfaceR(:,j,k)=0.0_wp
          reconstruction%qp_interfaceL(:,j,k)=0.0_wp; reconstruction%qp_interfaceR(:,j,k)=0.0_wp
          reconstruction%eta_interfaceL(j,k)=0.0_wp; reconstruction%eta_interfaceR(j,k)=0.0_wp
       END DO
    END DO
    DO k=1,ny+1
       DO j=1,nx
          IF (domain%solve_mask_y(j,k)) CYCLE
          reconstruction%q_interfaceB(:,j,k)=0.0_wp; reconstruction%q_interfaceT(:,j,k)=0.0_wp
          reconstruction%qp_interfaceB(:,j,k)=0.0_wp; reconstruction%qp_interfaceT(:,j,k)=0.0_wp
          reconstruction%eta_interfaceB(j,k)=0.0_wp; reconstruction%eta_interfaceT(j,k)=0.0_wp
       END DO
    END DO
  END SUBROUTINE normalize_inactive_face_payloads

  !> \brief Compare exact bit patterns to test immutable inputs and repeated public API evaluations.
  !> \param[in] label Assertion name.
  !> \param[in] actual,reference Flattened REAL fields to compare.
  SUBROUTINE same_bits(label,actual,reference)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual(:),reference(:)
    IF (ANY(TRANSFER(actual,[0_int64],SIZE(actual)) /= TRANSFER(reference,[0_int64],SIZE(reference)))) THEN
       WRITE(*,*) 'FAIL: ',label
       ERROR STOP 1
    END IF
  END SUBROUTINE same_bits
END PROGRAM test_spatial_operator
