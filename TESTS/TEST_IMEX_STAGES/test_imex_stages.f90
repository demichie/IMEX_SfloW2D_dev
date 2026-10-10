!> \brief Test-only thread-safe capture of raw production IMEX states.
MODULE imex_test_observer
  USE parameters_2d, ONLY : wp
  IMPLICIT NONE
  REAL(wp), ALLOCATABLE :: known(:,:,:,:), solved(:,:,:,:), raw_final(:,:,:)
  INTEGER, ALLOCATABLE :: seen(:,:,:,:), implicit_statuses(:,:,:)
CONTAINS
  !> \brief Copy one observed cell into exclusively owned test storage.
  !> \param[in] event Audit point: known stage (1), solved stage (2), raw final (3).
  !> \param[in] stage Stage index, zero at the final assembly.
  !> \param[in] j Cell x index.
  !> \param[in] k Cell y index.
  !> \param[in] q_cell Read-only conservative state.
  !> \param[in] implicit_status Local solve status, zero when not called.
  SUBROUTINE capture_state(event,stage,j,k,q_cell,implicit_status)
    INTEGER, INTENT(IN) :: event,stage,j,k,implicit_status
    REAL(wp), INTENT(IN) :: q_cell(:)
    SELECT CASE(event)
    CASE(1)
       known(:,j,k,stage)=q_cell
       seen(1,j,k,stage)=seen(1,j,k,stage)+1
    CASE(2)
       solved(:,j,k,stage)=q_cell
       implicit_statuses(j,k,stage)=implicit_status
       seen(2,j,k,stage)=seen(2,j,k,stage)+1
    CASE(3)
       raw_final(:,j,k)=q_cell
       seen(3,j,k,1)=seen(3,j,k,1)+1
    CASE DEFAULT
       ERROR STOP 'unknown IMEX observation event'
    END SELECT
  END SUBROUTINE capture_state
END MODULE imex_test_observer

!> \brief Exercise real IMEX stages, final repairs and Newton solves on generated fixtures.
PROGRAM test_imex_stages
  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE equation_metadata_2d, ONLY : equation_partition_type
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver, finalize_nonlinear_solver
  USE imex_test_observer
  USE omp_lib
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  TYPE(domain_type) :: domain
  TYPE(spatial_operator_type) :: spatial
  TYPE(time_integration_workspace_type) :: integration
  TYPE(equation_partition_type) :: partition
  REAL(wp), ALLOCATABLE :: q(:,:,:),q0(:,:,:),qp(:,:,:),Z(:,:),statistics(:,:),initial_mask(:,:)
  REAL(wp) :: dt,dt_bound,dt_requested,t,mass0,repair
  CHARACTER(LEN=32) :: argument
  INTEGER :: nx,ny,slope,curvature,limiter_id,input_unit,output_unit,j,k,i,step
  INTEGER :: threads,actual_threads,stages,steps,drag,observe,l

  CALL get_command_argument(1,argument); READ(argument,*) threads
  CALL get_command_argument(2,argument); READ(argument,*) stages
  CALL get_command_argument(3,argument); READ(argument,*) steps
  CALL get_command_argument(4,argument); READ(argument,*) drag
  CALL get_command_argument(5,argument); READ(argument,*) observe
  IF (stages<2 .OR. stages>4 .OR. steps<1) ERROR STOP 'invalid test controls'
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
  READ(input_unit,*) dt_requested
  CALL initialize_fixture
  READ(input_unit,*) B_vertex
  READ(input_unit,*) q
  CLOSE(input_unit)
  q0=q; qp=0.0_wp; qp(4,:,:)=T_ambient; Z=0.0_wp; t=0.0_wp
  ALLOCATE(known(n_vars,nx,ny,n_RK),solved(n_vars,nx,ny,n_RK),raw_final(n_vars,nx,ny))
  ALLOCATE(seen(3,nx,ny,n_RK),implicit_statuses(nx,ny,n_RK),statistics(15,steps),initial_mask(nx,ny))
  CALL refresh_topography_geometry
  mass0=SUM(q(1,:,:))
  DO step=1,steps
     CALL domain%check_solve(q,t,.FALSE.)
     IF (step==1) initial_mask=MERGE(1.0_wp,0.0_wp,domain%solve_mask)
     CALL spatial%compute_timestep(q,qp,t,dt_bound,domain)
     dt=dt_bound
     IF (dt_requested>0.0_wp) dt=dt_requested
     known=0.0_wp; solved=0.0_wp; raw_final=0.0_wp; seen=0; implicit_statuses=0
     IF (observe/=0) THEN
        CALL integration%advance(q,qp,t,dt,Z,partition,domain,spatial,observer=capture_state)
     ELSE
        CALL integration%advance(q,qp,t,dt,Z,partition,domain,spatial)
     END IF
     statistics(:,step)=0.0_wp
     statistics(1,step)=dt; statistics(2,step)=dt_bound
     IF (.NOT.ALL(ieee_is_finite(q))) ERROR STOP 'nonfinite advanced state'
     IF (observe/=0) THEN
        repair=0.0_wp
        DO l=1,domain%solve_cells
           j=domain%j_cent(l); k=domain%k_cent(l)
           IF (ANY(seen(1:2,j,k,:)/=1) .OR. seen(3,j,k,1)/=1) &
                ERROR STOP 'missing or repeated stage observation'
           IF (.NOT.ALL(ieee_is_finite(known(:,j,k,:))) .OR. &
                .NOT.ALL(ieee_is_finite(solved(:,j,k,:))) .OR. &
                .NOT.ALL(ieee_is_finite(raw_final(:,j,k)))) ERROR STOP 'nonfinite raw IMEX state'
           repair=MAX(repair,MAXVAL(ABS(q(:,j,k)-raw_final(:,j,k))))
        END DO
        statistics(3:4,step)=HUGE(1.0_wp)
        DO i=1,n_RK
           statistics(3,step)=MIN(statistics(3,step),MINVAL(known(1,:,:,i),MASK=domain%solve_mask)/rho_l)
           statistics(4,step)=MIN(statistics(4,step),MINVAL(solved(1,:,:,i),MASK=domain%solve_mask)/rho_l)
        END DO
        statistics(5,step)=MINVAL(raw_final(1,:,:),MASK=domain%solve_mask)/rho_l
        statistics(6,step)=repair
        statistics(7,step)=REAL(COUNT(implicit_statuses<0),wp)
        statistics(8,step)=REAL(COUNT(implicit_statuses>0),wp)
        statistics(9,step)=SUM(raw_final(1,:,:))-mass0
        statistics(13,step)=T_ambient
        statistics(14,step)=1.0_wp; statistics(15,step)=1.0_wp
        DO l=1,domain%solve_cells
           j=domain%j_cent(l); k=domain%k_cent(l)
           DO i=1,n_RK
              CALL admissibility_values(known(:,j,k,i),statistics(13:15,step))
              CALL admissibility_values(solved(:,j,k,i),statistics(13:15,step))
           END DO
           CALL admissibility_values(raw_final(:,j,k),statistics(13:15,step))
        END DO
     END IF
     statistics(10,step)=SUM(q(1,:,:))-mass0
     statistics(11,step)=REAL(domain%solve_cells,wp)
     statistics(12,step)=t+dt
     t=t+dt
  END DO
  OPEN(NEWUNIT=output_unit,FILE='result.bin',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE')
  WRITE(output_unit) q0,q,qp,initial_mask,known,solved,raw_final
  WRITE(output_unit) REAL(implicit_statuses,wp),statistics
  CLOSE(output_unit)
  CALL integration%finalize
  CALL spatial%finalize
  CALL partition%finalize
  CALL domain%finalize
  CALL finalize_nonlinear_solver
  WRITE(*,*) 'PASS: raw IMEX observation, finite states and actual OpenMP team'
CONTAINS
  !> \brief Reduce raw temperature and component/carrier fractions without primitive conversion.
  !> \param[in] cell Conservative liquid/zero-solid cell state.
  !> \param[in,out] minima Minimum temperature, solid fraction and carrier fraction for this step.
  SUBROUTINE admissibility_values(cell,minima)
    REAL(wp), INTENT(IN) :: cell(:)
    REAL(wp), INTENT(INOUT) :: minima(3)
    IF (cell(1)>0.0_wp) THEN
       minima(1)=MIN(minima(1),cell(4)/(sp_heat_l*cell(1)))
       minima(2)=MIN(minima(2),cell(5)/cell(1))
       minima(3)=MIN(minima(3),1.0_wp-cell(5)/cell(1))
    END IF
  END SUBROUTINE admissibility_values

  !> \brief Initialize a liquid model and all geometry/source arrays consumed by real IMEX stages.
  SUBROUTINE initialize_fixture
    INTEGER :: i
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    cell_size=dx; dx2=0.5_wp*dx; dy2=0.5_wp*dy
    one_by_dx=1.0_wp/dx; one_by_dy=1.0_wp/dy
    n_vars=5; n_eqns=5; n_solid=1; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    n_RK=stages; idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
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
    ALLOCATE(Z(nx,ny))
    ALLOCATE(B_vertex(nx+1,ny+1),B_cent(nx,ny),B_face_x(nx+1,ny),B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny),B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny),B_second_xy_geom(nx,ny),B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny),grav_coeff_stag_x(nx+1,ny),grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(source_cell(nx,ny),sourceW(nx,ny),sourceE(nx,ny),sourceS(nx,ny),sourceN(nx,ny))
    source_cell=0; sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
    ALLOCATE(B_nodata(nx,ny)); B_nodata=.FALSE.
    ALLOCATE(x_comp(nx),y_comp(ny),x_stag(nx+1),y_stag(ny+1))
    DO i=1,nx+1
       x_stag(i)=REAL(i-1,wp)*dx
       IF (i<=nx) x_comp(i)=(REAL(i,wp)-0.5_wp)*dx
    END DO
    DO i=1,ny+1
       y_stag(i)=REAL(i-1,wp)*dy
       IF (i<=ny) y_comp(i)=(REAL(i,wp)-0.5_wp)*dy
    END DO
    ALLOCATE(cell_source_fractions(nx,ny),cell_fissure_fractions(nx,ny,0))
    ALLOCATE(cell_arc_perim(nx,ny),cell_arc_n_x(nx,ny),cell_arc_n_y(nx,ny))
    cell_source_fractions=0.0_wp; cell_arc_perim=0.0_wp
    cell_arc_n_x=0.0_wp; cell_arc_n_y=0.0_wp
    IF (drag/=0) THEN
       rheology_flag=.TRUE.; rheology_model=6; friction_factor=1.0_wp
    END IF
    CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.])
    CALL initialize_nonlinear_solver
    CALL integration%initialize
    CALL domain%initialize
    CALL spatial%initialize
  END SUBROUTINE initialize_fixture

END PROGRAM test_imex_stages

