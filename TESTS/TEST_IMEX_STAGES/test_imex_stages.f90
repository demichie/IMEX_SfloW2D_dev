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

!> \brief Exercise real IMEX stages, final repairs and implicit-group solves on generated fixtures.
PROGRAM test_imex_stages
  USE parameters_2d
  USE constitutive_parameters_2d
  USE rheology_parameters_2d, ONLY : xi, T_ref
  USE rheology_voellmy_2d, ONLY : mu
  USE rheology_quadratic_2d, ONLY : friction_factor
  USE rheology_lava_2d, ONLY : nu_ref, visc_par, tau0
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE spatial_operator_2d, ONLY : spatial_operator_type
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE equation_metadata_2d, ONLY : equation_partition_type, IM_GROUP_METHOD_ANALYTIC, &
       IM_GROUP_METHOD_NEWTON
  USE nonlinear_solver_2d, ONLY : initialize_nonlinear_solver, finalize_nonlinear_solver
  USE state_conversion_2d, ONLY : qc_to_qp
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type
  USE imex_test_observer
  USE omp_lib
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  TYPE(domain_type) :: domain
  TYPE(spatial_operator_type) :: spatial
  TYPE(time_integration_workspace_type) :: integration
  TYPE(equation_partition_type) :: partition
  TYPE(reconstruction_workspace_type) :: audit_reconstruction
  TYPE(hyperbolic_workspace_type) :: audit_hyperbolic
  REAL(wp), ALLOCATABLE :: q(:,:,:),q0(:,:,:),qp(:,:,:),Z(:,:),statistics(:,:),initial_mask(:,:)
  REAL(wp) :: dt,dt_bound,dt_requested,t,mass0,repair
  REAL(wp), ALLOCATABLE :: composition_audit(:,:), audit_qp(:,:,:), audit_rhs(:,:,:), inventory0(:)
  CHARACTER(LEN=32) :: argument
  INTEGER :: nx,ny,slope,curvature,limiter_id,input_unit,output_unit,j,k,i,step
  INTEGER :: threads,actual_threads,stages,steps,drag,observe,l
  INTEGER :: composition_mode,boundary_test_mode,source_test_model,analytic_mode
  INTEGER :: mixed_group_mode

  CALL get_command_argument(1,argument); READ(argument,*) threads
  CALL get_command_argument(2,argument); READ(argument,*) stages
  CALL get_command_argument(3,argument); READ(argument,*) steps
  CALL get_command_argument(4,argument); READ(argument,*) drag
  CALL get_command_argument(5,argument); READ(argument,*) observe
  ! The optional sixth control extends the observer, not the production API.
  ! Omitted/zero retains the original liquid payload and numerical fixtures.
  composition_mode=0
  CALL get_command_argument(6,argument)
  IF (LEN_TRIM(argument)>0) READ(argument,*) composition_mode
  IF (composition_mode<0 .OR. composition_mode>3) ERROR STOP 'invalid composition fixture mode'
  boundary_test_mode=0
  CALL get_command_argument(7,argument)
  IF (LEN_TRIM(argument)>0) READ(argument,*) boundary_test_mode
  IF (boundary_test_mode<0 .OR. boundary_test_mode>1) ERROR STOP 'invalid boundary test mode'
  ! Optional D02 source-only integration fixtures; omitted preserves all
  ! historical driver behavior. No alternate production source is introduced.
  source_test_model=0
  CALL get_command_argument(8,argument)
  IF (LEN_TRIM(argument)>0) READ(argument,*) source_test_model
  IF (source_test_model/=0 .AND. source_test_model/=1 .AND. source_test_model/=3) &
       ERROR STOP 'invalid D02 source model'
  analytic_mode=0
  argument=''
  CALL get_command_argument(9,argument)
  IF (LEN_TRIM(argument)>0) READ(argument,*) analytic_mode
  IF (analytic_mode<0 .OR. analytic_mode>1) ERROR STOP 'invalid analytic fixture mode'
  IF (analytic_mode==1 .AND. (drag==0 .OR. composition_mode/=0 .OR. source_test_model/=0)) &
       ERROR STOP 'analytic fixture requires single-phase quadratic drag'
  mixed_group_mode=0
  argument=''
  CALL get_command_argument(10,argument)
  IF (LEN_TRIM(argument)>0) READ(argument,*) mixed_group_mode
  IF (mixed_group_mode<0 .OR. mixed_group_mode>2) ERROR STOP 'invalid mixed-group fixture mode'
  IF (mixed_group_mode>0 .AND. (source_test_model/=3 .OR. composition_mode/=0 .OR. &
       drag/=0 .OR. analytic_mode/=0)) ERROR STOP 'mixed-group fixture requires lava source model'
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
  CALL refresh_topography_geometry
  IF (boundary_test_mode==1) THEN
     CALL run_boundary_reconstruction_tests
     STOP 0
  END IF
  ALLOCATE(known(n_vars,nx,ny,n_RK),solved(n_vars,nx,ny,n_RK),raw_final(n_vars,nx,ny))
  ALLOCATE(seen(3,nx,ny,n_RK),implicit_statuses(nx,ny,n_RK),statistics(15,steps),initial_mask(nx,ny))
  mass0=SUM(q(1,:,:))
  IF (composition_mode>0) THEN
     ALLOCATE(composition_audit(2*(n_vars+1)+5,steps),inventory0(n_vars+1))
     ALLOCATE(audit_qp(n_vars+2,nx,ny),audit_rhs(n_eqns,nx,ny))
     composition_audit=0.0_wp
     DO i=1,n_vars
        inventory0(i)=SUM(q(i,:,:))
     END DO
     inventory0(n_vars+1)=SUM(q(1,:,:))-SUM(q(5:n_vars,:,:))
     CALL audit_reconstruction%initialize
     CALL audit_hyperbolic%initialize
  END IF
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
     IF (composition_mode>0) CALL audit_composition_step
     t=t+dt
  END DO
  OPEN(NEWUNIT=output_unit,FILE='result.bin',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE')
  WRITE(output_unit) q0,q,qp,initial_mask,known,solved,raw_final
  WRITE(output_unit) REAL(implicit_statuses,wp),statistics
  CLOSE(output_unit)
  IF (composition_mode>0) THEN
     OPEN(NEWUNIT=output_unit,FILE='composition.bin',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE')
     WRITE(output_unit) composition_audit
     ! Last-stage traces are retained without primitive postprocessing in Python.
     WRITE(output_unit) audit_reconstruction%qp_cellW,audit_reconstruction%qp_cellE
     WRITE(output_unit) audit_reconstruction%qp_cellS,audit_reconstruction%qp_cellN
     WRITE(output_unit) audit_reconstruction%q_interfaceL,audit_reconstruction%q_interfaceR
     WRITE(output_unit) audit_reconstruction%q_interfaceB,audit_reconstruction%q_interfaceT
     WRITE(output_unit) audit_rhs
     CLOSE(output_unit)
     CALL audit_hyperbolic%finalize
     CALL audit_reconstruction%finalize
  END IF
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
    REAL(wp) :: capacity, residual_carrier
    IF (cell(1)>0.0_wp .AND. composition_mode>0) THEN
       residual_carrier=cell(1)-SUM(cell(5:n_vars))
       capacity=DOT_PRODUCT(cell(5:6),sp_heat_s)
       IF (gas_flag) THEN
          capacity=capacity+cell(7)*sp_heat_g(1)+residual_carrier*sp_heat_a
          IF (liquid_flag) capacity=capacity+cell(8)*sp_heat_l
       ELSE
          capacity=capacity+residual_carrier*sp_heat_l
       END IF
       IF (capacity<=0.0_wp) ERROR STOP 'nonpositive raw mixture heat capacity'
       minima(1)=MIN(minima(1),cell(4)/capacity)
       minima(2)=MIN(minima(2),MINVAL(cell(5:n_vars))/cell(1))
       minima(3)=MIN(minima(3),residual_carrier/cell(1))
    ELSEIF (cell(1)>0.0_wp) THEN
       minima(1)=MIN(minima(1),cell(4)/(sp_heat_l*cell(1)))
       minima(2)=MIN(minima(2),cell(5)/cell(1))
       minima(3)=MIN(minima(3),1.0_wp-cell(5)/cell(1))
    END IF
  END SUBROUTINE admissibility_values

  !> \brief Exercise production linear boundary slopes, source traces and collapsed axes.
  SUBROUTINE run_boundary_reconstruction_tests
    USE equation_terms_2d, ONLY : eval_source_bdry
    TYPE(reconstruction_workspace_type) :: boundary_reconstruction
    REAL(wp), ALLOCATABLE :: boundary_qp(:,:,:)
    INTEGER, ALLOCATABLE :: boundary_j(:),boundary_k(:)
    REAL(wp) :: x_gradient,y_gradient,center_value,boundary_value,slope_value,expected
    REAL(wp) :: source_temperature,source_state(n_vars+2)
    INTEGER :: cell_j,cell_k,cell_l,side,boundary_flag,side_index,unit

    ALLOCATE(boundary_qp(n_vars+2,nx,ny),boundary_j(nx*ny),boundary_k(nx*ny))
    x_gradient=2.0_wp
    y_gradient=3.0_wp
    DO cell_k=1,ny
       DO cell_j=1,nx
          boundary_qp(:,cell_j,cell_k)=0.0_wp
          boundary_qp(1,cell_j,cell_k)=1.0_wp
          boundary_qp(2,cell_j,cell_k)=0.3_wp
          boundary_qp(3,cell_j,cell_k)=-0.2_wp
          boundary_qp(idx_T,cell_j,cell_k)=T_ambient+x_gradient*x_comp(cell_j)+ &
               y_gradient*y_comp(cell_k)
          boundary_qp(idx_solid_first,cell_j,cell_k)=0.2_wp
          boundary_qp(idx_u,cell_j,cell_k)=0.3_wp
          boundary_qp(idx_v,cell_j,cell_k)=-0.2_wp
          cell_l=cell_j+(cell_k-1)*nx
          boundary_j(cell_l)=cell_j
          boundary_k(cell_l)=cell_k
       END DO
    END DO

    limiter(idx_T)=5
    CALL boundary_reconstruction%initialize
    source_cell=0
    sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
    lateral_source_flag=.FALSE.; radial_source_flag=.FALSE.

    DO side=1,4
       IF ((side<=2 .AND. nx<=1) .OR. (side>=3 .AND. ny<=1)) CYCLE
       DO boundary_flag=0,2
          bcW%flag=1; bcE%flag=1; bcS%flag=1; bcN%flag=1
          bcW%value=0.0_wp; bcE%value=0.0_wp
          bcS%value=0.0_wp; bcN%value=0.0_wp
          SELECT CASE(side)
          CASE(1)
             side_index=MAX(1,MIN(ny,(ny+1)/2))
             boundary_value=T_ambient+y_gradient*y_comp(side_index)
             bcW(idx_T)%flag=boundary_flag; bcW(idx_T)%value=boundary_value
             cell_j=1; cell_k=side_index; slope_value=x_gradient
          CASE(2)
             side_index=MAX(1,MIN(ny,(ny+1)/2))
             boundary_value=T_ambient+x_gradient*x_stag(nx+1)+y_gradient*y_comp(side_index)
             bcE(idx_T)%flag=boundary_flag; bcE(idx_T)%value=boundary_value
             cell_j=nx; cell_k=side_index; slope_value=x_gradient
          CASE(3)
             side_index=MAX(1,MIN(nx,(nx+1)/2))
             boundary_value=T_ambient+x_gradient*x_comp(side_index)
             bcS(idx_T)%flag=boundary_flag; bcS(idx_T)%value=boundary_value
             cell_j=side_index; cell_k=1; slope_value=y_gradient
          CASE(4)
             side_index=MAX(1,MIN(nx,(nx+1)/2))
             boundary_value=T_ambient+x_gradient*x_comp(side_index)+y_gradient*y_stag(ny+1)
             bcN(idx_T)%flag=boundary_flag; bcN(idx_T)%value=boundary_value
             cell_j=side_index; cell_k=ny; slope_value=y_gradient
          END SELECT
          IF (boundary_flag==1) THEN
             SELECT CASE(side)
             CASE(1)
                bcW(idx_T)%value=slope_value
             CASE(2)
                bcE(idx_T)%value=slope_value
             CASE(3)
                bcS(idx_T)%value=slope_value
             CASE(4)
                bcN(idx_T)%value=slope_value
             END SELECT
          END IF

          CALL boundary_reconstruction%reconstruct(boundary_qp,0.0_wp,nx*ny, &
               boundary_j,boundary_k,.TRUE.)
          center_value=boundary_qp(idx_T,cell_j,cell_k)
          SELECT CASE(side)
          CASE(1)
             expected=center_value-reconstr_coeff*dx2*slope_value
             IF (boundary_flag==0) expected=boundary_value
             CALL assert_boundary_value('west trace', &
                  boundary_reconstruction%qp_cellW(idx_T,cell_j,cell_k),expected)
             expected=center_value+reconstr_coeff*dx2*slope_value
             CALL assert_boundary_value('west interior trace', &
                  boundary_reconstruction%qp_cellE(idx_T,cell_j,cell_k),expected)
          CASE(2)
             expected=center_value+reconstr_coeff*dx2*slope_value
             IF (boundary_flag==0) expected=boundary_value
             CALL assert_boundary_value('east trace', &
                  boundary_reconstruction%qp_cellE(idx_T,cell_j,cell_k),expected)
             expected=center_value-reconstr_coeff*dx2*slope_value
             CALL assert_boundary_value('east interior trace', &
                  boundary_reconstruction%qp_cellW(idx_T,cell_j,cell_k),expected)
          CASE(3)
             expected=center_value-reconstr_coeff*dy2*slope_value
             IF (boundary_flag==0) expected=boundary_value
             CALL assert_boundary_value('south trace', &
                  boundary_reconstruction%qp_cellS(idx_T,cell_j,cell_k),expected)
             expected=center_value+reconstr_coeff*dy2*slope_value
             CALL assert_boundary_value('south interior trace', &
                  boundary_reconstruction%qp_cellN(idx_T,cell_j,cell_k),expected)
          CASE(4)
             expected=center_value+reconstr_coeff*dy2*slope_value
             IF (boundary_flag==0) expected=boundary_value
             CALL assert_boundary_value('north trace', &
                  boundary_reconstruction%qp_cellN(idx_T,cell_j,cell_k),expected)
             expected=center_value-reconstr_coeff*dy2*slope_value
             CALL assert_boundary_value('north interior trace', &
                  boundary_reconstruction%qp_cellS(idx_T,cell_j,cell_k),expected)
          END SELECT
       END DO
    END DO

    IF (nx==1) THEN
       CALL boundary_reconstruction%reconstruct(boundary_qp,0.0_wp,nx*ny, &
            boundary_j,boundary_k,.TRUE.)
       DO cell_k=1,ny
          CALL assert_boundary_value('collapsed x west', &
               boundary_reconstruction%qp_cellW(idx_T,1,cell_k), &
               boundary_qp(idx_T,1,cell_k))
          CALL assert_boundary_value('collapsed x east', &
               boundary_reconstruction%qp_cellE(idx_T,1,cell_k), &
               boundary_qp(idx_T,1,cell_k))
       END DO
    END IF
    IF (ny==1) THEN
       CALL boundary_reconstruction%reconstruct(boundary_qp,0.0_wp,nx*ny, &
            boundary_j,boundary_k,.TRUE.)
       DO cell_j=1,nx
          CALL assert_boundary_value('collapsed y south', &
               boundary_reconstruction%qp_cellS(idx_T,cell_j,1), &
               boundary_qp(idx_T,cell_j,1))
          CALL assert_boundary_value('collapsed y north', &
               boundary_reconstruction%qp_cellN(idx_T,cell_j,1), &
               boundary_qp(idx_T,cell_j,1))
       END DO
    END IF

    IF (nx>=3 .AND. ny>=3) THEN
       CALL allocate_boundary_source_vectors
       lateral_source_flag=.TRUE.; radial_source_flag=.FALSE.
       h_source=2.0_wp; vel_source=0.25_wp; T_source=350.0_wp
       xs_source(1:n_solid)=0.2_wp
       time_param=[1.0_wp,1.0_wp,0.0_wp,10.0_wp]
       source_temperature=T_source
       CALL eval_source_bdry(0.0_wp,1.0_wp,0.0_wp,source_state)
       CALL assert_boundary_value('source boundary helper temperature', &
            source_state(idx_T),source_temperature)
       DO side=1,2
          source_cell=0
          sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
          sourceW_vect_x=0.0_wp; sourceW_vect_y=0.0_wp
          sourceE_vect_x=0.0_wp; sourceE_vect_y=0.0_wp
          sourceS_vect_x=0.0_wp; sourceS_vect_y=0.0_wp
          sourceN_vect_x=0.0_wp; sourceN_vect_y=0.0_wp
          SELECT CASE(side)
          CASE(1)
             cell_j=1; cell_k=(ny+1)/2
             sourceW(cell_j,cell_k)=.TRUE.; sourceW_vect_x(cell_j,cell_k)=1.0_wp
          CASE(2)
             cell_j=nx; cell_k=(ny+1)/2
             sourceE(cell_j,cell_k)=.TRUE.; sourceE_vect_x(cell_j,cell_k)=-1.0_wp
          END SELECT
          source_cell(cell_j,cell_k)=2
          CALL boundary_reconstruction%reconstruct(boundary_qp,0.0_wp,nx*ny, &
               boundary_j,boundary_k,.TRUE.)
          SELECT CASE(side)
          CASE(1)
             expected=boundary_reconstruction%qp_interfaceR(idx_T,cell_j,cell_k)
          CASE(2)
             expected=boundary_reconstruction%qp_interfaceL(idx_T,cell_j+1,cell_k)
          END SELECT
          CALL assert_boundary_value('lateral source trace',expected,source_temperature)
       END DO

       lateral_source_flag=.FALSE.; radial_source_flag=.TRUE.
       DO side=1,4
          source_cell=0
          sourceW=.FALSE.; sourceE=.FALSE.; sourceS=.FALSE.; sourceN=.FALSE.
          SELECT CASE(side)
          CASE(1)
             cell_j=2; cell_k=(ny+1)/2; sourceW(cell_j,cell_k)=.TRUE.
          CASE(2)
             cell_j=nx-1; cell_k=(ny+1)/2; sourceE(cell_j,cell_k)=.TRUE.
          CASE(3)
             cell_j=(nx+1)/2; cell_k=2; sourceS(cell_j,cell_k)=.TRUE.
          CASE(4)
             cell_j=(nx+1)/2; cell_k=ny-1; sourceN(cell_j,cell_k)=.TRUE.
          END SELECT
          source_cell(cell_j,cell_k)=2
          CALL boundary_reconstruction%reconstruct(boundary_qp,0.0_wp,nx*ny, &
               boundary_j,boundary_k,.TRUE.)
          SELECT CASE(side)
          CASE(1)
             CALL assert_boundary_value('radial west momentum reflection', &
                  boundary_reconstruction%qp_interfaceL(2,cell_j,cell_k), &
                  -boundary_reconstruction%qp_interfaceR(2,cell_j,cell_k))
             CALL assert_boundary_value('radial west velocity reflection', &
                  boundary_reconstruction%qp_interfaceL(idx_u,cell_j,cell_k), &
                  -boundary_reconstruction%qp_interfaceR(idx_u,cell_j,cell_k))
          CASE(2)
             CALL assert_boundary_value('radial east momentum reflection', &
                  boundary_reconstruction%qp_interfaceR(2,cell_j+1,cell_k), &
                  -boundary_reconstruction%qp_interfaceL(2,cell_j+1,cell_k))
             CALL assert_boundary_value('radial east velocity reflection', &
                  boundary_reconstruction%qp_interfaceR(idx_u,cell_j+1,cell_k), &
                  -boundary_reconstruction%qp_interfaceL(idx_u,cell_j+1,cell_k))
          CASE(3)
             CALL assert_boundary_value('radial south momentum reflection', &
                  boundary_reconstruction%qp_interfaceB(3,cell_j,cell_k), &
                  -boundary_reconstruction%qp_interfaceT(3,cell_j,cell_k))
             CALL assert_boundary_value('radial south velocity reflection', &
                  boundary_reconstruction%qp_interfaceB(idx_v,cell_j,cell_k), &
                  -boundary_reconstruction%qp_interfaceT(idx_v,cell_j,cell_k))
          CASE(4)
             CALL assert_boundary_value('radial north momentum reflection', &
                  boundary_reconstruction%qp_interfaceT(3,cell_j,cell_k+1), &
                  -boundary_reconstruction%qp_interfaceB(3,cell_j,cell_k+1))
             CALL assert_boundary_value('radial north velocity reflection', &
                  boundary_reconstruction%qp_interfaceT(idx_v,cell_j,cell_k+1), &
                  -boundary_reconstruction%qp_interfaceB(idx_v,cell_j,cell_k+1))
          END SELECT
       END DO
    END IF

    OPEN(NEWUNIT=unit,FILE='boundary_traces.bin',ACCESS='STREAM', &
         FORM='UNFORMATTED',STATUS='REPLACE')
    WRITE(unit) boundary_reconstruction%qp_cellW,boundary_reconstruction%qp_cellE
    WRITE(unit) boundary_reconstruction%qp_cellS,boundary_reconstruction%qp_cellN
    WRITE(unit) boundary_reconstruction%qp_interfaceL,boundary_reconstruction%qp_interfaceR
    WRITE(unit) boundary_reconstruction%qp_interfaceB,boundary_reconstruction%qp_interfaceT
    CLOSE(unit)
    CALL boundary_reconstruction%finalize
    WRITE(*,*) 'PASS: W/E/S/N boundary modes, sources and collapsed-axis reconstruction'
  END SUBROUTINE run_boundary_reconstruction_tests

  SUBROUTINE allocate_boundary_source_vectors
    IF (.NOT.ALLOCATED(sourceW_vect_x)) THEN
       ALLOCATE(sourceW_vect_x(nx,ny),sourceW_vect_y(nx,ny))
       ALLOCATE(sourceE_vect_x(nx,ny),sourceE_vect_y(nx,ny))
       ALLOCATE(sourceS_vect_x(nx,ny),sourceS_vect_y(nx,ny))
       ALLOCATE(sourceN_vect_x(nx,ny),sourceN_vect_y(nx,ny))
    END IF
  END SUBROUTINE allocate_boundary_source_vectors

  SUBROUTINE assert_boundary_value(label,actual,expected)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual,expected
    REAL(wp) :: tolerance
    tolerance=512.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,ABS(expected))
    IF (ABS(actual-expected)>tolerance) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),actual,expected
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_boundary_value

  !> \brief Initialize a liquid model and all geometry/source arrays consumed by real IMEX stages.
  SUBROUTINE initialize_fixture
    INTEGER :: i
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    cell_size=dx; dx2=0.5_wp*dx; dy2=0.5_wp*dy
    one_by_dx=1.0_wp/dx; one_by_dy=1.0_wp/dy
    n_vars=5; n_eqns=5; n_solid=1; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    kt_trace_reconstruction = 'froude_blend'
    n_RK=stages; idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
    idx_solid_first=5; idx_solid_last=5; idx_solidEqn_first=5; idx_solidEqn_last=5
    idx_add_gas_first=6; idx_add_gas_last=5
    idx_addGasEqn_first=6; idx_addGasEqn_last=5
    idx_u=6; idx_v=7; idx_stoch=0; idx_pore=0
    liquid_flag=.TRUE.; gas_flag=.FALSE.; rheology_flag=.FALSE.; rheology_model_id=0
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
    IF (composition_mode>0) THEN
       ! Genuine current closures: two distinct solids, optionally additional
       ! gas and explicit liquid. The untransported remainder is liquid or air.
       n_solid=2; n_add_gas=0
       gas_flag=composition_mode>=2; liquid_flag=composition_mode/=2
       IF (gas_flag) n_add_gas=1
       n_vars=4+n_solid+n_add_gas
       IF (gas_flag .AND. liquid_flag) n_vars=n_vars+1
       n_eqns=n_vars; idx_u=n_vars+1; idx_v=n_vars+2
       idx_solid_last=6; idx_solidEqn_last=6
       idx_add_gas_first=7; idx_add_gas_last=6+n_add_gas
       idx_addGasEqn_first=7; idx_addGasEqn_last=6+n_add_gas
       idx_totMassEqn=1; idx_uEqn=2; idx_vEqn=3; idx_engyEqn=4
       DEALLOCATE(rho_s,inv_rho_s,sp_heat_s,sp_heat_g,sp_gas_const_g)
       ALLOCATE(rho_s(2),inv_rho_s(2),sp_heat_s(2),sp_heat_g(n_add_gas),sp_gas_const_g(n_add_gas))
       rho_s=[2500.0_wp,1800.0_wp]; inv_rho_s=1.0_wp/rho_s
       sp_heat_s=[900.0_wp,1100.0_wp]
       sp_heat_a=1005.0_wp; pres=101325.0_wp; inv_pres=1.0_wp/pres
       rho_a_amb=pres/(sp_gas_const_a*T_ambient)
       sp_heat_g=1850.0_wp; sp_gas_const_g=461.5_wp
    END IF
    ALLOCATE(bcW(n_vars),bcE(n_vars),bcS(n_vars),bcN(n_vars))
    bcW%flag=1; bcW%value=0.0_wp; bcE%flag=1; bcE%value=0.0_wp
    bcS%flag=1; bcS%value=0.0_wp; bcN%flag=1; bcN%value=0.0_wp
    IF (composition_mode>0) THEN
       bcW(2)%flag=0; bcE(2)%flag=0
       bcS(3)%flag=0; bcN(3)%flag=0
    END IF
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
       rheology_flag=.TRUE.; rheology_model_id=6; friction_factor=1.0_wp
    END IF
    IF (source_test_model/=0) THEN
       rheology_flag=.TRUE.; rheology_model_id=source_test_model
       mu=0.3_wp; xi=500.0_wp; tau0=20.0_wp
       nu_ref=0.02_wp; visc_par=0.003_wp;        T_ref=900.0_wp; T_env=300.0_wp
       radiative_term_coeff=3.0E-8_wp; convective_term_coeff=4.0_wp
       IF (mixed_group_mode>0) THEN
          visc_par=0.3_wp
          T_ref=1000.0_wp
          convective_term_coeff=1.0E4_wp
       END IF
       IF (source_test_model==3) kt_trace_reconstruction='h_x'
    END IF
    IF (composition_mode==0) THEN
       IF (analytic_mode==1) THEN
          CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.], &
               [0,1,1,0,0], [IM_GROUP_METHOD_ANALYTIC])
       ELSE IF (mixed_group_mode==1) THEN
          CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.], &
               [0,1,1,2,0], [IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON])
       ELSE IF (mixed_group_mode==2) THEN
          CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.], &
               [0,2,2,1,0], [IM_GROUP_METHOD_NEWTON, IM_GROUP_METHOD_ANALYTIC])
       ELSE
          CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.])
       END IF
       IF (source_test_model==3 .AND. mixed_group_mode==0) &
            CALL partition%initialize([.FALSE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.])
    ELSE
       BLOCK
         LOGICAL :: mask(n_vars)
         mask=.FALSE.; mask(2:3)=.TRUE.
         CALL partition%initialize(mask)
       END BLOCK
    END IF
    CALL initialize_nonlinear_solver
    CALL integration%initialize
    CALL domain%initialize
    CALL spatial%initialize
  END SUBROUTINE initialize_fixture

  !> \brief Audit all raw-stage inventories and independent reconstructed traces for one step.
  !> \details Test-only workspace is allocated once. Serial diagnostics do not
  !>          write production stage/source storage or change the active workset.
  SUBROUTINE audit_composition_step
    REAL(wp) :: pressure_dynamic, face_error, boundary_flux
    INTEGER :: stage_id, cell_j, cell_k
    IF (observe==0 .OR. domain%solve_cells/=nx*ny) ERROR STOP 'composition audit requires wet observed cells'
    DO stage_id=1,n_RK
       CALL audit_inventory(known(:,:,:,stage_id),1)
       CALL audit_inventory(solved(:,:,:,stage_id),1)
       DO cell_k=1,ny
          DO cell_j=1,nx
             CALL qc_to_qp(solved(:,cell_j,cell_k,stage_id),audit_qp(:,cell_j,cell_k),pressure_dynamic)
          END DO
       END DO
       audit_rhs=0.0_wp
       CALL audit_hyperbolic%evaluate_terms(audit_reconstruction,audit_qp,audit_rhs,t, &
            domain%solve_cells,domain%j_cent,domain%k_cent,domain%solve_interfaces_x, &
            domain%j_stag_x,domain%k_stag_x,domain%solve_interfaces_y,domain%j_stag_y,domain%k_stag_y)
       ! Test reconstructed fractions BEFORE final_hp_state closure, plus
       ! conservative face admissibility after the actual production conversion.
       face_error=0.0_wp
       DO cell_k=1,ny
          DO cell_j=1,nx
             CALL audit_face(audit_reconstruction%qp_cellW(:,cell_j,cell_k),cell_j,cell_k,1,face_error)
             CALL audit_face(audit_reconstruction%qp_cellE(:,cell_j,cell_k),cell_j,cell_k,1,face_error)
             CALL audit_face(audit_reconstruction%qp_cellS(:,cell_j,cell_k),cell_j,cell_k,2,face_error)
             CALL audit_face(audit_reconstruction%qp_cellN(:,cell_j,cell_k),cell_j,cell_k,2,face_error)
          END DO
       END DO
       composition_audit(2*(n_vars+1)+1,step)=MAX(composition_audit(2*(n_vars+1)+1,step),face_error)
       boundary_flux=0.0_wp
       IF (nx>1) boundary_flux=MAX(MAXVAL(ABS(audit_hyperbolic%G_interface_xL(1,1,:))), &
            MAXVAL(ABS(audit_hyperbolic%G_interface_xR(1,nx+1,:))))
       IF (ny>1) boundary_flux=MAX(boundary_flux,MAXVAL(ABS(audit_hyperbolic%G_interface_yB(1,:,1))), &
            MAXVAL(ABS(audit_hyperbolic%G_interface_yT(1,:,ny+1))))
       composition_audit(2*(n_vars+1)+2,step)=MAX(composition_audit(2*(n_vars+1)+2,step),boundary_flux)
       composition_audit(2*(n_vars+1)+3,step)=MAX(composition_audit(2*(n_vars+1)+3,step), &
            MAXVAL(ABS(audit_rhs(2:3,:,:))))
       DO cell_k=1,ny
          DO cell_j=1,nx+1
             CALL audit_conservative_face(audit_reconstruction%q_interfaceL(:,cell_j,cell_k))
             CALL audit_conservative_face(audit_reconstruction%q_interfaceR(:,cell_j,cell_k))
          END DO
       END DO
       DO cell_k=1,ny+1
          DO cell_j=1,nx
             CALL audit_conservative_face(audit_reconstruction%q_interfaceB(:,cell_j,cell_k))
             CALL audit_conservative_face(audit_reconstruction%q_interfaceT(:,cell_j,cell_k))
          END DO
       END DO
    END DO
    CALL audit_inventory(raw_final,1)
    CALL audit_inventory(q,2)
  END SUBROUTINE audit_composition_step

  !> \brief Record maximum raw/final signed inventory errors including the untransported carrier.
  !> \param[in] state Conservative state at one audit point.
  !> \param[in] block_id One for raw stages/final; two for the assembled final state.
  SUBROUTINE audit_inventory(state,block_id)
    REAL(wp), INTENT(IN) :: state(:,:,:)
    INTEGER, INTENT(IN) :: block_id
    REAL(wp) :: total(n_vars+1)
    INTEGER :: component, first
    DO component=1,n_vars
       total(component)=SUM(state(component,:,:))
    END DO
    total(n_vars+1)=SUM(state(1,:,:))-SUM(state(5:n_vars,:,:))
    first=(block_id-1)*(n_vars+1)+1
    composition_audit(first:first+n_vars,step)=MAX(composition_audit(first:first+n_vars,step),ABS(total-inventory0))
  END SUBROUTINE audit_inventory

  !> \brief Measure positivity, carrier closure and local TVD bounds of pre-closure face fractions.
  !> \param[in] physical Raw cell-owned primitive face trace.
  !> \param[in] cell_j Owner x index.
  !> \param[in] cell_k Owner y index.
  !> \param[in] axis Normal coordinate, one for x and two for y.
  !> \param[in,out] error Maximum absolute fraction-bound violation over faces.
  SUBROUTINE audit_face(physical,cell_j,cell_k,axis,error)
    REAL(wp), INTENT(IN) :: physical(:)
    INTEGER, INTENT(IN) :: cell_j,cell_k,axis
    REAL(wp), INTENT(INOUT) :: error
    REAL(wp) :: low,high
    INTEGER :: component
    IF (.NOT.ALL(ieee_is_finite(physical))) ERROR STOP 'nonfinite primitive face'
    error=MAX(error,-physical(1),-MINVAL(physical(5:n_vars)),SUM(physical(5:n_vars))-1.0_wp)
    DO component=5,n_vars
       IF (axis==1) THEN
          low=MINVAL(audit_qp(component,MAX(1,cell_j-1):MIN(nx,cell_j+1),cell_k))
          high=MAXVAL(audit_qp(component,MAX(1,cell_j-1):MIN(nx,cell_j+1),cell_k))
       ELSE
          low=MINVAL(audit_qp(component,cell_j,MAX(1,cell_k-1):MIN(ny,cell_k+1)))
          high=MAXVAL(audit_qp(component,cell_j,MAX(1,cell_k-1):MIN(ny,cell_k+1)))
       END IF
       error=MAX(error,low-physical(component),physical(component)-high)
    END DO
  END SUBROUTINE audit_face

  !> \brief Record fraction and temperature admissibility of conservative face states.
  !> \param[in] state Conservative face after production primitive closure/conversion.
  SUBROUTINE audit_conservative_face(state)
    REAL(wp), INTENT(IN) :: state(:)
    REAL(wp) :: minima(3)
    IF (.NOT.ALL(ieee_is_finite(state))) ERROR STOP 'nonfinite conservative face'
    minima=[HUGE(1.0_wp),1.0_wp,1.0_wp]
    CALL admissibility_values(state,minima)
    composition_audit(2*(n_vars+1)+4,step)=MAX(composition_audit(2*(n_vars+1)+4,step), &
         -MIN(minima(2),minima(3)),-state(1)/MAX(1.0_wp,MAXVAL(q0(1,:,:))))
    composition_audit(2*(n_vars+1)+5,step)=MAX(composition_audit(2*(n_vars+1)+5,step), &
         MAX(0.0_wp,273.0_wp-minima(1)))
  END SUBROUTINE audit_conservative_face

END PROGRAM test_imex_stages
