!> \brief N8-B prescribed-rate refinement and seeded repeated production geometry checks.
PROGRAM test_projection_stress
  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE state_2d, ONLY : state_type
  USE domain_2d, ONLY : domain_type
  USE source_registry_2d, ONLY : source_registry_type,source_group_type
  USE source_catalog_2d, ONLY : SOURCE_METHOD_EXPLICIT
  USE lie_source_executor_2d, ONLY : execute_lie_source_groups
  USE state_conversion_2d, ONLY : qc_to_qp
  USE omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite,ieee_value,ieee_quiet_nan
  IMPLICIT NONE
  INTEGER, PARAMETER :: ns=2, nq=6, rough_steps=32
  REAL(wp), PARAMETER :: tolerance=2048.0_wp*EPSILON(1.0_wp)
  INTEGER :: nx,ny,dimension,level,step,unit,requested,actual,j,k
  INTEGER(int64) :: seed
  CHARACTER(LEN=16) :: arg
  REAL(wp), ALLOCATABLE :: q(:,:,:),qp(:,:,:),rate(:,:),nodes(:,:),weights(:,:)
  TYPE(domain_type) :: domain
  TYPE(state_type) :: state

  CALL get_command_argument(1,arg)
  READ(arg,*) requested
  actual=0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual=omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual/=requested) ERROR STOP 'incorrect actual team'
  WRITE(*,*) 'Actual OpenMP team:',actual
  OPEN(NEWUNIT=unit,FILE='projection_snapshot.bin',ACCESS='STREAM', &
       FORM='UNFORMATTED',STATUS='REPLACE')
  CALL initialize_physics
  DO dimension=1,2
     DO level=0,2
        nx=16*2**level; ny=1
        IF (dimension==2) ny=12*2**level
        CALL allocate_grid
        CALL reset_state(.FALSE.)
        IF (level==0) THEN
           CALL check_lie_exchange_transaction
           CALL check_lie_rejected_transaction
        END IF
        DO k=1,ny
           DO j=1,nx
              rate(j,k)=-0.001_wp*(2.0_wp+COS(2.0_wp*ACOS(-1.0_wp)*(REAL(j,wp)-0.5_wp)*dx))
              IF (dimension==2) rate(j,k)=-0.001_wp*(2.0_wp+ &
                   COS(2.0_wp*ACOS(-1.0_wp)*(REAL(j,wp)-0.5_wp)*dx)* &
                   COS(2.0_wp*ACOS(-1.0_wp)*(REAL(k,wp)-0.5_wp)*dy))
           END DO
        END DO
        CALL projection_update('smooth',1)
        CALL reset_state(.TRUE.)
        rate=0.0_wp
        WHERE(q(1,:,:)>0.0_wp) rate=-0.001_wp
        CALL projection_update('wet_dry',1)

        ! Signed, cell-scale random rates exercise the public production
        ! projector without pretending that a rough field must converge.
        CALL reset_state(.FALSE.)
        seed=104729_int64
        DO step=1,rough_steps
           CALL random_rates(.TRUE.)
           CALL projection_update('rough',step)
        END DO

        ! Independently seeded repeated geometry updates for the N8-B oracle.
        CALL reset_state(.FALSE.)
        seed=130363_int64
        DO step=1,64
           CALL random_rates(.TRUE.)
           CALL projection_update('repeated',step)
        END DO
        CALL release_grid
     END DO
  END DO
  CLOSE(unit)
  WRITE(*,*) 'PASS: Lie-compatible prescribed projection refinement and geometry checks'
CONTAINS

  !> \brief Configure the existing two-solid thermal-energy liquid exchange law.
  SUBROUTINE initialize_physics
    n_solid=ns; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    kt_trace_reconstruction = 'froude_blend'
    n_vars=nq; n_eqns=nq
    CALL configure_canonical_state_indices()
    gas_flag=.FALSE.; liquid_flag=.TRUE.; rheology_flag=.FALSE.; rheology_model='none'
    stochastic_flag=.FALSE.; stoch_transport_flag=.FALSE.; pore_pressure_flag=.FALSE.
    gas_loss_flag=.FALSE.; liquid_vaporization_flag=.FALSE.; entrainment_flag=.FALSE.
    sutherland_flag=.FALSE.; slope_correction_flag=.TRUE.; topo_change_flag=.TRUE.
    n_fissures=0
    erodible_deposit_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.; radial_source_flag=.FALSE.
    erosion_coeff=1000.0_wp; settling_flag=.FALSE.; alphastot_min=0.0_wp
    erodible_porosity=0.25_wp; coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
    grav=9.81_wp; rho_a_amb=1.2_wp; T_ambient=300.0_wp; T_erodible=317.0_wp
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; rho_c_sub=rho_l
    sp_heat_l=4180.0_wp; sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    pres=101300.0_wp; inv_pres=1.0_wp/pres
    kin_visc_c=1.0E-6_wp; hydraulic_permeability=0.0_wp; maximum_solid_packing=0.6_wp
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4; verbose_level=0
    ALLOCATE(rho_s(ns),inv_rho_s(ns),sp_heat_s(ns),diam_s(ns),erodible_fract(ns))
    ALLOCATE(sp_heat_g(0),sp_gas_const_g(0),loss_rate)
    rho_s=[2000.0_wp,3000.0_wp]; inv_rho_s=1.0_wp/rho_s
    sp_heat_s=[900.0_wp,1100.0_wp]; diam_s=[1.0E-4_wp,2.0E-4_wp]
    erodible_fract=[0.6_wp,0.4_wp]; loss_rate=0.0_wp
  END SUBROUTINE initialize_physics

  !> \brief Allocate the rectangular production grid and independent oracle buffers.
  SUBROUTINE allocate_grid
    comp_cells_x=nx; comp_cells_y=ny; comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    dx=1.0_wp/REAL(nx,wp); dy=1.0_wp/REAL(ny,wp); cell_size=dx
    ALLOCATE(q(nq,nx,ny),qp(nq+2,nx,ny),rate(nx,ny),nodes(nx+1,ny+1),weights(nx+1,ny+1))
    ALLOCATE(B_vertex(nx+1,ny+1),B_cent(nx,ny),B_face_x(nx+1,ny),B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny),B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny),B_second_xy_geom(nx,ny),B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny),grav_coeff_stag_x(nx+1,ny),grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(state%deposit(nx,ny,ns),state%erosion(nx,ny,ns), &
         state%erodible(ns,nx,ny),state%lava_deposit(nx,ny))
    ALLOCATE(state%deposit_tot(nx,ny),state%erosion_tot(nx,ny))
    state%deposit_tot=0.0_wp;state%erosion_tot=0.0_wp
    ALLOCATE(B_zone(nx,ny),cell_source_fractions(nx,ny))
    ALLOCATE(cell_fissure_fractions(nx,ny,n_fissures))
    ALLOCATE(cell_arc_perim(nx,ny),cell_arc_n_x(nx,ny),cell_arc_n_y(nx,ny))
    cell_fissure_fractions=0.0_wp
    cell_arc_perim=0.0_wp;cell_arc_n_x=0.0_wp;cell_arc_n_y=0.0_wp
    CALL domain%initialize
  END SUBROUTINE allocate_grid

  !> \brief Reset a known bed/flow and select a full or wet-plus-halo workset.
  !> \param[in] wet_dry Restrict flow to the left half or southwest quadrant.
  SUBROUTINE reset_state(wet_dry)
    LOGICAL, INTENT(IN) :: wet_dry
    INTEGER :: a,b,l
    REAL(wp) :: x,y,carrier_mass
    state%deposit=0.0_wp; state%erosion=0.0_wp; state%erodible=0.0_wp; state%lava_deposit=0.0_wp
    B_zone=0; cell_source_fractions=0.0_wp
    qp=0.0_wp; qp(4,:,:)=T_ambient
    l=0
    DO b=1,ny
       DO a=1,nx
          carrier_mass=rho_l*0.86_wp
          q(5:6,a,b)=rho_s*[0.10_wp,0.04_wp]
          q(1,a,b)=carrier_mass+SUM(q(5:6,a,b))
          q(2,a,b)=q(1,a,b); q(3,a,b)=-0.3_wp*q(1,a,b)
          q(4,a,b)=350.0_wp*(carrier_mass*sp_heat_l+DOT_PRODUCT(q(5:6,a,b),sp_heat_s))
          IF (wet_dry) THEN
             IF (a>nx/2 .OR. (dimension==2 .AND. b>ny/2)) q(:,a,b)=0.0_wp
             IF (a>nx/2+1 .OR. (dimension==2 .AND. b>ny/2+1)) CYCLE
          END IF
          l=l+1; domain%j_cent(l)=a; domain%k_cent(l)=b
       END DO
    END DO
    domain%solve_cells=l
    DO b=1,ny+1
       y=REAL(b-1,wp)*dy
       DO a=1,nx+1
          x=REAL(a-1,wp)*dx
          B_vertex(a,b)=2.0_wp+0.03_wp*x+0.002_wp*x*x
          IF (dimension==2) B_vertex(a,b)=B_vertex(a,b)-0.02_wp*y+0.001_wp*x*y+0.003_wp*y*y
       END DO
    END DO
    CALL refresh_topography_geometry
    CALL check_geometry
  END SUBROUTINE reset_state

  !> Lie group-0 exchange: verify complete conservative/inventory accounting.
  !> Use the production rate evaluator, executor and geometry update.
  SUBROUTINE check_lie_exchange_transaction
    TYPE(source_registry_type) :: registry
    TYPE(source_group_type) :: groups(1)
    REAL(wp) :: q_before(nq,nx,ny),qp_before(nq+2,nx,ny)
    REAL(wp) :: deposited_before(nx,ny,ns),eroded_before(nx,ny,ns)
    REAL(wp) :: bed_before(nx+1,ny+1),available_before(ns,nx,ny)
    REAL(wp) :: pressure,fluid_delta,solid_delta,inventory_delta
    INTEGER :: a,b,l
    LOGICAL :: failed,accepted
    CHARACTER(LEN=256) :: message,error
    groups(1)%id=0
    groups(1)%method=SOURCE_METHOD_EXPLICIT
    groups(1)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion']
    groups(1)%equation_ids=[(l,l=1,nq)]
    CALL registry%initialize(nq,0,groups,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Lie exchange registry initialization'
    DO b=1,ny
       DO a=1,nx
          CALL qc_to_qp(q(:,a,b),qp(:,a,b),pressure,'Lie exchange test',accepted)
          IF (.NOT.accepted) ERROR STOP 'Lie exchange initial state conversion'
          state%erodible(:,a,b)=[0.002_wp,0.001_wp]
       END DO
    END DO
    q_before=q;qp_before=qp
    deposited_before=state%deposit;eroded_before=state%erosion
    bed_before=B_vertex;available_before=state%erodible
    CALL execute_lie_source_groups(state,q,qp,1.0_wp,domain,registry, &
         q_original=q_before,qp_original=qp_before, &
         step_failed=failed,failure_message=message)
    IF (failed) THEN
       PRINT *, TRIM(message)
       ERROR STOP 'Lie exchange transaction rejected'
    END IF
    fluid_delta=SUM(q(1,:,:)-q_before(1,:,:))
    solid_delta=SUM(q(5:6,:,:)-q_before(5:6,:,:))
    inventory_delta=SUM(state%erosion-eroded_before)
    CALL require_small('Lie solid mass gain / state%erosion inventory', &
         ABS(solid_delta-DOT_PRODUCT(rho_s, &
         [SUM(state%erosion(:,:,1)-eroded_before(:,:,1)), &
          SUM(state%erosion(:,:,2)-eroded_before(:,:,2))])), &
         MAX(1.0_wp,ABS(solid_delta)))
    IF (ANY(state%erodible<0.0_wp)) ERROR STOP 'Lie erodible inventory negative'
    IF (ANY(state%erodible>available_before)) ERROR STOP 'Lie erodible inventory increased'
    DO l=1,ns
       CALL require_small('Lie erodible transfer balance', &
            MAXVAL(ABS(state%erodible(l,:,:)+state%erosion(:,:,l)-eroded_before(:,:,l) &
            -available_before(l,:,:))),1.0_wp)
    END DO
    IF (.NOT.ALL(ieee_is_finite(q)) .OR. .NOT.ALL(ieee_is_finite(B_vertex))) &
         ERROR STOP 'Lie exchange non-finite state'
    IF (fluid_delta<0.0_wp .OR. solid_delta<0.0_wp .OR. inventory_delta<0.0_wp) &
         ERROR STOP 'Lie erosion transfer changed sign'
    CALL check_geometry
    CALL registry%finalize()
    CALL reset_state(.FALSE.)
  END SUBROUTINE check_lie_exchange_transaction

  !> A later-cell failure must not commit any earlier-cell exchange.
  !> The caller uses trial fluid state and the geometry transaction rollback,
  !> matching the Lie time driver. Reject without invoking fatal_error.
  SUBROUTINE check_lie_rejected_transaction
    TYPE(source_registry_type) :: registry
    TYPE(source_group_type) :: groups(1)
    REAL(wp) :: q_before(nq,nx,ny),qp_before(nq+2,nx,ny)
    REAL(wp) :: q_trial(nq,nx,ny),qp_trial(nq+2,nx,ny)
    REAL(wp) :: qp_invalid(nq+2,nx,ny)
    REAL(wp) :: bed_before(nx+1,ny+1)
    REAL(wp) :: deposited_before(nx,ny,ns),eroded_before(nx,ny,ns)
    REAL(wp) :: erodible_before(ns,nx,ny),lava_before(nx,ny)
    REAL(wp) :: centers_before(nx,ny),slope_before(nx,ny)
    REAL(wp) :: pressure
    INTEGER :: a,b,l
    LOGICAL :: failed,accepted
    CHARACTER(LEN=256) :: error,message
    CALL reset_state(.FALSE.)
    groups(1)%id=0
    groups(1)%method=SOURCE_METHOD_EXPLICIT
    groups(1)%source_names=[CHARACTER(LEN=64)::'mass_exchange_erosion']
    groups(1)%equation_ids=[(l,l=1,nq)]
    CALL registry%initialize(nq,0,groups,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Lie rejection registry initialization'
    DO b=1,ny
       DO a=1,nx
          CALL qc_to_qp(q(:,a,b),qp(:,a,b),pressure,'Lie rejection test',accepted)
          IF (.NOT.accepted) ERROR STOP 'Lie rejection primitive conversion'
          state%erodible(:,a,b)=[0.002_wp,0.001_wp]
       END DO
    END DO
    q_before=q;qp_before=qp
    bed_before=B_vertex
    deposited_before=state%deposit;eroded_before=state%erosion
    erodible_before=state%erodible;lava_before=state%lava_deposit
    centers_before=B_cent;slope_before=B_prime_x_geom
    q_trial=q_before;qp_trial=qp_before
    ! Cell 2 is rejected after cell 1 has been considered for transfer.
    ! An independent immutable pre-spatial input avoids aliasing an INTENT(INOUT)
    ! candidate with the executor's INTENT(IN) GROUP-0 reference.
    qp_invalid=qp_before
    qp_invalid(:,2,1)=ieee_value(0.0_wp,ieee_quiet_nan)
    CALL begin_geometry_transaction
    CALL execute_lie_source_groups(state,q_trial,qp_trial,1.0_wp,domain,registry, &
         q_original=q_before,qp_original=qp_invalid, &
         step_failed=failed,failure_message=message)
    IF (.NOT.failed) ERROR STOP 'Lie exchange rejection not reported'
    IF (LEN_TRIM(message)==0) ERROR STOP 'Lie rejection missing diagnostic'
    CALL rollback_geometry_transaction
    ! Fluid workspaces are not published by the Lie driver on rejection.
    q_trial=q_before;qp_trial=qp_before
    IF (ANY(q_trial/=q_before) .OR. ANY(qp_trial/=qp_before)) &
         ERROR STOP 'Rejected Lie fluid state changed'
    IF (ANY(B_vertex/=bed_before) .OR. ANY(B_cent/=centers_before) .OR. &
         ANY(B_prime_x_geom/=slope_before)) ERROR STOP 'Lie bed rollback mismatch'
    IF (ANY(state%deposit/=deposited_before) .OR. ANY(state%erosion/=eroded_before) .OR. &
         ANY(state%erodible/=erodible_before) .OR. ANY(state%lava_deposit/=lava_before)) &
         ERROR STOP 'Lie inventory rollback mismatch'
    CALL registry%finalize()
    CALL reset_state(.FALSE.)
  END SUBROUTINE check_lie_rejected_transaction

  !> \brief Generate reproducible rates with a portable nonoverflowing integer generator.
  !> \param[in] signed_rate Allow either deposition-like or erosion-like bed proposals.
  SUBROUTINE random_rates(signed_rate)
    LOGICAL, INTENT(IN) :: signed_rate
    INTEGER :: a,b
    REAL(wp) :: value
    DO b=1,ny
       DO a=1,nx
          seed=MODULO(16807_int64*seed,2147483647_int64)
          value=REAL(seed,wp)/2147483647.0_wp
          rate(a,b)=-0.001_wp*(0.5_wp+value)
          IF (signed_rate) rate(a,b)=0.002_wp*(value-0.5_wp)
       END DO
    END DO
  END SUBROUTINE random_rates

  !> \brief Independently scatter AC/4 contributions, including all dry cells and boundaries.
  SUBROUTINE reference_nodes
    INTEGER :: a,b
    nodes=0.0_wp; weights=0.0_wp
    DO b=1,ny
       DO a=1,nx
          nodes(a:a+1,b:b+1)=nodes(a:a+1,b:b+1)+0.25_wp*dx*dy*rate(a,b)
          weights(a:a+1,b:b+1)=weights(a:a+1,b:b+1)+0.25_wp*dx*dy
       END DO
    END DO
    nodes=nodes/weights
  END SUBROUTINE reference_nodes

  !> \brief Advance the public signed-rate projector and refresh actual production geometry.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.
  SUBROUTINE projection_update(name,iteration)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: iteration
    REAL(wp) :: actual_nodes(nx+1,ny+1),bed0(nx+1,ny+1),q0(nq,nx,ny)
    CALL reference_nodes
    bed0=B_vertex; q0=q
    CALL project_cell_field_to_vertices(rate,actual_nodes)
    CALL require_small('signed production projection',MAXVAL(ABS(actual_nodes-nodes)),1.0_wp)
    B_vertex=B_vertex+actual_nodes
    CALL refresh_topography_geometry
    CALL require_small('projector does not change flow',MAXVAL(ABS(q-q0)),1.0_wp)
    CALL check_update(name,iteration,bed0)
  END SUBROUTINE projection_update

  !> \brief Apply saturated erosion with an independent complete conservative transaction oracle.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.


  !> \brief Verify bed volume, front localization and all refreshed caches after each update.
  !> \param[in] name Scenario label.
  !> \param[in] iteration Update number.
  !> \param[in] bed0 Authoritative nodal bed before the transaction.
  SUBROUTINE check_update(name,iteration,bed0)
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: iteration
    REAL(wp), INTENT(IN) :: bed0(nx+1,ny+1)
    REAL(wp) :: delta(nx+1,ny+1),geom(nx,ny),mismatch(nx,ny),metrics(6),x,y
    LOGICAL :: front_band
    INTEGER :: a,b
    delta=B_vertex-bed0
    CALL require_small('limited nodal proposal',MAXVAL(ABS(delta-nodes)),MAXVAL(ABS(bed0)))
    CALL require_small('weighted global nodal volume',ABS(SUM(weights*delta)-dx*dy*SUM(rate)),1.0_wp)
    DO b=1,ny
       DO a=1,nx
          geom(a,b)=0.25_wp*SUM(delta(a:a+1,b:b+1))
       END DO
    END DO
    mismatch=geom-rate
    metrics=[dx*dy*SUM(rate),dx*dy*SUM(geom),dx*dy*SUM(mismatch), &
         dx*dy*SUM(ABS(mismatch)),MAXVAL(ABS(mismatch)),0.0_wp]
    CALL require_small('signed global mismatch',ABS(metrics(3)),1.0_wp)
    IF (name=='wet_dry') THEN
       DO b=1,ny
          y=(REAL(b,wp)-0.5_wp)*dy
          DO a=1,nx
             x=(REAL(a,wp)-0.5_wp)*dx
             front_band=ABS(x-0.5_wp)<=dx
             IF (dimension==2) front_band=(front_band .AND. y<=0.5_wp+dy) .OR. &
                  (ABS(y-0.5_wp)<=dy .AND. x<=0.5_wp+dx)
             IF (.NOT.front_band) metrics(6)=MAX(metrics(6),ABS(mismatch(a,b)))
          END DO
       END DO
       CALL require_small('front-localized mismatch',metrics(6),1.0_wp)
       IF (ABS(mismatch(nx/2+1,1))<=tolerance) ERROR STOP 'missing dry geometric redistribution'
    END IF
    IF (metrics(4)<=tolerance) ERROR STOP 'mismatch stress not exercised'
    CALL check_geometry
    WRITE(*,'(A,1X,A,4(1X,I0),6(1X,ES24.16))') 'N8B_RECORD',name,dimension,nx,ny,iteration,metrics
    ! Capture intermediate fields too: final-only equality could hide a transient race.
    WRITE(unit) q,qp,B_vertex,B_cent,B_face_x,B_face_y,B_prime_x_geom,B_prime_y_geom, &
         B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
         grav_coeff_stag_x,grav_coeff_stag_y,state%deposit,state%erosion,state%erodible
  END SUBROUTINE check_update

  !> \brief Independently differentiate Q1 cell means with retained five-point LS kernels.
  SUBROUTINE check_geometry
    REAL(wp), PARAMETER :: c1(5)=[-2.0_wp,-1.0_wp,0.0_wp,1.0_wp,2.0_wp]
    REAL(wp), PARAMETER :: c2(5)=[2.0_wp,-1.0_wp,-2.0_wp,-1.0_wp,2.0_wp]
    REAL(wp) :: centers(nx,ny),bx,by,bxx,byy,bxy,G,bed_scale
    INTEGER :: a,b,ac,bc,s,t
    centers=0.25_wp*(B_vertex(1:nx,1:ny)+B_vertex(2:nx+1,1:ny)+ &
         B_vertex(1:nx,2:ny+1)+B_vertex(2:nx+1,2:ny+1))
    ! Freeze cancellation-aware operation scales before refinement runs:
    ! the epsilon multiplier is unchanged, while derivatives divide bed
    ! elevations by dx/dx**2. Do not scale by a nearly vanishing derivative.
    bed_scale=MAXVAL(ABS(centers))
    CALL require_small('Q1 centers',MAXVAL(ABS(B_cent-centers)),MAXVAL(ABS(centers)))
    CALL require_small('Q1 x faces',MAXVAL(ABS(B_face_x- &
         0.5_wp*(B_vertex(:,1:ny)+B_vertex(:,2:ny+1)))),MAXVAL(ABS(B_vertex)))
    CALL require_small('Q1 y faces',MAXVAL(ABS(B_face_y- &
         0.5_wp*(B_vertex(1:nx,:)+B_vertex(2:nx+1,:)))),MAXVAL(ABS(B_vertex)))
    DO b=1,ny
       bc=1
       IF (dimension==2) bc=MAX(3,MIN(ny-2,b))
       DO a=1,nx
          ac=MAX(3,MIN(nx-2,a))
          bx=DOT_PRODUCT(c1,centers(ac-2:ac+2,bc))/(10.0_wp*dx)
          bxx=DOT_PRODUCT(c2,centers(ac-2:ac+2,bc))/(7.0_wp*dx*dx)
          by=0.0_wp; byy=0.0_wp; bxy=0.0_wp
          IF (dimension==2) THEN
             by=DOT_PRODUCT(c1,centers(ac,bc-2:bc+2))/(10.0_wp*dy)
             byy=DOT_PRODUCT(c2,centers(ac,bc-2:bc+2))/(7.0_wp*dy*dy)
             DO t=1,5
                DO s=1,5
                   bxy=bxy+c1(s)*c1(t)*centers(ac+s-3,bc+t-3)
                END DO
             END DO
             bxy=bxy/(100.0_wp*dx*dy)
          END IF
          G=1.0_wp/(1.0_wp+bx*bx+by*by)
          CALL require_small('filtered x slope',ABS(B_prime_x_geom(a,b)-bx),bed_scale*SUM(ABS(c1))/(10.0_wp*dx))
          CALL require_small('filtered y slope',ABS(B_prime_y_geom(a,b)-by),bed_scale*SUM(ABS(c1))/(10.0_wp*dy))
          CALL require_small('filtered xx Hessian',ABS(B_second_xx_geom(a,b)-bxx),bed_scale*SUM(ABS(c2))/(7.0_wp*dx**2))
          CALL require_small('filtered xy Hessian',ABS(B_second_xy_geom(a,b)-bxy), &
               bed_scale*SUM(ABS(c1))**2/(100.0_wp*dx*dy))
          CALL require_small('filtered yy Hessian',ABS(B_second_yy_geom(a,b)-byy),bed_scale*SUM(ABS(c2))/(7.0_wp*dy**2))
          CALL require_small('center G',ABS(grav_coeff(a,b)-G),1.0_wp)
       END DO
    END DO
    CALL require_small('interior x-face G',MAXVAL(ABS(grav_coeff_stag_x(2:nx,:)- &
         0.5_wp*(grav_coeff(1:nx-1,:)+grav_coeff(2:nx,:)))),1.0_wp)
    IF (ny>1) CALL require_small('interior y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,2:ny)- &
         0.5_wp*(grav_coeff(:,1:ny-1)+grav_coeff(:,2:ny)))),1.0_wp)
    CALL require_small('west G',MAXVAL(ABS(grav_coeff_stag_x(1,:)-grav_coeff(1,:))),1.0_wp)
    CALL require_small('east G',MAXVAL(ABS(grav_coeff_stag_x(nx+1,:)-grav_coeff(nx,:))),1.0_wp)
    CALL require_small('south G',MAXVAL(ABS(grav_coeff_stag_y(:,1)-grav_coeff(:,1))),1.0_wp)
    CALL require_small('north G',MAXVAL(ABS(grav_coeff_stag_y(:,ny+1)-grav_coeff(:,ny))),1.0_wp)
  END SUBROUTINE check_geometry

  !> \brief Reject nonfinite or excessive errors with the unchanged N8 roundoff threshold.
  !> \param[in] label Assertion label.
  !> \param[in] error Absolute error.
  !> \param[in] scale Reference magnitude.
  SUBROUTINE require_small(label,error,scale)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: error,scale
    IF (.NOT.ieee_is_finite(error) .OR. error>tolerance*MAX(1.0_wp,scale)) THEN
       WRITE(*,*) 'FAIL: ',label,error,tolerance*MAX(1.0_wp,scale)
       ERROR STOP 1
    END IF
  END SUBROUTINE require_small

  !> \brief Release test-grid allocations before changing resolution or dimension.
  SUBROUTINE release_grid
    CALL domain%finalize
    DEALLOCATE(q,qp,rate,nodes,weights,B_vertex,B_cent,B_face_x,B_face_y)
    DEALLOCATE(B_prime_x_geom,B_prime_y_geom,B_second_xx_geom,B_second_xy_geom,B_second_yy_geom)
    DEALLOCATE(grav_coeff,grav_coeff_stag_x,grav_coeff_stag_y,B_zone,cell_source_fractions)
    DEALLOCATE(cell_fissure_fractions,cell_arc_perim,cell_arc_n_x,cell_arc_n_y)
    DEALLOCATE(state%deposit,state%erosion,state%erodible,state%lava_deposit, &
         state%deposit_tot,state%erosion_tot)
  END SUBROUTINE release_grid
END PROGRAM test_projection_stress
