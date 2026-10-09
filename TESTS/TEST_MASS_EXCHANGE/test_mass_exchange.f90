!> \brief N8 production exchange budgets, limited nodal updates and true checkpoint round trips.
PROGRAM test_mass_exchange

  USE parameters_2d
  USE constitutive_parameters_2d
  USE geometry_2d
  USE domain_2d, ONLY : domain_type
  USE state_2d, ONLY : state_type
  USE model_layout_2d, ONLY : model_layout_type
  USE runtime_2d, ONLY : runtime_state_type
  USE stochastic_module, ONLY : stochastic_workspace_type, sym_noise, noise_pow_val
  USE state_conversion_2d, ONLY : qc_to_qp, settling_velocity
  USE mass_exchange_2d, ONLY : update_erosion_deposition_cell, release_topography_workspace
  USE inpout_2d, ONLY : write_restart_file, read_restart_file, output_idx, &
       t_output, t_runout, t_probes
  USE omp_lib
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

  IMPLICIT NONE

  INTEGER, PARAMETER :: nx=8, ny=7, ns=2, nq=6
  REAL(wp), PARAMETER :: tol=2048.0_wp*EPSILON(1.0_wp)
  TYPE(state_type) :: state
  TYPE(domain_type) :: domain
  TYPE(model_layout_type) :: layout
  TYPE(runtime_state_type) :: runtime
  TYPE(stochastic_workspace_type) :: stochastic
  CHARACTER(LEN=64) :: selected, argument
  INTEGER :: actual_threads, requested_threads, status, snapshot_unit

  CALL get_command_argument(1, argument)
  READ(argument,*) requested_threads
  actual_threads = 0
  !$OMP PARALLEL
  !$OMP SINGLE
  actual_threads = omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual_threads /= requested_threads) ERROR STOP 'wrong actual OpenMP team'
  WRITE(*,*) 'Actual OpenMP team:', actual_threads
  CALL get_environment_variable('N8_CASE', selected, STATUS=status)
  IF (status /= 0) selected = 'all'
  CALL initialize_fixture
  OPEN(NEWUNIT=snapshot_unit, FILE='snapshot.bin', ACCESS='STREAM',      &
       FORM='UNFORMATTED', STATUS='REPLACE')

  SELECT CASE(TRIM(selected))
  CASE('all')
     CALL run_case('deposition')
     CALL run_case('erosion')
     CALL run_case('combined')
     CALL run_case('masked')
     CALL run_case('masked_fissural')
     CALL run_case('cutoff')
     CALL run_case('loss_only')
     CALL run_case('loss_reserve')
     CALL run_case('thermal')
     CALL run_restart
  CASE('deposition','erosion','combined','masked','masked_fissural','cutoff','loss_only','loss_reserve','thermal')
     CALL run_case(TRIM(selected))
  CASE('restart')
     CALL run_restart
  CASE DEFAULT
     ERROR STOP 'unknown N8_CASE'
  END SELECT

  CLOSE(snapshot_unit)
  CALL release_topography_workspace
  CALL release_topography_workspace
  CALL stochastic%finalize
  CALL domain%finalize
  CALL state%finalize
  CALL layout%finalize
  WRITE(*,*) 'PASS: production inventory limits, bed volume, geometry and restart'

CONTAINS

  !> \brief Allocate a rectangular two-solid liquid fixture and the real restart dependencies.
  SUBROUTINE initialize_fixture
    INTEGER :: j,k,l
    comp_cells_x=nx; comp_cells_y=ny
    comp_interfaces_x=nx+1; comp_interfaces_y=ny+1
    dx=0.75_wp; dy=1.25_wp; cell_size=dx
    n_solid=ns; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
    n_vars=nq; n_eqns=nq; n_RK=2
    n_thickness_levels=1; n_dyn_pres_levels=1
    idx_h=1; idx_hu=2; idx_hv=3; idx_T=4
    idx_solid_first=5; idx_solid_last=6
    idx_solidEqn_first=5; idx_solidEqn_last=6
    idx_add_gas_first=7; idx_add_gas_last=6
    idx_stoch=0; idx_pore=0; idx_u=7; idx_v=8
    gas_flag=.FALSE.; liquid_flag=.TRUE.; rheology_flag=.FALSE.
    rheology_model=0; stochastic_flag=.FALSE.; stoch_transport_flag=.FALSE.
    pore_pressure_flag=.FALSE.; gas_loss_flag=.FALSE.
    liquid_vaporization_flag=.FALSE.; sutherland_flag=.FALSE.
    entrainment_flag=.FALSE.; slope_correction_flag=.TRUE.
    topo_change_flag=.TRUE.; radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    grav=9.81_wp; T_ambient=300.0_wp; rho_a_amb=1.2_wp
    rho_l=1000.0_wp; inv_rho_l=1.0_wp/rho_l; sp_heat_l=4180.0_wp
    sp_heat_a=998.0_wp; sp_gas_const_a=287.05_wp
    pres=101300.0_wp; inv_pres=1.0_wp/pres
    rho_c_sub=rho_l; kin_visc_c=1.0E-6_wp; hydraulic_permeability=0.0_wp
    eps_sing=1.0E-8_wp; eps_sing4=eps_sing**4
    maximum_solid_packing=0.6_wp
    ALLOCATE(rho_s(ns), inv_rho_s(ns), sp_heat_s(ns), diam_s(ns), erodible_fract(ns))
    ALLOCATE(sp_heat_g(0), sp_gas_const_g(0), loss_rate)
    rho_s=[2000.0_wp,3000.0_wp]; inv_rho_s=1.0_wp/rho_s
    sp_heat_s=[900.0_wp,1100.0_wp]; diam_s=[1.0E-4_wp,2.0E-4_wp]
    erodible_fract=[0.6_wp,0.4_wp]
    ALLOCATE(B_vertex(nx+1,ny+1), B_cent(nx,ny), B_face_x(nx+1,ny), B_face_y(nx,ny+1))
    ALLOCATE(B_prime_x_geom(nx,ny), B_prime_y_geom(nx,ny))
    ALLOCATE(B_second_xx_geom(nx,ny), B_second_xy_geom(nx,ny), B_second_yy_geom(nx,ny))
    ALLOCATE(grav_coeff(nx,ny), grav_coeff_stag_x(nx+1,ny), grav_coeff_stag_y(nx,ny+1))
    ALLOCATE(deposit(nx,ny,ns), erosion(nx,ny,ns), erodible(ns,nx,ny))
    ALLOCATE(B_zone(nx,ny), cell_source_fractions(nx,ny))
    CALL layout%initialize(1,nq,nq,ns,0,0,0,.FALSE.)
    CALL state%initialize(layout)
    CALL domain%initialize
    ! All cells except one wet inactive corner are on the production workset.
    ! Two dry cells remain on the list as a realistic solve/halo boundary.
    l=0
    DO k=1,ny
       DO j=1,nx
          IF (j==nx .AND. k==ny) CYCLE
          l=l+1; domain%j_cent(l)=j; domain%k_cent(l)=k
       END DO
    END DO
    domain%solve_cells=l
    sym_noise=0.0_wp; noise_pow_val=1.0_wp
    CALL stochastic%initialize
    output_idx=0; t_output=0.0_wp; t_runout=0.0_wp; t_probes=0.0_wp
  END SUBROUTINE initialize_fixture

  !> \brief Reset conservative inventories and an analytic initial bed before each independent case.
  SUBROUTINE reset_fixture
    INTEGER :: j,k
    REAL(wp) :: h, phi(ns), temperature, x,y
    erosion_coeff=0.0_wp; settling_flag=.FALSE.; loss_rate=0.0_wp
    alphastot_min=0.0_wp; erodible_porosity=0.25_wp
    coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
    erodible_deposit_flag=.FALSE.; bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.
    T_erodible=317.0_wp; verbose_level=2
    deposit=0.0_wp; erosion=0.0_wp; erodible=0.0_wp
    B_zone=0; cell_source_fractions=0.0_wp
    DO k=1,ny+1
       y=REAL(k-1,wp)*dy
       DO j=1,nx+1
          x=REAL(j-1,wp)*dx
          B_vertex(j,k)=2.0_wp+0.03_wp*x-0.02_wp*y+0.002_wp*x*x+0.001_wp*x*y
       END DO
    END DO
    CALL refresh_topography_geometry
    DO k=1,ny
       DO j=1,nx
          h=0.5_wp+0.02_wp*REAL(j+k,wp)
          phi=[0.10_wp,0.04_wp]
          temperature=300.0_wp+2.0_wp*REAL(j,wp)+3.0_wp*REAL(k,wp)
          CALL set_cell(j,k,h,phi,temperature)
          erodible(:,j,k)=[0.002_wp,0.001_wp]*(1.0_wp+0.01_wp*REAL(j+k,wp))
       END DO
    END DO
    state%q(:,5,4)=0.0_wp; state%q(:,7,6)=0.0_wp
  END SUBROUTINE reset_fixture

  !> \brief Construct conserved masses directly from known phase volumes, avoiding inverse test logic.
  !> \param[in] j,k Cell indices.
  !> \param[in] h Thickness [m].
  !> \param[in] phi Solid volume fractions.
  !> \param[in] temperature Flow temperature [K].
  SUBROUTINE set_cell(j,k,h,phi,temperature)
    INTEGER, INTENT(IN) :: j,k
    REAL(wp), INTENT(IN) :: h, phi(ns), temperature
    REAL(wp) :: liquid_mass
    liquid_mass=rho_l*h*(1.0_wp-SUM(phi))
    state%q(5:6,j,k)=rho_s*h*phi
    state%q(1,j,k)=liquid_mass+SUM(state%q(5:6,j,k))
    state%q(2,j,k)=state%q(1,j,k)*(1.0_wp+0.01_wp*REAL(j,wp))
    state%q(3,j,k)=state%q(1,j,k)*(-0.3_wp+0.01_wp*REAL(k,wp))
    state%q(4,j,k)=temperature*(liquid_mass*sp_heat_l+DOT_PRODUCT(state%q(5:6,j,k),sp_heat_s))
  END SUBROUTINE set_cell

  !> \brief Select a saturation/guard case and apply the real production transaction.
  !> \param[in] name Case name, also used by the diagnostic comparator.
  SUBROUTINE run_case(name)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp) :: dt
    INTEGER :: j,k
    CALL reset_fixture
    dt=1.0_wp
    SELECT CASE(name)
    CASE('deposition')
       settling_flag=.TRUE.; dt=1.0E6_wp
    CASE('erosion')
       erosion_coeff=100.0_wp
    CASE('combined','masked','masked_fissural')
       erosion_coeff=100.0_wp; settling_flag=.TRUE.
       erodible_deposit_flag=.TRUE.; dt=1.0E6_wp
       IF (name=='masked' .OR. name=='masked_fissural') THEN
          bottom_radial_source_flag=name=='masked'
          bottom_fissural_source_flag=name=='masked_fissural'
          cell_source_fractions(1,1)=1.0_wp
          cell_source_fractions(2,1)=0.5_wp
          cell_source_fractions(3,1)=0.25_wp
       END IF
    CASE('cutoff')
       erosion_coeff=100.0_wp; alphastot_min=0.3_wp
    CASE('loss_only')
       loss_rate=1000.0_wp
    CASE('loss_reserve')
       loss_rate=1000.0_wp
       ! Activate the exchange evaluator through erosion as well, with an
       ! empty substrate: the carrier reserve must never become a gain term.
       erosion_coeff=1.0_wp; erodible=0.0_wp
       erodible_porosity=0.7_wp
       coeff_porosity=erodible_porosity/(1.0_wp-erodible_porosity)
       DO k=1,ny
          DO j=1,nx
             IF (state%q(1,j,k)>0.0_wp) CALL set_cell(j,k,0.8_wp,[0.38_wp,0.17_wp],350.0_wp)
          END DO
       END DO
    CASE('thermal')
       erosion_coeff=100.0_wp; erodible_deposit_flag=.TRUE.
    END SELECT
    CALL check_update(name,dt)
  END SUBROUTINE run_case

  !> \brief Compute independently limited exchange increments for the liquid two-solid fixture.
  !> \param[in] q0 Conservative base state.
  !> \param[in] inventory Available solid volume per area [m].
  !> \param[in] bx,by Filtered slopes used in the existing erosion law.
  !> \param[in] dt Step duration [s].
  !> \param[out] dep,ers Limited solid volume increments [m].
  !> \param[out] loss Limited carrier volume loss [m].
  !> \param[out] h,u,v,temperature Physical base state recovered analytically.
  SUBROUTINE reference_exchange(q0,inventory,bx,by,dt,dep,ers,loss,h,u,v,temperature)
    REAL(wp), INTENT(IN) :: q0(nq),inventory(ns),bx,by,dt
    REAL(wp), INTENT(OUT) :: dep(ns),ers(ns),loss,h,u,v,temperature
    REAL(wp) :: solid_volume(ns), phi(ns), carrier_mass, speed, settling, reserve
    INTEGER :: i
    dep=0.0_wp; ers=0.0_wp; loss=0.0_wp
    h=0.0_wp; u=0.0_wp; v=0.0_wp; temperature=T_ambient
    IF (q0(1)<=0.0_wp) RETURN
    solid_volume=q0(5:6)/rho_s
    carrier_mass=q0(1)-SUM(q0(5:6))
    h=SUM(solid_volume)+carrier_mass/rho_l
    phi=solid_volume/h
    u=q0(2)/q0(1); v=q0(3)/q0(1)
    temperature=q0(4)/(carrier_mass*sp_heat_l+DOT_PRODUCT(q0(5:6),sp_heat_s))
    IF (SUM(phi)<=alphastot_min) RETURN
    speed=SQRT(u*u+v*v+(u*bx+v*by)**2)
    ers=MAX(0.0_wp,MIN(inventory,dt*erosion_coeff*speed*h*(1.0_wp-SUM(phi)) &
         *(1.0_wp-erodible_porosity)*erodible_fract))
    IF (settling_flag) THEN
       DO i=1,ns
          settling=settling_velocity(diam_s(i),rho_s(i),rho_l,1.0_wp/kin_visc_c)
          dep(i)=MIN(solid_volume(i),dt*phi(i)*settling                    &
               *(1.0_wp-MIN(1.0_wp,SUM(phi)/0.6_wp))**4.65_wp)
       END DO
    END IF
    reserve=carrier_mass/rho_l-coeff_porosity*(SUM(solid_volume)-SUM(dep))
    loss=MAX(0.0_wp,MIN(dt*loss_rate+coeff_porosity*SUM(dep),              &
         h*MAX(0.0_wp,maximum_solid_packing-SUM(phi)),reserve))
  END SUBROUTINE reference_exchange

  !> \brief Assert the complete production transaction and projection after all cell-local limits.
  !> \param[in] name Case label for failures and diagnostic output.
  !> \param[in] dt Duration of the production update [s].
  SUBROUTINE check_update(name,dt)
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: q0(nq,nx,ny), dep0(nx,ny,ns), ers0(nx,ny,ns), inv0(ns,nx,ny)
    REAL(wp) :: bed0(nx+1,ny+1), q_expected(nq,nx,ny), inv_expected(ns,nx,ny)
    REAL(wp) :: dep_expected(nx,ny,ns), ers_expected(nx,ny,ns), cell_delta(nx,ny)
    REAL(wp) :: vertex_delta(nx+1,ny+1), geom_delta(nx,ny), mismatch(nx,ny)
    REAL(wp) :: dep(ns),ers(ns),loss,h,u,v,temperature,substrate_temperature,mask
    REAL(wp) :: outgoing_mass, entering_mass, weight, nodal_volume, diagnostics(5)
    INTEGER :: j,k,l,j0,j1,k0,k1
    q0=state%q; bed0=B_vertex; dep0=deposit; ers0=erosion; inv0=erodible
    q_expected=q0; dep_expected=dep0; ers_expected=ers0; inv_expected=inv0
    cell_delta=0.0_wp
    DO l=1,domain%solve_cells
       j=domain%j_cent(l); k=domain%k_cent(l)
       CALL reference_exchange(q0(:,j,k),inv0(:,j,k),B_prime_x_geom(j,k), &
            B_prime_y_geom(j,k),dt,dep,ers,loss,h,u,v,temperature)
       IF (j==1 .AND. k==1) THEN
          SELECT CASE(name)
          CASE('deposition','combined','masked','masked_fissural')
             CALL assert_small(name//' saturated deposition',MAXVAL(ABS(dep-q0(5:6,j,k)/rho_s)),1.0_wp)
          END SELECT
          SELECT CASE(name)
          CASE('erosion','combined','masked','masked_fissural','thermal')
             CALL assert_small(name//' saturated erosion',MAXVAL(ABS(ers-inv0(:,j,k))),1.0_wp)
          END SELECT
          IF (name=='loss_only' .AND. loss<=tol) ERROR STOP 'carrier loss not exercised'
          IF (name=='loss_reserve' .AND. loss/=0.0_wp) ERROR STOP 'carrier reserve not exhausted'
       END IF
       mask=1.0_wp
       IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) &
            mask=1.0_wp-cell_source_fractions(j,k)
       dep=dep*mask; ers=ers*mask; loss=loss*mask
       dep_expected(j,k,:)=dep0(j,k,:)+dep
       ers_expected(j,k,:)=ers0(j,k,:)+ers
       inv_expected(:,j,k)=inv0(:,j,k)-ers
       IF (erodible_deposit_flag) inv_expected(:,j,k)=inv_expected(:,j,k)+dep
       outgoing_mass=DOT_PRODUCT(rho_s,dep)+rho_l*loss
       entering_mass=DOT_PRODUCT(rho_s,ers)+rho_c_sub*coeff_porosity*SUM(ers)
       q_expected(1,j,k)=q0(1,j,k)+entering_mass-outgoing_mass
       q_expected(2:3,j,k)=q0(2:3,j,k)-outgoing_mass*[u,v]
       substrate_temperature=T_erodible
       IF (erodible_deposit_flag) substrate_temperature=temperature
       q_expected(4,j,k)=q0(4,j,k)                                      &
            -temperature*(DOT_PRODUCT(rho_s*sp_heat_s,dep)+rho_l*sp_heat_l*loss) &
            +substrate_temperature*(DOT_PRODUCT(rho_s*sp_heat_s,ers)     &
            +rho_c_sub*sp_heat_l*coeff_porosity*SUM(ers))
       q_expected(5:6,j,k)=q0(5:6,j,k)+rho_s*(ers-dep)
       cell_delta(j,k)=SUM(dep-ers)/(1.0_wp-erodible_porosity)
    END DO
    DO k=1,ny+1
       DO j=1,nx+1
          j0=MAX(1,j-1); j1=MIN(nx,j); k0=MAX(1,k-1); k1=MIN(ny,k)
          vertex_delta(j,k)=SUM(cell_delta(j0:j1,k0:k1))/REAL((j1-j0+1)*(k1-k0+1),wp)
       END DO
    END DO
    CALL update_erosion_deposition_cell(state%q,state%qp,dt,domain)
    CALL assert_small(name//' conservative transaction',MAXVAL(ABS(state%q-q_expected)),MAXVAL(ABS(q0)))
    CALL assert_small(name//' deposit',MAXVAL(ABS(deposit-dep_expected)),1.0_wp)
    CALL assert_small(name//' erosion',MAXVAL(ABS(erosion-ers_expected)),1.0_wp)
    CALL assert_small(name//' erodible',MAXVAL(ABS(erodible-inv_expected)),1.0_wp)
    CALL assert_small(name//' read-only substrate temperature',ABS(T_erodible-317.0_wp),317.0_wp)
    IF (MINVAL(erodible)<-tol .OR. MINVAL(state%q(5:6,:,:))<-tol) ERROR STOP 'negative inventory'
    CALL assert_small(name//' limited nodal proposal',MAXVAL(ABS(B_vertex-bed0-vertex_delta)),MAXVAL(ABS(bed0)))
    ! Sum the solid flow budget and the actual accumulated inventory exchange.
    DO j=1,ns
       CALL assert_small(name//' solid budget',                          &
            ABS(SUM(state%q(4+j,:,:)-q0(4+j,:,:))+rho_s(j)*              &
            SUM(deposit(:,:,j)-dep0(:,:,j)-erosion(:,:,j)+ers0(:,:,j))), &
            SUM(ABS(q0(4+j,:,:))))
    END DO
    nodal_volume=0.0_wp
    DO k=1,ny+1
       DO j=1,nx+1
          weight=dx*dy
          IF (j==1 .OR. j==nx+1) weight=0.5_wp*weight
          IF (k==1 .OR. k==ny+1) weight=0.5_wp*weight
          nodal_volume=nodal_volume+weight*(B_vertex(j,k)-bed0(j,k))
       END DO
    END DO
    CALL assert_small(name//' Q1 weighted volume',ABS(nodal_volume-dx*dy*SUM(cell_delta)), &
         MAX(1.0_wp,dx*dy*SUM(ABS(cell_delta))))
    DO k=1,ny
       DO j=1,nx
          geom_delta(j,k)=0.25_wp*SUM(vertex_delta(j:j+1,k:k+1))
       END DO
    END DO
    mismatch=geom_delta-cell_delta
    IF (SUM(ABS(cell_delta))>tol) THEN
       IF (ABS(mismatch(5,4))<=tol .OR. ABS(mismatch(nx,ny))<=tol) &
            ERROR STOP 'missing dry/inactive geometric mismatch'
    END IF
    diagnostics=[dx*dy*SUM(cell_delta),dx*dy*SUM(geom_delta),dx*dy*SUM(mismatch), &
         dx*dy*SUM(ABS(mismatch)),MAXVAL(ABS(mismatch))]/dt
    WRITE(*,'(A,1X,A,5(1X,ES24.16))') 'N8_EXPECTED_DIAGNOSTICS',name,diagnostics
    CALL check_geometry
    WRITE(snapshot_unit) state%q,state%qp,B_vertex,B_cent,B_face_x,B_face_y, &
         B_prime_x_geom,B_prime_y_geom,B_second_xx_geom,B_second_xy_geom,   &
         B_second_yy_geom,grav_coeff,grav_coeff_stag_x,grav_coeff_stag_y,   &
         deposit,erosion,erodible
    WRITE(*,*) 'PASS: ',name
  END SUBROUTINE check_update

  !> \brief Independently evaluate Q1 identities and the retained five-point LS derivative filter.
  SUBROUTINE check_geometry
    REAL(wp), PARAMETER :: c1(5)=[-2.0_wp,-1.0_wp,0.0_wp,1.0_wp,2.0_wp]
    REAL(wp), PARAMETER :: c2(5)=[2.0_wp,-1.0_wp,-2.0_wp,-1.0_wp,2.0_wp]
    REAL(wp) :: centers(nx,ny), bx,by,bxx,byy,bxy,G
    INTEGER :: j,k,jc,kc,a,b
    DO k=1,ny
       DO j=1,nx
          centers(j,k)=0.25_wp*SUM(B_vertex(j:j+1,k:k+1))
       END DO
    END DO
    CALL assert_small('fresh Q1 centers',MAXVAL(ABS(B_cent-centers)),MAXVAL(ABS(centers)))
    CALL assert_small('fresh Q1 x faces',MAXVAL(ABS(B_face_x-             &
         0.5_wp*(B_vertex(:,1:ny)+B_vertex(:,2:ny+1)))),MAXVAL(ABS(B_vertex)))
    CALL assert_small('fresh Q1 y faces',MAXVAL(ABS(B_face_y-             &
         0.5_wp*(B_vertex(1:nx,:)+B_vertex(2:nx+1,:)))),MAXVAL(ABS(B_vertex)))
    DO k=1,ny
       kc=MAX(3,MIN(ny-2,k))
       DO j=1,nx
          jc=MAX(3,MIN(nx-2,j))
          bx=DOT_PRODUCT(c1,centers(jc-2:jc+2,kc))/(10.0_wp*dx)
          by=DOT_PRODUCT(c1,centers(jc,kc-2:kc+2))/(10.0_wp*dy)
          bxx=DOT_PRODUCT(c2,centers(jc-2:jc+2,kc))/(7.0_wp*dx*dx)
          byy=DOT_PRODUCT(c2,centers(jc,kc-2:kc+2))/(7.0_wp*dy*dy)
          bxy=0.0_wp
          DO b=1,5
             DO a=1,5
                bxy=bxy+c1(a)*c1(b)*centers(jc+a-3,kc+b-3)
             END DO
          END DO
          bxy=bxy/(100.0_wp*dx*dy)
          G=1.0_wp/(1.0_wp+bx*bx+by*by)
          CALL assert_small('fresh slope x',ABS(B_prime_x_geom(j,k)-bx),1.0_wp)
          CALL assert_small('fresh slope y',ABS(B_prime_y_geom(j,k)-by),1.0_wp)
          CALL assert_small('fresh curvature xx',ABS(B_second_xx_geom(j,k)-bxx),1.0_wp)
          CALL assert_small('fresh curvature yy',ABS(B_second_yy_geom(j,k)-byy),1.0_wp)
          CALL assert_small('fresh curvature xy',ABS(B_second_xy_geom(j,k)-bxy),1.0_wp)
          CALL assert_small('fresh center G',ABS(grav_coeff(j,k)-G),1.0_wp)
       END DO
    END DO
    CALL assert_small('fresh interior x-face G',MAXVAL(ABS(grav_coeff_stag_x(2:nx,:)- &
         0.5_wp*(grav_coeff(1:nx-1,:)+grav_coeff(2:nx,:)))),1.0_wp)
    CALL assert_small('fresh interior y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,2:ny)- &
         0.5_wp*(grav_coeff(:,1:ny-1)+grav_coeff(:,2:ny)))),1.0_wp)
    CALL assert_small('fresh boundary x-face G',MAXVAL(ABS(grav_coeff_stag_x(1,:)-grav_coeff(1,:))),1.0_wp)
    CALL assert_small('fresh boundary y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,1)-grav_coeff(:,1))),1.0_wp)
    CALL assert_small('fresh east x-face G',MAXVAL(ABS(grav_coeff_stag_x(nx+1,:)-grav_coeff(nx,:))),1.0_wp)
    CALL assert_small('fresh north y-face G',MAXVAL(ABS(grav_coeff_stag_y(:,ny+1)-grav_coeff(:,ny))),1.0_wp)
  END SUBROUTINE check_geometry

  !> \brief Restore a real checkpoint after corrupting caches, then compare a nonzero second update.
  SUBROUTINE run_restart
    REAL(wp) :: q1(nq,nx,ny),bed1(nx+1,ny+1),dep1(nx,ny,ns),ers1(nx,ny,ns),inv1(ns,nx,ny)
    REAL(wp) :: q2(nq,nx,ny),bed2(nx+1,ny+1),dep2(nx,ny,ns),ers2(nx,ny,ns),inv2(ns,nx,ny)
    CALL reset_fixture
    settling_flag=.TRUE.; erosion_coeff=0.001_wp; erodible_deposit_flag=.TRUE.
    CALL check_update('restart_first',0.2_wp)
    q1=state%q; bed1=B_vertex; dep1=deposit; ers1=erosion; inv1=erodible
    IF (SUM(deposit)+SUM(erosion)<=0.0_wp) ERROR STOP 'restart bed did not evolve'
    runtime%t=0.2_wp; runtime%dt=0.2_wp
    CALL write_restart_file('checkpoint.bin',runtime,stochastic,state,domain)
    CALL check_update('restart_continuous',0.2_wp)
    q2=state%q; bed2=B_vertex; dep2=deposit; ers2=erosion; inv2=erodible
    IF (ALL(bed2==bed1)) ERROR STOP 'restart continuation is a no-op'
    state%q=-999.0_wp; state%qp=-999.0_wp; B_vertex=-999.0_wp
    deposit=-999.0_wp; erosion=-999.0_wp; erodible=-999.0_wp
    B_cent=-999.0_wp; B_face_x=-999.0_wp; B_face_y=-999.0_wp
    B_prime_x_geom=-999.0_wp; B_prime_y_geom=-999.0_wp
    B_second_xx_geom=-999.0_wp; B_second_xy_geom=-999.0_wp; B_second_yy_geom=-999.0_wp
    grav_coeff=-999.0_wp; grav_coeff_stag_x=-999.0_wp; grav_coeff_stag_y=-999.0_wp
    CALL read_restart_file('checkpoint.bin',runtime,stochastic,state,domain)
    CALL same_bits('restart q',RESHAPE(state%q,[SIZE(q1)]),RESHAPE(q1,[SIZE(q1)]))
    CALL same_bits('restart B_vertex',RESHAPE(B_vertex,[SIZE(bed1)]),RESHAPE(bed1,[SIZE(bed1)]))
    CALL same_bits('restart deposit',RESHAPE(deposit,[SIZE(dep1)]),RESHAPE(dep1,[SIZE(dep1)]))
    CALL same_bits('restart erosion',RESHAPE(erosion,[SIZE(ers1)]),RESHAPE(ers1,[SIZE(ers1)]))
    CALL same_bits('restart erodible',RESHAPE(erodible,[SIZE(inv1)]),RESHAPE(inv1,[SIZE(inv1)]))
    CALL check_geometry
    ! qp is intentionally not checkpointed; the main program normally rebuilds
    ! it from q. Reset its inactive cache before the continuation fingerprint.
    state%qp=0.0_wp; state%qp(4,:,:)=T_ambient
    CALL check_update('restart_resumed',0.2_wp)
    CALL same_bits('continued q',RESHAPE(state%q,[SIZE(q2)]),RESHAPE(q2,[SIZE(q2)]))
    CALL same_bits('continued bed',RESHAPE(B_vertex,[SIZE(bed2)]),RESHAPE(bed2,[SIZE(bed2)]))
    CALL same_bits('continued deposit',RESHAPE(deposit,[SIZE(dep2)]),RESHAPE(dep2,[SIZE(dep2)]))
    CALL same_bits('continued erosion',RESHAPE(erosion,[SIZE(ers2)]),RESHAPE(ers2,[SIZE(ers2)]))
    CALL same_bits('continued inventory',RESHAPE(erodible,[SIZE(inv2)]),RESHAPE(inv2,[SIZE(inv2)]))
  END SUBROUTINE run_restart

  !> \brief Reject nonfinite values and scaled errors beyond the specified roundoff tolerance.
  !> \param[in] label Assertion name.
  !> \param[in] error Absolute error.
  !> \param[in] scale Reference magnitude for the tolerance.
  SUBROUTINE assert_small(label,error,scale)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: error,scale
    IF (.NOT.ieee_is_finite(error) .OR. error>tol*MAX(1.0_wp,scale)) THEN
       WRITE(*,*) 'FAIL: ',label,error,tol*MAX(1.0_wp,scale)
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_small

  !> \brief Compare serialized REAL fields including their exact floating-point bit patterns.
  !> \param[in] label Assertion name.
  !> \param[in] actual,reference Flattened REAL fields to compare.
  SUBROUTINE same_bits(label,actual,reference)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: actual(:),reference(:)
    IF (ANY(TRANSFER(actual,[0_int64],SIZE(actual)) /=                   &
         TRANSFER(reference,[0_int64],SIZE(reference)))) THEN
       WRITE(*,*) 'FAIL: ',label
       ERROR STOP 1
    END IF
  END SUBROUTINE same_bits

END PROGRAM test_mass_exchange
