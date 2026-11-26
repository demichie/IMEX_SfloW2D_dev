!********************************************************************************
!> \brief Hyperbolic flux evaluation
!>
!> This module owns characteristic speeds and semidiscrete numerical interface
!> fluxes. Active-cell and active-interface lists are supplied explicitly by
!> the solver driver.
!>
!> Owns one wave-speed bound per face and oriented HP-PCCU face values. Active-cell and face index
!> lists restrict evaluation to flow plus the required reconstruction halo.
!********************************************************************************

MODULE hyperbolic_2d

    USE parameters_2d, ONLY : wp, n_eqns, n_vars, idx_u, idx_v
    USE parameters_2d, ONLY : dry_thickness_tolerance
    USE parameters_2d, ONLY : idx_solidEqn_first, idx_solidEqn_last
    USE parameters_2d, ONLY : idx_solid_first, idx_solid_last
    USE parameters_2d, ONLY : kt_trace_reconstruction
    USE parameters_2d, ONLY : lateral_source_flag
    USE parameters_2d, ONLY : bcW, bcE, bcS, bcN
    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
    USE geometry_2d, ONLY : grav_coeff_stag_x, grav_coeff_stag_y
    USE geometry_2d, ONLY : B_face_x, B_face_y
    USE geometry_2d, ONLY : one_by_dx, one_by_dy
    USE geometry_2d, ONLY : source_cell

    USE reconstruction_2d, ONLY : reconstruction_workspace_type
    IMPLICIT NONE

    PRIVATE

    PUBLIC :: add_interlayer_terms, eval_total_depth_face_speed_bounds

    !> \brief Persistent wave speeds, oriented face contributions and cell paths.
    !> \details Face arrays distinguish the contribution seen by either adjacent
    !>          cell; hydrostatic pressure belongs to the path terms, not the
    !>          inertial endpoint flux.
    TYPE, PUBLIC :: hyperbolic_workspace_type
        ! One common characteristic bound per face and sign, not per equation.
        !> Negative-going x-normal wave-speed bound at each x face [m/s].
        REAL(wp), ALLOCATABLE :: a_interface_xNeg(:,:)
        !> Positive-going x-normal wave-speed bound at each x face [m/s].
        REAL(wp), ALLOCATABLE :: a_interface_xPos(:,:)
        !> Negative-going y-normal wave-speed bound at each y face [m/s].
        REAL(wp), ALLOCATABLE :: a_interface_yNeg(:,:)
        !> Positive-going y-normal wave-speed bound at each y face [m/s].
        REAL(wp), ALLOCATABLE :: a_interface_yPos(:,:)
        ! Hydrostatic coefficients derived while evaluating the same final
        ! endpoint traces for characteristic speeds. Owned by the HP workspace.
        !> Left/right hydrostatic density-gravity coefficients at x faces.
        REAL(wp), ALLOCATABLE :: gamma_interface_xL(:,:),gamma_interface_xR(:,:)
        !> Bottom/top hydrostatic density-gravity coefficients at y faces.
        REAL(wp), ALLOCATABLE :: gamma_interface_yB(:,:),gamma_interface_yT(:,:)
        !> Hydrostatic path contributions on left-oriented x-face traces.
        REAL(wp), ALLOCATABLE :: G_interface_xL(:,:,:)
        !> Hydrostatic path contributions on right-oriented x-face traces.
        REAL(wp), ALLOCATABLE :: G_interface_xR(:,:,:)
        !> Hydrostatic path contributions on bottom-oriented y-face traces.
        REAL(wp), ALLOCATABLE :: G_interface_yB(:,:,:)
        !> Hydrostatic path contributions on top-oriented y-face traces.
        REAL(wp), ALLOCATABLE :: G_interface_yT(:,:,:)
        !> Cell-centered x-direction nonconservative pressure-path residuals.
        REAL(wp), ALLOCATABLE :: P_cell_x(:,:,:)
        !> Cell-centered y-direction nonconservative pressure-path residuals.
        REAL(wp), ALLOCATABLE :: P_cell_y(:,:,:)
        !> For L2 only, suppress negative reduced gravity in hydrostatic pressure.
        LOGICAL :: clamp_buoyant_pressure=.FALSE.
        !> Enable optional spatial reconstruction/flux timing.
        LOGICAL :: profile_enabled=.FALSE.
        !> Number of sampled hyperbolic evaluation calls.
        INTEGER :: profile_calls=0
        !> Cumulative reconstruction, flux and assembly timings [s].
        REAL(wp) :: profile_reconstruction=0.0_wp,profile_flux=0.0_wp,profile_assembly=0.0_wp
    CONTAINS
        PROCEDURE :: initialize => initialize_hyperbolic
        PROCEDURE :: finalize => finalize_hyperbolic
        PROCEDURE :: evaluate_terms => eval_hyperbolic_terms
        PROCEDURE :: evaluate_speeds => eval_speeds
        PROCEDURE :: evaluate_total_depth_speeds => eval_total_depth_speeds
        PROCEDURE :: diagnose_cell => diagnose_hyperbolic_cell
    END TYPE hyperbolic_workspace_type

CONTAINS

    !> @brief Estimate a shared two-layer wave envelope using the total depth.
    !> @details Each physically present layer uses the single-layer local-speed
    !> kernel with h=h1+h2. This recovers the single-layer speeds for equal
    !> material and velocity and encloses real gravitational speeds only under
    !> the stable-density assumptions of docs/TWO_LAYER_HYPERBOLICITY.md.
    !> Complex roots cannot be detected or repaired by this envelope.
    !> @param[in] q1_left Lower-layer primitive trace at the negative face endpoint.
    !> @param[in] q2_left Upper-layer primitive trace at the negative face endpoint.
    !> @param[in] q1_right Lower-layer primitive trace at the positive face endpoint.
    !> @param[in] q2_right Upper-layer primitive trace at the positive face endpoint.
    !> @param[in] metric Large-slope gravity factor at the shared face.
    !> @param[in] normal Normal velocity index (idx_u or idx_v).
    !> @param[out] a_minus Lower shared characteristic-speed envelope.
    !> @param[out] a_plus Upper shared characteristic-speed envelope.
    SUBROUTINE eval_total_depth_face_speed_bounds(q1_left,q2_left,q1_right,q2_right, &
        metric,normal,a_minus,a_plus)
        USE equation_terms_2d, ONLY : eval_local_speeds_direction
        REAL(wp), INTENT(IN) :: q1_left(n_vars),q2_left(n_vars),q1_right(n_vars),q2_right(n_vars)
        REAL(wp), INTENT(IN) :: metric
        INTEGER, INTENT(IN) :: normal
        REAL(wp), INTENT(OUT) :: a_minus,a_plus
        REAL(wp) :: q(n_vars),depth,lo,hi
        INTEGER :: endpoint,layer

        a_minus=0.0_wp
        a_plus=0.0_wp

        DO endpoint=1,2
            IF (endpoint==1) THEN
                depth=q1_left(1)+q2_left(1)
            ELSE
                depth=q1_right(1)+q2_right(1)
            END IF

            IF (depth<=dry_thickness_tolerance) CYCLE

            DO layer=1,2
                IF (endpoint==1) THEN
                    IF (layer==1) THEN
                        q=q1_left
                    ELSE
                        q=q2_left
                    END IF
                ELSE
                    IF (layer==1) THEN
                        q=q1_right
                    ELSE
                        q=q2_right
                    END IF
                END IF

                ! A dry layer cannot contribute an independent wave.
                IF (q(1)<=dry_thickness_tolerance) CYCLE
                q(1)=depth
                CALL eval_local_speeds_direction(q,metric,normal,lo,hi)
                a_minus=MIN(a_minus,lo)
                a_plus=MAX(a_plus,hi)
            END DO
        END DO
    END SUBROUTINE eval_total_depth_face_speed_bounds

    !> @brief Populate canonical face-speed storage using the total-depth envelope.
    !> @details Both reconstruction workspaces must be current on the shared OR
    !> face workset. The single-layer path remains unchanged. See
    !> docs/TWO_LAYER_HYPERBOLICITY.md for the conditional speed bound.
    !> @param[in,out] this Existing HP-PCCU face-speed workspace.
    !> @param[in] lower Prepared lower-layer reconstruction workspace.
    !> @param[in] upper Prepared upper-layer reconstruction workspace.
    !> @param[in] nxfaces Number of active x-direction faces.
    !> @param[in] jx X-direction face column indices.
    !> @param[in] kx X-direction face row indices.
    !> @param[in] nyfaces Number of active y-direction faces.
    !> @param[in] jy Y-direction face column indices.
    !> @param[in] ky Y-direction face row indices.
    SUBROUTINE eval_total_depth_speeds(this,lower,upper,nxfaces,jx,kx,nyfaces,jy,ky)
        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: lower,upper
        INTEGER, INTENT(IN) :: nxfaces,nyfaces,jx(:),kx(:),jy(:),ky(:)
        REAL(wp) :: q1a(n_vars),q1b(n_vars),q2a(n_vars),q2b(n_vars)
        INTEGER :: l,j,k

        IF (lower%layer/=1 .OR. upper%layer/=2) &
            ERROR STOP 'Total-depth bounds require both prepared layer reconstructions'

        IF (comp_cells_x>1) THEN
            !$OMP PARALLEL DO PRIVATE(l,j,k,q1a,q1b,q2a,q2b)
            DO l=1,nxfaces
                j=jx(l)
                k=kx(l)
                q1a=lower%qp_adv_interfaceL(:,j,k)
                q1b=lower%qp_adv_interfaceR(:,j,k)
                q2a=upper%qp_adv_interfaceL(:,j,k)
                q2b=upper%qp_adv_interfaceR(:,j,k)
                q1a(1)=lower%h_selected_interfaceL(j,k)
                q1b(1)=lower%h_selected_interfaceR(j,k)
                q2a(1)=upper%h_selected_interfaceL(j,k)
                q2b(1)=upper%h_selected_interfaceR(j,k)
                CALL eval_total_depth_face_speed_bounds(q1a,q2a,q1b,q2b, &
                    grav_coeff_stag_x(j,k),idx_u, &
                    this%a_interface_xNeg(j,k),this%a_interface_xPos(j,k))
            END DO
            !$OMP END PARALLEL DO
        END IF

        IF (comp_cells_y>1) THEN
            !$OMP PARALLEL DO PRIVATE(l,j,k,q1a,q1b,q2a,q2b)
            DO l=1,nyfaces
                j=jy(l)
                k=ky(l)
                q1a=lower%qp_adv_interfaceB(:,j,k)
                q1b=lower%qp_adv_interfaceT(:,j,k)
                q2a=upper%qp_adv_interfaceB(:,j,k)
                q2b=upper%qp_adv_interfaceT(:,j,k)
                q1a(1)=lower%h_selected_interfaceB(j,k)
                q1b(1)=lower%h_selected_interfaceT(j,k)
                q2a(1)=upper%h_selected_interfaceB(j,k)
                q2b(1)=upper%h_selected_interfaceT(j,k)
                CALL eval_total_depth_face_speed_bounds(q1a,q2a,q1b,q2b, &
                    grav_coeff_stag_y(j,k),idx_v, &
                    this%a_interface_yNeg(j,k),this%a_interface_yPos(j,k))
            END DO
            !$OMP END PARALLEL DO
        END IF
    END SUBROUTINE eval_total_depth_speeds

    !> @brief Accumulate coupled interlayer face and within-cell pressure in both residuals.
    !> @details Reuses the common pressure kernel and shared total-depth face-speed
    !> envelope on the OR workset. Subtract within-cell paths with
    !> the exact Gamma2 selected by the existing L2 hydrostatic cell kernel.
    !> Reconstructed layer traces on the shared workset.
    !> @param[in] hp2 Upper-layer pressure coefficients from the same traces.
    !> Existing layer spatial residuals.
    !> @param[in] active Shared active-cell mask.
    !> Shared active cell count and coordinates.
    !> @param[in] qp2 Upper-layer primitive state for dry-cell treatment.
    !> Shared face counts by Cartesian direction.
    !> Shared compact face coordinates.
    !> Shared total-depth face-speed bounds.
    !> @param[in] l1 Lower-layer reconstructed face traces.
    !> @param[in] l2 Upper-layer reconstructed face traces.
    !> @param[in,out] r1 Lower-layer residual accumulated in place.
    !> @param[in,out] r2 Upper-layer residual accumulated in place.
    !> @param[in] nc Count of shared active cells.
    !> @param[in] jc X indices of shared active cells.
    !> @param[in] kc Y indices of shared active cells.
    !> @param[in] nx Count of active x-normal faces.
    !> @param[in] jx X-face indices of active x-normal faces.
    !> @param[in] kx Y-cell indices of active x-normal faces.
    !> @param[in] ny Count of active y-normal faces.
    !> @param[in] jy X-cell indices of active y-normal faces.
    !> @param[in] ky Y-face indices of active y-normal faces.
    !> @param[in] axn Negative-going x-normal wave speeds [m/s].
    !> @param[in] axp Positive-going x-normal wave speeds [m/s].
    !> @param[in] ayn Negative-going y-normal wave speeds [m/s].
    !> @param[in] ayp Positive-going y-normal wave speeds [m/s].
    SUBROUTINE add_interlayer_terms(l1,l2,hp2,r1,r2,active,nc,jc,kc,qp2, &
        nx,jx,kx,ny,jy,ky,axn,axp,ayn,ayp)
        USE nonconservative_2d, ONLY : PATH_DIR_X,PATH_DIR_Y
        USE pccu_2d, ONLY : eval_interlayer_oriented_face,eval_interlayer_pressure_correction
        CLASS(reconstruction_workspace_type), INTENT(IN) :: l1,l2
        CLASS(hyperbolic_workspace_type), INTENT(IN) :: hp2
        REAL(wp), INTENT(INOUT) :: r1(n_eqns,comp_cells_x,comp_cells_y),r2(n_eqns,comp_cells_x,comp_cells_y)
        LOGICAL, INTENT(IN) :: active(comp_cells_x,comp_cells_y)
        INTEGER, INTENT(IN) :: nc,jc(:),kc(:)
        REAL(wp), INTENT(IN) :: qp2(n_vars,comp_cells_x,comp_cells_y)
        INTEGER, INTENT(IN) :: nx,ny,jx(:),kx(:),jy(:),ky(:)
        REAL(wp), INTENT(IN) :: axn(:,:),axp(:,:),ayn(:,:),ayp(:,:)
        REAL(wp) :: fl1(n_eqns),fr1(n_eqns),fl2(n_eqns),fr2(n_eqns)
        REAL(wp) :: path1(n_eqns),path2(n_eqns),hydro_path(n_eqns)
        REAL(wp) :: qp_a(n_vars),qp_b(n_vars)
        REAL(wp) :: h1a,h1b,h2a,h2b,gamma_a,gamma_b,eta_a,eta_b,metric,am,ap,spacing
        REAL(wp) :: metric_a,metric_b,selected_a,selected_b,cached_a,cached_b
        INTEGER :: direction,l,j,k,prev_j,prev_k,nfaces

        IF (l1%layer/=1 .OR. l2%layer/=2) ERROR STOP 'Invalid interlayer face ordering'
        IF (.NOT.hp2%clamp_buoyant_pressure) ERROR STOP 'Interlayer face requires L2 pressure'

        DO direction=PATH_DIR_X,PATH_DIR_Y
            nfaces=nx
            spacing=one_by_dx
            IF (direction==PATH_DIR_Y) THEN
                nfaces=ny
                spacing=one_by_dy
            END IF

            DO l=1,nfaces
                IF (direction==PATH_DIR_X) THEN
                    j=jx(l)
                    k=kx(l)
                    h1a=l1%h_selected_interfaceL(j,k)
                    h1b=l1%h_selected_interfaceR(j,k)
                    h2a=l2%h_selected_interfaceL(j,k)
                    h2b=l2%h_selected_interfaceR(j,k)
                    gamma_a=hp2%gamma_interface_xL(j,k)
                    gamma_b=hp2%gamma_interface_xR(j,k)
                    eta_a=l2%eta_interfaceL(j,k)
                    eta_b=l2%eta_interfaceR(j,k)
                    metric=grav_coeff_stag_x(j,k)
                    am=axn(j,k)
                    ap=axp(j,k)
                    prev_j=j-1
                    prev_k=k
                ELSE
                    j=jy(l)
                    k=ky(l)
                    h1a=l1%h_selected_interfaceB(j,k)
                    h1b=l1%h_selected_interfaceT(j,k)
                    h2a=l2%h_selected_interfaceB(j,k)
                    h2b=l2%h_selected_interfaceT(j,k)
                    gamma_a=hp2%gamma_interface_yB(j,k)
                    gamma_b=hp2%gamma_interface_yT(j,k)
                    eta_a=l2%eta_interfaceB(j,k)
                    eta_b=l2%eta_interfaceT(j,k)
                    metric=grav_coeff_stag_y(j,k)
                    am=ayn(j,k)
                    ap=ayp(j,k)
                    prev_j=j
                    prev_k=k-1
                END IF

                CALL eval_interlayer_oriented_face(direction,h1a,h2a,gamma_a,eta_a, &
                    h1b,h2b,gamma_b,eta_b,metric,am,ap,fl1,fr1,fl2,fr2)
                IF (prev_j>=1 .AND. prev_k>=1) THEN
                    IF (active(prev_j,prev_k)) THEN
                        r1(:,prev_j,prev_k)=r1(:,prev_j,prev_k)+fl1*spacing
                        r2(:,prev_j,prev_k)=r2(:,prev_j,prev_k)+fl2*spacing
                    END IF
                END IF
                IF (j<=comp_cells_x .AND. k<=comp_cells_y) THEN
                    IF (active(j,k)) THEN
                        r1(:,j,k)=r1(:,j,k)-fr1*spacing
                        r2(:,j,k)=r2(:,j,k)-fr2*spacing
                    END IF
                END IF
            END DO
        END DO

        ! Keep the same west/east and south/north endpoint ordering as P_cell.
        ! The Gamma2 outputs include the L2 dry-endpoint regularization.
        DO l=1,nc
            j=jc(l)
            k=kc(l)

            IF (.NOT.active(j,k)) ERROR STOP 'Interlayer cell outside shared workset'
            IF (qp2(1,j,k)<=dry_thickness_tolerance) CYCLE

            DO direction=PATH_DIR_X,PATH_DIR_Y
                IF (direction==PATH_DIR_X .AND. comp_cells_x<=1) CYCLE
                IF (direction==PATH_DIR_Y .AND. comp_cells_y<=1) CYCLE

                IF (direction==PATH_DIR_X) THEN
                    h1a=MAX(l1%eta_cellW(j,k)-B_face_x(j,k),0.0_wp)
                    h1b=MAX(l1%eta_cellE(j,k)-B_face_x(j+1,k),0.0_wp)
                    h2a=MAX(l2%eta_cellW(j,k)-l1%eta_cellW(j,k),0.0_wp)
                    h2b=MAX(l2%eta_cellE(j,k)-l1%eta_cellE(j,k),0.0_wp)
                    eta_a=l2%eta_cellW(j,k)
                    eta_b=l2%eta_cellE(j,k)
                    metric_a=grav_coeff_stag_x(j,k)
                    metric_b=grav_coeff_stag_x(j+1,k)
                    selected_a=l2%h_selected_interfaceR(j,k)
                    selected_b=l2%h_selected_interfaceL(j+1,k)
                    gamma_a=hp2%gamma_interface_xR(j,k)
                    gamma_b=hp2%gamma_interface_xL(j+1,k)
                    qp_a=l2%qp_adv_interfaceR(:,j,k)
                    qp_b=l2%qp_adv_interfaceL(:,j+1,k)
                    spacing=one_by_dx
                ELSE
                    h1a=MAX(l1%eta_cellS(j,k)-B_face_y(j,k),0.0_wp)
                    h1b=MAX(l1%eta_cellN(j,k)-B_face_y(j,k+1),0.0_wp)
                    h2a=MAX(l2%eta_cellS(j,k)-l1%eta_cellS(j,k),0.0_wp)
                    h2b=MAX(l2%eta_cellN(j,k)-l1%eta_cellN(j,k),0.0_wp)
                    eta_a=l2%eta_cellS(j,k)
                    eta_b=l2%eta_cellN(j,k)
                    metric_a=grav_coeff_stag_y(j,k)
                    metric_b=grav_coeff_stag_y(j,k+1)
                    selected_a=l2%h_selected_interfaceT(j,k)
                    selected_b=l2%h_selected_interfaceB(j,k+1)
                    gamma_a=hp2%gamma_interface_yT(j,k)
                    gamma_b=hp2%gamma_interface_yB(j,k+1)
                    qp_a=l2%qp_adv_interfaceT(:,j,k)
                    qp_b=l2%qp_adv_interfaceB(:,j,k+1)
                    spacing=one_by_dy
                END IF

                qp_a(1)=h2a
                qp_b(1)=h2b
                IF (h2a<=dry_thickness_tolerance) qp_a(2:3)=0.0_wp
                IF (h2b<=dry_thickness_tolerance) qp_b(2:3)=0.0_wp

                ! Return exactly the pressure coefficients used by the L2 path.
                cached_a=gamma_a
                cached_b=gamma_b
                CALL eval_cell_hydrostatic_path(qp_a,eta_a,metric_a,qp_b,eta_b,metric_b, &
                    direction,hydro_path,cached_a,selected_a,cached_b,selected_b, &
                    hp2%clamp_buoyant_pressure,gamma_minus_out=gamma_a,gamma_plus_out=gamma_b)
                CALL eval_interlayer_pressure_correction(direction, &
                    h1a,h2a,gamma_a,eta_a,metric_a, &
                    h1b,h2b,gamma_b,eta_b,metric_b,path1,path2)
                r1(:,j,k)=r1(:,j,k)-path1*spacing
                r2(:,j,k)=r2(:,j,k)-path2*spacing
            END DO
        END DO
    END SUBROUTINE add_interlayer_terms

    !> @brief Report the resolved face contributions to one conservative residual.
    !> @param[in] this Face-oriented HP-PCCU workspace.
    !> @param[in] recon Reconstructed primitive traces on each face.
    !> @param[in] j Horizontal index of the failing cell.
    !> @param[in] k Vertical index of the failing cell.
    !> @param[in] eq Conservative equation index.
    SUBROUTINE diagnose_hyperbolic_cell(this,recon,j,k,eq)
        CLASS(hyperbolic_workspace_type), INTENT(IN) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        INTEGER, INTENT(IN) :: j,k,eq
        REAL(wp) :: x_right,x_left,y_top,y_bottom,path_x,path_y
        x_right=0.0_wp
        x_left=0.0_wp
        y_top=0.0_wp
        y_bottom=0.0_wp
        path_x=0.0_wp
        path_y=0.0_wp
        IF (comp_cells_x>1) THEN
            x_right=this%G_interface_xL(eq,j+1,k)
            x_left=this%G_interface_xR(eq,j,k)
            path_x=this%P_cell_x(eq,j,k)
        END IF
        IF (comp_cells_y>1) THEN
            y_top=this%G_interface_yB(eq,j,k+1)
            y_bottom=this%G_interface_yT(eq,j,k)
            path_y=this%P_cell_y(eq,j,k)
        END IF
        WRITE(*,'(A,I0,A,6(1X,ES23.15))') 'LIE_FACE_FAILURE equation=',eq, &
            ' right/left/top/bottom/path_x/path_y=', &
            x_right,x_left,y_top,y_bottom,path_x,path_y
        WRITE(*,'(A,I0,A,3(1X,ES23.15))') 'LIE_FACE_FAILURE equation=',eq, &
            ' div_x/div_y/total=', &
            (x_right-x_left-path_x)*one_by_dx, &
            (y_top-y_bottom-path_y)*one_by_dy, &
            (x_right-x_left-path_x)*one_by_dx+(y_top-y_bottom-path_y)*one_by_dy
        IF (eq>=idx_solidEqn_first .AND. eq<=idx_solidEqn_last) THEN
            WRITE(*,'(A,4(1X,ES23.15))') 'LIE_TRACE_FAILURE y_upper solidB/solidT/hB/hT=', &
                recon%qp_adv_interfaceB(idx_solid_first+eq-idx_solidEqn_first,j,k+1), &
                recon%qp_adv_interfaceT(idx_solid_first+eq-idx_solidEqn_first,j,k+1), &
                recon%h_selected_interfaceB(j,k+1),recon%h_selected_interfaceT(j,k+1)
            WRITE(*,'(A,4(1X,ES23.15))') 'LIE_TRACE_FAILURE y_upper vB/vT/amin/aplus=', &
                recon%qp_adv_interfaceB(idx_v,j,k+1),recon%qp_adv_interfaceT(idx_v,j,k+1), &
                this%a_interface_yNeg(j,k+1),this%a_interface_yPos(j,k+1)
        END IF
    END SUBROUTINE diagnose_hyperbolic_cell

    !> \brief Allocate and zero wave-speed and oriented KT work arrays.
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
    !> @param[in] clamp_buoyant_pressure Suppress negative hydrostatic pressure for the dilute upper layer.

    SUBROUTINE initialize_hyperbolic( this, clamp_buoyant_pressure )

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        LOGICAL, INTENT(IN), OPTIONAL :: clamp_buoyant_pressure

        BLOCK
            CHARACTER(LEN=8) :: value
            INTEGER :: status
            value=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_PROFILE_LIE',value,STATUS=status)
            this%profile_enabled=status==0 .AND. TRIM(value)=='1'
        END BLOCK
        this%clamp_buoyant_pressure=.FALSE.
        IF (PRESENT(clamp_buoyant_pressure)) this%clamp_buoyant_pressure=clamp_buoyant_pressure
        ALLOCATE( this%a_interface_xNeg(comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%a_interface_xPos(comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%a_interface_yNeg(comp_cells_x,comp_interfaces_y) )
        ALLOCATE( this%a_interface_yPos(comp_cells_x,comp_interfaces_y) )
        ALLOCATE( this%gamma_interface_xL(comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%gamma_interface_xR(comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%gamma_interface_yB(comp_cells_x,comp_interfaces_y) )
        ALLOCATE( this%gamma_interface_yT(comp_cells_x,comp_interfaces_y) )

        ALLOCATE( this%G_interface_xL(n_eqns,comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%G_interface_xR(n_eqns,comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%G_interface_yB(n_eqns,comp_cells_x,comp_interfaces_y) )
        ALLOCATE( this%G_interface_yT(n_eqns,comp_cells_x,comp_interfaces_y) )
        ALLOCATE( this%P_cell_x(n_eqns,comp_cells_x,comp_cells_y) )
        ALLOCATE( this%P_cell_y(n_eqns,comp_cells_x,comp_cells_y) )

        this%a_interface_xNeg = 0.0_wp
        this%a_interface_xPos = 0.0_wp
        this%a_interface_yNeg = 0.0_wp
        this%a_interface_yPos = 0.0_wp
        this%gamma_interface_xL = 0.0_wp
        this%gamma_interface_xR = 0.0_wp
        this%gamma_interface_yB = 0.0_wp
        this%gamma_interface_yT = 0.0_wp
        this%G_interface_xL = 0.0_wp
        this%G_interface_xR = 0.0_wp
        this%G_interface_yB = 0.0_wp
        this%G_interface_yT = 0.0_wp
        this%P_cell_x = 0.0_wp
        this%P_cell_y = 0.0_wp

    END SUBROUTINE initialize_hyperbolic

    !> \brief Assemble the complete HP-PCCU spatial term on active cells.
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
    !> \param[in,out] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
    !> \param[in] qp_expl Cell-centered physical states at the explicit stage, including the
    !> required
    !>                    stencil halo.
    !> \param[out] divFlux_iRK Complete stage spatial term, including face and cell paths, indexed
    !>                         (equation,x,y).
    !> \param[in] t Current simulation or stage time [s].
    !> \param[in] solve_cells Number of entries in the active-cell index lists.
    !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
    !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
    !> \param[in] solve_interfaces_x Number of active x-normal faces.
    !> \param[in] j_stag_x X-face indices paired with k_stag_x.
    !> \param[in] k_stag_x Y-cell indices of active x-normal faces.
    !> \param[in] solve_interfaces_y Number of active y-normal faces.
    !> \param[in] j_stag_y X-cell indices of active y-normal faces.
    !> \param[in] k_stag_y Y-face indices paired with j_stag_y.
    !> \param[in] source_boundary_enabled Optional override for boundary-source traces.

    !> @param[in] reuse_cfl_traces Flag allowing reuse of CFL reconstruction traces.
    !> @param[in] lower_recon Lower-layer reconstruction.
    !> @param[in] h_lower Lower-layer depth.
    SUBROUTINE eval_hyperbolic_terms( this, recon, qp_expl,                                                            &
            divFlux_iRK, t, solve_cells, j_cent, k_cent, solve_interfaces_x,                                           &
            j_stag_x, k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y,                                                &
            source_boundary_enabled,reuse_cfl_traces,lower_recon,h_lower )

    IMPLICIT NONE

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: recon
        REAL(wp), INTENT(IN) :: qp_expl(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(OUT) :: divFlux_iRK(n_eqns,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: t
        INTEGER, INTENT(IN) :: solve_cells
        INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
        INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
        INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
        INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)
        LOGICAL, INTENT(IN), OPTIONAL :: source_boundary_enabled,reuse_cfl_traces
        CLASS(reconstruction_workspace_type), INTENT(IN), OPTIONAL :: lower_recon
        REAL(wp), INTENT(IN), OPTIONAL :: h_lower(:,:)
        LOGICAL :: use_cached_traces

        INTEGER :: l , i, j, k      !< loop counters
        INTEGER :: tick,tock,tick_rate
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tick,tick_rate)
            this%profile_calls=this%profile_calls+1
        END IF

        ! The CFL pass has already reconstructed the same unchanged state and
        ! evaluated these face speeds. Only production explicitly opts in.
        use_cached_traces=.FALSE.
        IF (PRESENT(reuse_cfl_traces)) use_cached_traces=reuse_cfl_traces
        IF (PRESENT(lower_recon) .NEQV. PRESENT(h_lower)) &
            ERROR STOP 'L2 residual requires lower reconstruction and thickness'
        IF (PRESENT(lower_recon)) THEN
            IF (recon%layer/=2 .OR. lower_recon%layer/=1) ERROR STOP 'Invalid L1/L2 HP pairing'
        END IF
        IF (.NOT.use_cached_traces) THEN
            IF (PRESENT(lower_recon)) THEN
                CALL recon%reconstruct(qp_expl,t,solve_cells,j_cent,k_cent,.TRUE., &
                    source_boundary_enabled=source_boundary_enabled,lower_recon=lower_recon,h_lower=h_lower)
                CALL recon%reconstruct_kt_eta(qp_expl,lower_recon=lower_recon)
            ELSE
                CALL recon%reconstruct(qp_expl,t,solve_cells,j_cent,k_cent,.TRUE.,                                       &
                    source_boundary_enabled=source_boundary_enabled )
                CALL recon%reconstruct_kt_eta(qp_expl)
            END IF
            CALL this%evaluate_speeds( recon, solve_interfaces_x, j_stag_x,                                            &
                k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )
        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_reconstruction=this%profile_reconstruction+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        ! The production spatial operator is the single HP-PCCU baseline.
        CALL eval_flux_PCCU( this, recon, qp_expl,                                                                     &
            .TRUE., t,                                                                                                 &
            solve_interfaces_x, j_stag_x,                                                                              &
            k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y, lower_recon )

        CALL eval_cell_hydrostatic_paths( this, recon, qp_expl,                                                        &
            .TRUE., solve_cells, j_cent, k_cent, lower_recon )

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_flux=this%profile_flux+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            tick=tock
        END IF
        ! Assemble each cell from its outward-oriented face values, then subtract
        ! the within-cell hydrostatic path. This cancellation is essential for
        ! lake-at-rest balance; adding a separate bed-force source would count it twice.
        !$OMP PARALLEL DO private(l,j,k,i)

        cells_loop:DO l = 1,solve_cells

            j = j_cent(l)
            k = k_cent(l)

            DO i=1,n_eqns

                divFlux_iRK(i,j,k) = 0.0_wp

                IF ( comp_cells_x .GT. 1 ) THEN

                    divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                                                          &
                        ( this%G_interface_xL(i,j+1,k)                                                                 &
                        - this%G_interface_xR(i,j,k)                                                                   &
                        - this%P_cell_x(i,j,k) ) * one_by_dx

                END IF

                IF ( comp_cells_y .GT. 1 ) THEN

                    divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                                                          &
                        ( this%G_interface_yB(i,j,k+1)                                                                 &
                        - this%G_interface_yT(i,j,k)                                                                   &
                        - this%P_cell_y(i,j,k) ) * one_by_dy

                END IF

            END DO

            ! Only a completely stationary stencil may discard hydrostatic
            ! cancellation noise. Never classify a nonzero physical velocity as
            ! rest: even tiny motion must retain the unmodified spatial operator.
            ! A pressure-scaled roundoff residual otherwise creates a tiny discharge
            ! which disables the exact-rest scalar-flux guard at the next stage.
            IF ( exactly_resting_stencil(recon,qp_expl(:,j,k),j,k) ) THEN
                IF ( comp_cells_x .GT. 1 ) CALL cancel_resting_roundoff(                                               &
                    this,recon,j,k,1,divFlux_iRK(2,j,k) )
                IF ( comp_cells_y .GT. 1 ) CALL cancel_resting_roundoff(                                               &
                    this,recon,j,k,2,divFlux_iRK(3,j,k) )

            END IF

        END DO cells_loop

        !$OMP END PARALLEL DO
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_assembly=this%profile_assembly+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
        RETURN
    END SUBROUTINE eval_hyperbolic_terms

    !> \brief Require exact stationarity of a cell and every adjacent face endpoint.
    !> \param[in] recon Final reconstructed physical interface states.
    !> \param[in] qp_cell Physical cell state; entries 2:3 are volumetric momenta.
    !> \param[in] j Cell x index.
    !> \param[in] k Cell y index.
    !> \return True only when all participating momenta are exactly zero.
    PURE FUNCTION exactly_resting_stencil(recon,qp_cell,j,k) RESULT(resting)
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        REAL(wp), INTENT(IN) :: qp_cell(:)
        INTEGER, INTENT(IN) :: j,k
        LOGICAL :: resting

        resting = ALL(qp_cell(2:3) .EQ. 0.0_wp)
        IF ( .NOT. resting ) RETURN

        IF ( comp_cells_x .GT. 1 ) THEN
            resting = ALL(recon%qp_adv_interfaceL(2:3,j:j+1,k) .EQ. 0.0_wp) .AND.                                      &
                ALL(recon%qp_adv_interfaceR(2:3,j:j+1,k) .EQ. 0.0_wp)
            IF ( .NOT. resting ) RETURN

        END IF

        IF ( comp_cells_y .GT. 1 ) THEN
            resting = ALL(recon%qp_adv_interfaceB(2:3,j,k:k+1) .EQ. 0.0_wp) .AND.                                      &
                ALL(recon%qp_adv_interfaceT(2:3,j,k:k+1) .EQ. 0.0_wp)

        END IF
    END FUNCTION exactly_resting_stencil

    !> \brief Remove only pressure-scaled arithmetic noise from stationary momentum balance.
    !> \details Called only for an exactly stationary stencil. The fixed 64-epsilon
    !>          allowance covers closure, path quadrature and the three-term
    !>          cancellation; the scale retains sensitivity to absolute eta in
    !>          path differences. It has no absolute force or velocity floor.
    !>          Scalar fluxes and moving states are never modified here.
    !> \param[in] recon Final face states and free-surface traces.
    !> \param[in] j Cell x index.
    !> \param[in] k Cell y index.
    !> \param[in] axis One for x momentum, two for y momentum.
    !> \param[in,out] residual Assembled normal-momentum spatial term.
    !> @param[in] this Owning object.
    SUBROUTINE cancel_resting_roundoff(this,recon,j,k,axis,residual)
        CLASS(hyperbolic_workspace_type), INTENT(IN) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        INTEGER, INTENT(IN) :: j,k,axis
        REAL(wp), INTENT(INOUT) :: residual
        REAL(wp) :: scale,inverse_spacing

        ! An already vanishing momentum residual requires no thermodynamics.
        IF (residual==0.0_wp) THEN
            residual=0.0_wp
            RETURN
        END IF

        IF (axis .EQ. 1) THEN
            scale = MAX(                                                                                               &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceL(:,j,k),                                           &
                recon%eta_interfaceL(j,k),grav_coeff_stag_x(j,k),                                                     &
                this%gamma_interface_xL(j,k),recon%h_selected_interfaceL(j,k),                                         &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceR(:,j,k),                                           &
                recon%eta_interfaceR(j,k),grav_coeff_stag_x(j,k),                                                     &
                this%gamma_interface_xR(j,k),recon%h_selected_interfaceR(j,k),                                         &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceL(:,j+1,k),                                         &
                recon%eta_interfaceL(j+1,k),grav_coeff_stag_x(j+1,k),                                                 &
                this%gamma_interface_xL(j+1,k),recon%h_selected_interfaceL(j+1,k),                                     &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceR(:,j+1,k),                                         &
                recon%eta_interfaceR(j+1,k),grav_coeff_stag_x(j+1,k),                                                 &
                this%gamma_interface_xR(j+1,k),recon%h_selected_interfaceR(j+1,k),                                     &
                this%clamp_buoyant_pressure) )
            inverse_spacing = one_by_dx
        ELSE
            scale = MAX(                                                                                               &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceB(:,j,k),                                           &
                recon%eta_interfaceB(j,k),grav_coeff_stag_y(j,k),                                                     &
                this%gamma_interface_yB(j,k),recon%h_selected_interfaceB(j,k),                                         &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceT(:,j,k),                                           &
                recon%eta_interfaceT(j,k),grav_coeff_stag_y(j,k),                                                     &
                this%gamma_interface_yT(j,k),recon%h_selected_interfaceT(j,k),                                         &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceB(:,j,k+1),                                         &
                recon%eta_interfaceB(j,k+1),grav_coeff_stag_y(j,k+1),                                                 &
                this%gamma_interface_yB(j,k+1),recon%h_selected_interfaceB(j,k+1),                                     &
                this%clamp_buoyant_pressure),                                                                          &
                hydrostatic_arithmetic_scale(recon%qp_adv_interfaceT(:,j,k+1),                                         &
                recon%eta_interfaceT(j,k+1),grav_coeff_stag_y(j,k+1),                                                 &
                this%gamma_interface_yT(j,k+1),recon%h_selected_interfaceT(j,k+1),                                     &
                this%clamp_buoyant_pressure) )
            inverse_spacing = one_by_dy

        END IF

        IF ( ABS(residual) .LE. 64.0_wp*EPSILON(1.0_wp)*scale*inverse_spacing )                                        &
            residual = 0.0_wp
    END SUBROUTINE cancel_resting_roundoff

    !> Estimate the hydrostatic arithmetic scale for one unmodified face trace.
    !> The cached coefficient comes from the selected-depth CFL trace and uses
    !> identical temperature and mass fractions. Only the depth can differ.
    !> If it changes the wet/dry branch, retain the original evaluation.
    !> @param[in] qp_face Final primitive face state.
    !> @param[in] eta Reconstructed free-surface elevation.
    !> @param[in] gravity_factor Terrain large-slope metric.
    !> @param[in] cached_gamma Hydrostatic coefficient from the CFL face trace.
    !> @param[in] selected_h Depth at which cached_gamma was evaluated.
    !> @param[in] clamp_buoyant_pressure Apply the upper-layer nonnegative-pressure condition.
    !> @return Absolute pressure scale used only for stationary-roundoff cancellation.
    !> @brief Estimate hydrostatic pressure scale for roundoff cancellation.
    FUNCTION hydrostatic_arithmetic_scale(qp_face,eta,gravity_factor,cached_gamma,selected_h,                          &
        clamp_buoyant_pressure) RESULT(scale)
    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
        REAL(wp), INTENT(IN) :: qp_face(:),eta,gravity_factor,cached_gamma,selected_h
        LOGICAL, INTENT(IN) :: clamp_buoyant_pressure
        REAL(wp) :: scale,gamma,reduced_gravity,h

        h = MAX(0.0_wp,qp_face(1))
        scale = 0.0_wp
        IF (h .EQ. 0.0_wp) RETURN
        IF (h<=EPSILON(1.0_wp)) THEN
            gamma=0.0_wp
        ELSEIF (selected_h>EPSILON(1.0_wp)) THEN
            gamma=cached_gamma
        ELSE
            CALL eval_hydrostatic_coefficient(qp_face,reduced_gravity,gamma)
        END IF
        IF (clamp_buoyant_pressure) gamma=MAX(gamma,0.0_wp)
        scale = ABS(gravity_factor*gamma)*h*MAX(h,ABS(eta))
    END FUNCTION hydrostatic_arithmetic_scale

    !******************************************************************************
    !> \brief Build both oriented KT values at every active face.
    !>
    !> The conservative endpoint flux is purely inertial.  Hydrostatic pressure
    !> and bed geometry enter only through the path contribution used to form the
    !> two oriented values at each face.
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
    !> \param[in] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
    !> \param[in] qp_center Cell-centered physical state used for face fluxes.
    !> \param[in] use_dry_donor Enable the dry-donor face treatment.
    !> \param[in] t Current simulation or stage time [s].
    !> \param[in] solve_interfaces_x Number of active x-normal faces.
    !> \param[in] j_stag_x X-face indices paired with k_stag_x.
    !> \param[in] k_stag_x Y-cell indices of active x-normal faces.
    !> \param[in] solve_interfaces_y Number of active y-normal faces.
    !> \param[in] j_stag_y X-cell indices of active y-normal faces.
    !> \param[in] k_stag_y Y-face indices paired with j_stag_y.
    !******************************************************************************

    !> @brief Evaluate oriented PCCU numerical fluxes for active faces.
    !> @param[in] lower_recon Lower-layer reconstruction.
    SUBROUTINE eval_flux_PCCU( this, recon, qp_center, use_dry_donor, t, solve_interfaces_x, j_stag_x,                 &
            k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y, lower_recon )

    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

    IMPLICIT NONE

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        REAL(wp), INTENT(IN) :: qp_center(n_vars,comp_cells_x,comp_cells_y)
        LOGICAL, INTENT(IN) :: use_dry_donor
        REAL(wp), INTENT(IN) :: t
        INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
        INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
        INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)
        CLASS(reconstruction_workspace_type), INTENT(IN), OPTIONAL :: lower_recon

        REAL(wp) :: a_minus, a_plus
        REAL(wp) :: source_inflow, other_inflow, face_inflow
        LOGICAL :: dry_donor_left,dry_donor_right,front_trace_enabled,trace_front_face
        INTEGER :: j, k, l, diagnostic_status
        CHARACTER(LEN=1) :: diagnostic_flag

        diagnostic_flag=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_FRONT_MASS_DIAGNOSTICS',diagnostic_flag,STATUS=diagnostic_status)
        front_trace_enabled=diagnostic_status==0 .AND. diagnostic_flag=='1'

        ! X- and Y-face kernels write disjoint oriented arrays. Keep one
        ! OpenMP team and let workers begin Y faces as soon as X work ends.
        !$OMP PARALLEL IF(comp_cells_x>1 .OR. comp_cells_y>1) &
        !$OMP& PRIVATE(l,j,k,a_minus,a_plus,dry_donor_left,dry_donor_right,trace_front_face)
        IF ( comp_cells_x .GT. 1 ) THEN

            !$OMP DO
            DO l = 1, solve_interfaces_x

                j = j_stag_x(l)
                k = k_stag_x(l)

                a_minus = this%a_interface_xNeg(j,k)
                a_plus = this%a_interface_xPos(j,k)
                dry_donor_left = recon%qp_adv_interfaceL(1,j,k) .LE. dry_thickness_tolerance
                dry_donor_right = recon%qp_adv_interfaceR(1,j,k) .LE. dry_thickness_tolerance
                IF ( j .GT. 1 ) dry_donor_left =                                                                       &
                    qp_center(1,j-1,k) .LE. dry_thickness_tolerance
                IF ( j .LE. comp_cells_x ) dry_donor_right =                                                           &
                    qp_center(1,j,k) .LE. dry_thickness_tolerance
                trace_front_face=.FALSE.
                IF (front_trace_enabled .AND. recon%layer==1) THEN
                    IF (j>1 .AND. j<=comp_cells_x) THEN
                        trace_front_face=qp_center(1,j-1,k)>dry_thickness_tolerance .AND. &
                            qp_center(1,j,k)<=dry_thickness_tolerance
                    END IF
                END IF
                IF (trace_front_face) THEN
                    WRITE(*,'(A,2(I0,1X),A,ES24.16)') &
                        'FRONT_KERNEL enter face(j,k)=',j,k,' t=',t
                END IF
                IF (PRESENT(lower_recon)) THEN
                    CALL eval_pccu_interface( recon%eta_interfaceL(j,k), recon%eta_interfaceR(j,k), B_face_x(j,k),         &
                        recon%h_selected_interfaceL(j,k),                                                              &
                        recon%h_selected_interfaceR(j,k),                                                              &
                    recon%qp_adv_interfaceL(:,j,k),                                                                    &
                    recon%qp_adv_interfaceR(:,j,k), grav_coeff_stag_x(j,k),                                            &
                    a_minus, a_plus, this%gamma_interface_xL(j,k),this%gamma_interface_xR(j,k),                      &
                    PATH_DIR_X, 2, j, k, use_dry_donor,                                               &
                    dry_donor_left,dry_donor_right,t,                                                                  &
                    this%G_interface_xL(:,j,k), this%G_interface_xR(:,j,k), &
                    bed_face_left=lower_recon%eta_interfaceL(j,k), &
                    bed_face_right=lower_recon%eta_interfaceR(j,k),diagnose_front=trace_front_face )
                ELSE
                    CALL eval_pccu_interface( recon%eta_interfaceL(j,k), recon%eta_interfaceR(j,k), B_face_x(j,k),         &
                        recon%h_selected_interfaceL(j,k),                                                              &
                        recon%h_selected_interfaceR(j,k),                                                              &
                    recon%qp_adv_interfaceL(:,j,k),                                                                    &
                    recon%qp_adv_interfaceR(:,j,k), grav_coeff_stag_x(j,k),                                            &
                    a_minus, a_plus, this%gamma_interface_xL(j,k),this%gamma_interface_xR(j,k),                      &
                    PATH_DIR_X, 2, j, k, use_dry_donor,                                               &
                    dry_donor_left,dry_donor_right,t,                                                                  &
                    this%G_interface_xL(:,j,k), this%G_interface_xR(:,j,k), &
                    diagnose_front=trace_front_face )
                END IF
                IF (trace_front_face) THEN
                    WRITE(*,'(A,2(ES24.16,1X))') 'FRONT_KERNEL exit left,right=', &
                        this%G_interface_xL(1,j,k),this%G_interface_xR(1,j,k)
                END IF

                ! Apply the boundary's one-way condition to the final numerical flux.
                ! V2 limited the reconstructed normal momentum to outward-only values
                ! (west/south <= 0; east/north >= 0). HP reconstruction follows the
                ! primitive boundary trace and can otherwise reintroduce inflow, so
                ! remove incoming mass-carried transport here. Keep normal momentum
                ! flux, which contains the wall pressure contribution. Source cells
                ! are exempt, as in the legacy reconstruction.
                IF ( ( j .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),                                      &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcW(2,recon%layer)%flag .NE. 0 ) .AND.                                                           &
                    ( this%G_interface_xL(1,j,k) .GT. 0.0_wp ) ) THEN
                    this%G_interface_xL(1,j,k) = 0.0_wp
                    this%G_interface_xR(1,j,k) = 0.0_wp
                    this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
                    this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( j .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),                                  &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcW(2,recon%layer)%flag .EQ. 0 ) .AND.                                                           &
                    ( bcW(2,recon%layer)%value .EQ. 0.0_wp ) ) THEN
                    this%G_interface_xL(1,j,k) = 0.0_wp
                    this%G_interface_xR(1,j,k) = 0.0_wp
                    this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
                    this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( j .EQ. comp_interfaces_x ) .AND.                                                            &
                    ( source_cell(MAX(1,MIN(j-1,SIZE(source_cell,1))),                                                 &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcE(2,recon%layer)%flag .NE. 0 ) .AND.                                                           &
                    ( this%G_interface_xL(1,j,k) .LT. 0.0_wp ) ) THEN
                    this%G_interface_xL(1,j,k) = 0.0_wp
                    this%G_interface_xR(1,j,k) = 0.0_wp
                    this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
                    this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( j .EQ. comp_interfaces_x ) .AND.                                                            &
                    ( source_cell(MAX(1,MIN(j-1,SIZE(source_cell,1))),                                                 &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcE(2,recon%layer)%flag .EQ. 0 ) .AND.                                                           &
                    ( bcE(2,recon%layer)%value .EQ. 0.0_wp ) ) THEN
                    this%G_interface_xL(1,j,k) = 0.0_wp
                    this%G_interface_xR(1,j,k) = 0.0_wp
                    this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
                    this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp

                END IF

            END DO
            !$OMP END DO NOWAIT

        END IF

        IF ( comp_cells_y .GT. 1 ) THEN

            !$OMP DO
            DO l = 1, solve_interfaces_y

                j = j_stag_y(l)
                k = k_stag_y(l)

                a_minus = this%a_interface_yNeg(j,k)
                a_plus = this%a_interface_yPos(j,k)
                dry_donor_left = recon%qp_adv_interfaceB(1,j,k) .LE. dry_thickness_tolerance
                dry_donor_right = recon%qp_adv_interfaceT(1,j,k) .LE. dry_thickness_tolerance
                IF ( k .GT. 1 ) dry_donor_left =                                                                       &
                    qp_center(1,j,k-1) .LE. dry_thickness_tolerance
                IF ( k .LE. comp_cells_y ) dry_donor_right =                                                           &
                    qp_center(1,j,k) .LE. dry_thickness_tolerance
                IF (PRESENT(lower_recon)) THEN
                    CALL eval_pccu_interface( recon%eta_interfaceB(j,k), recon%eta_interfaceT(j,k), B_face_y(j,k),         &
                        recon%h_selected_interfaceB(j,k),                                                              &
                        recon%h_selected_interfaceT(j,k),                                                              &
                    recon%qp_adv_interfaceB(:,j,k),                                                                    &
                    recon%qp_adv_interfaceT(:,j,k), grav_coeff_stag_y(j,k),                                            &
                    a_minus, a_plus, this%gamma_interface_yB(j,k),this%gamma_interface_yT(j,k),                      &
                    PATH_DIR_Y, 3, j, k, use_dry_donor,                                               &
                    dry_donor_left,dry_donor_right,t,                                                                  &
                    this%G_interface_yB(:,j,k), this%G_interface_yT(:,j,k), &
                    bed_face_left=lower_recon%eta_interfaceB(j,k), &
                    bed_face_right=lower_recon%eta_interfaceT(j,k) )
                ELSE
                    CALL eval_pccu_interface( recon%eta_interfaceB(j,k), recon%eta_interfaceT(j,k), B_face_y(j,k),         &
                        recon%h_selected_interfaceB(j,k),                                                              &
                        recon%h_selected_interfaceT(j,k),                                                              &
                    recon%qp_adv_interfaceB(:,j,k),                                                                    &
                    recon%qp_adv_interfaceT(:,j,k), grav_coeff_stag_y(j,k),                                            &
                    a_minus, a_plus, this%gamma_interface_yB(j,k),this%gamma_interface_yT(j,k),                      &
                    PATH_DIR_Y, 3, j, k, use_dry_donor,                                               &
                    dry_donor_left,dry_donor_right,t,                                                                  &
                    this%G_interface_yB(:,j,k), this%G_interface_yT(:,j,k) )
                END IF

                IF ( ( k .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),                                      &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcS(3,recon%layer)%flag .NE. 0 ) .AND.                                                           &
                    ( this%G_interface_yB(1,j,k) .GT. 0.0_wp ) ) THEN
                    this%G_interface_yB(1:2,j,k) = 0.0_wp
                    this%G_interface_yT(1:2,j,k) = 0.0_wp
                    this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
                    this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( k .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),                                  &
                    MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                                                         &
                    ( bcS(3,recon%layer)%flag .EQ. 0 ) .AND.                                                           &
                    ( bcS(3,recon%layer)%value .EQ. 0.0_wp ) ) THEN
                    this%G_interface_yB(1:2,j,k) = 0.0_wp
                    this%G_interface_yT(1:2,j,k) = 0.0_wp
                    this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
                    this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( k .EQ. comp_interfaces_y ) .AND.                                                            &
                    ( source_cell(MIN(j,SIZE(source_cell,1)),                                                          &
                    MAX(1,MIN(k-1,SIZE(source_cell,2)))) .NE. 2 ) .AND.                                                &
                    ( bcN(3,recon%layer)%flag .NE. 0 ) .AND.                                                           &
                    ( this%G_interface_yB(1,j,k) .LT. 0.0_wp ) ) THEN
                    this%G_interface_yB(1:2,j,k) = 0.0_wp
                    this%G_interface_yT(1:2,j,k) = 0.0_wp
                    this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
                    this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
                ELSEIF ( ( k .EQ. comp_interfaces_y ) .AND.                                                            &
                    ( source_cell(MIN(j,SIZE(source_cell,1)),                                                          &
                    MAX(1,MIN(k-1,SIZE(source_cell,2)))) .NE. 2 ) .AND.                                                &
                    ( bcN(3,recon%layer)%flag .EQ. 0 ) .AND.                                                           &
                    ( bcN(3,recon%layer)%value .EQ. 0.0_wp ) ) THEN
                    this%G_interface_yB(1:2,j,k) = 0.0_wp
                    this%G_interface_yT(1:2,j,k) = 0.0_wp
                    this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
                    this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp

                END IF

            END DO
            !$OMP END DO

        END IF
        !$OMP END PARALLEL

        ! Optional stage-level boundary budget for source-mass CI diagnostics.
        ! The oriented mass contributions are measured on the interior side of
        ! each domain boundary, after the one-way boundary treatment.
        IF (lateral_source_flag .AND. t .GE. 0.0999_wp .AND. t .LE. 0.1001_wp) THEN
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_SOURCE_BOUNDARY_BUDGET', diagnostic_flag, STATUS=diagnostic_status)

            IF (diagnostic_status .EQ. 0 .AND. diagnostic_flag .EQ. '1') THEN
                source_inflow = 0.0_wp
                other_inflow = 0.0_wp

                DO l = 1, solve_interfaces_x
                    j = j_stag_x(l)
                    k = k_stag_x(l)

                    IF (j .EQ. 1) THEN
                        face_inflow = this%G_interface_xR(1,j,k) / one_by_dy

                        IF (source_cell(1,k) .EQ. 2) THEN
                            source_inflow = source_inflow + face_inflow
                        ELSE
                            other_inflow = other_inflow + face_inflow

                        END IF
                    ELSEIF (j .EQ. comp_interfaces_x) THEN
                        face_inflow = -this%G_interface_xL(1,j,k) / one_by_dy

                        IF (source_cell(comp_cells_x,k) .EQ. 2) THEN
                            source_inflow = source_inflow + face_inflow
                        ELSE
                            other_inflow = other_inflow + face_inflow

                        END IF

                    END IF

                END DO

                DO l = 1, solve_interfaces_y
                    j = j_stag_y(l)
                    k = k_stag_y(l)

                    IF (k .EQ. 1) THEN
                        face_inflow = this%G_interface_yT(1,j,k) / one_by_dx

                        IF (source_cell(j,1) .EQ. 2) THEN
                            source_inflow = source_inflow + face_inflow
                        ELSE
                            other_inflow = other_inflow + face_inflow

                        END IF
                    ELSEIF (k .EQ. comp_interfaces_y) THEN
                        face_inflow = -this%G_interface_yB(1,j,k) / one_by_dx

                        IF (source_cell(j,comp_cells_y) .EQ. 2) THEN
                            source_inflow = source_inflow + face_inflow
                        ELSE
                            other_inflow = other_inflow + face_inflow

                        END IF

                    END IF

                END DO

                WRITE (*, '(A,ES23.15,A,ES23.15,A,ES23.15)')                                                           &
                    'SOURCE_BOUNDARY_BUDGET t=', t,                                                                    &
                    ' source_kg_s=', source_inflow,                                                                    &
                    ' other_kg_s=', other_inflow

            END IF

        END IF

    END SUBROUTINE eval_flux_PCCU

    !> \brief Evaluate both oriented PCCU values for one interface.
    !>
    !> The operation order is shared by x and y faces; the caller supplies their
    !> orientation, normal momentum slot, endpoint traces and output arrays.
    !> \param[in] eta_left Free-surface elevation at the negative-side endpoint.
    !> \param[in] eta_right Free-surface elevation at the positive-side endpoint.
    !> \param[in] bed_face Unique bed elevation shared by both states at this face.
    !> \param[in] qp_adv_left Physical face state, including advective velocities.
    !> \param[in] qp_adv_right Physical face state, including advective velocities.
    !> \param[in] grav_coeff Large-slope gravity factor at the interface.
    !> \param[in] a_minus Negative face speed bound.
    !> \param[in] a_plus Positive face speed bound.
    !> \param[in] direction Cartesian face orientation.
    !> \param[in] normal_momentum_index Normal momentum slot in physical state.
    !> \param[in] face_i X interface index.
    !> \param[in] face_j Y interface index.
    !> \param[in] use_dry_donor Enable the dry-donor face treatment.
    !> \param[in] dry_donor_left Whether the negative side uses a dry-donor trace.
    !> \param[in] dry_donor_right Whether the positive side uses a dry-donor trace.
    !> \param[in] t Current simulation or stage time [s].
    !> \param[out] value_left PCCU contribution assigned to the negative-side cell.
    !> @param[in] diagnose_front Enable opt-in diagnostic for a wet-to-dry x-face.
    !> \param[out] value_right PCCU contribution assigned to the positive-side cell.
    !> @param[in] h_selected_left Selected left-trace depth [m].
    !> @param[in] h_selected_right Selected right-trace depth [m].
    !> @param[in] gamma_left_cached Hydrostatic coefficient cached from the left CFL trace.
    !> @param[in] gamma_right_cached Hydrostatic coefficient cached from the right CFL trace.
    !> @param[in] bed_face_left Optional left-oriented face bed elevation [m].
    !> @param[in] bed_face_right Optional right-oriented face bed elevation [m].
    SUBROUTINE eval_pccu_interface(eta_left,eta_right,bed_face, h_selected_left,h_selected_right,                      &
            qp_adv_left,qp_adv_right,                                                                                  &
            grav_coeff,a_minus,a_plus,gamma_left_cached,gamma_right_cached,                                          &
            direction,normal_momentum_index,face_i,                                                                    &
            face_j,use_dry_donor,dry_donor_left,dry_donor_right,t,                                                     &
            value_left,value_right,bed_face_left,bed_face_right,diagnose_front)

    USE equation_terms_2d, ONLY : eval_inertial_flux
    USE equation_terms_2d, ONLY : limit_component_mass_flux
    USE state_conversions, ONLY : primitive_to_conservative
    USE diagnostics_2d, ONLY : fatal_error
    USE pccu_2d, ONLY : eval_hydrostatic_path
    USE pccu_2d, ONLY : eval_oriented_pccu_pair
    USE pccu_2d, ONLY : eval_pccu_advective_flux
    USE nonconservative_2d, ONLY : PATH_DIR_X

    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: eta_left, eta_right, bed_face
        REAL(wp), INTENT(IN), OPTIONAL :: bed_face_left,bed_face_right
        REAL(wp), INTENT(IN) :: h_selected_left,h_selected_right
        REAL(wp), INTENT(IN) :: qp_adv_left(n_vars), qp_adv_right(n_vars)
        REAL(wp), INTENT(IN) :: grav_coeff, a_minus, a_plus
        REAL(wp), INTENT(IN) :: gamma_left_cached,gamma_right_cached
        INTEGER, INTENT(IN) :: direction, normal_momentum_index
        INTEGER, INTENT(IN) :: face_i, face_j
        LOGICAL, INTENT(IN) :: use_dry_donor
        LOGICAL, INTENT(IN) :: dry_donor_left,dry_donor_right
        LOGICAL, INTENT(IN), OPTIONAL :: diagnose_front
        REAL(wp), INTENT(IN) :: t
        REAL(wp), INTENT(OUT) :: value_left(n_eqns), value_right(n_eqns)

        REAL(wp) :: q_adv_left_use(n_vars), q_adv_right_use(n_vars)
        REAL(wp) :: qp_adv_left_use(n_vars), qp_adv_right_use(n_vars)
        REAL(wp) :: flux_left(n_eqns), flux_right(n_eqns)
        REAL(wp) :: path_contribution(n_eqns), interface_flux(n_eqns)
        REAL(wp) :: qp_path_left(n_vars), qp_path_right(n_vars)
        REAL(wp) :: gamma_left, gamma_right
        REAL(wp) :: eta_scale, eta_roundoff
        REAL(wp) :: u_normal_left,u_normal_right
        REAL(wp) :: eta_path_left, eta_path_right
        REAL(wp) :: h_eta_left,h_eta_right
        LOGICAL :: source_trace_left,source_trace_right
        LOGICAL :: face_conversion_ok,trace_front
        ! Source-side flags remain necessary for hydrostatic path endpoints.

        trace_front=.FALSE.
        IF (PRESENT(diagnose_front)) trace_front=diagnose_front

        qp_adv_left_use=qp_adv_left
        qp_adv_right_use=qp_adv_right
        CALL lateral_source_trace_sides(direction,face_i,face_j, source_trace_left,source_trace_right)
        h_eta_left=h_selected_left
        h_eta_right=h_selected_right
        ! Reconstruct u,v separately at cell centres and form the face
        ! momenta as h_face*u_face using the selected shared depth traces.
        qp_adv_left_use(1)=h_eta_left
        qp_adv_right_use(1)=h_eta_right

        IF ( h_eta_left .LE. dry_thickness_tolerance ) THEN
            qp_adv_left_use(1:3)=0.0_wp
            qp_adv_left_use(idx_u:idx_v)=0.0_wp
        ELSE


        END IF

        IF ( h_eta_right .LE. dry_thickness_tolerance ) THEN
            qp_adv_right_use(1:3)=0.0_wp
            qp_adv_right_use(idx_u:idx_v)=0.0_wp
        ELSE


        END IF

        CALL primitive_to_conservative(qp_adv_left_use,q_adv_left_use,face_conversion_ok)
        IF (.NOT.face_conversion_ok) CALL fatal_error('KT interface left advective trace: invalid composition')
        CALL primitive_to_conservative(qp_adv_right_use,q_adv_right_use,face_conversion_ok)
        IF (.NOT.face_conversion_ok) CALL fatal_error('KT interface right advective trace: invalid composition')

        CALL eval_inertial_flux(q_adv_left_use,qp_adv_left_use,direction,flux_left)
        CALL eval_inertial_flux(q_adv_right_use,qp_adv_right_use,direction,flux_right)

        ! Build jump/path states from the same eta-based face depths. Both sides use
        ! the single shared bed value at the face.
        ! The selected-depth conservative and primitive states are already ready.
        ! Avoid copying four whole face vectors merely to name the flux endpoints.

        ! The advective and jump states share the selected depth traces. Keep
        ! separate primitive states for the hydrostatic path contribution.
        ! Factor the central-upwind transport flux into wave-speed differences
        ! times endpoint states. This is algebraically identical to the usual
        ! flux formula, but avoids subtracting nearly equal convective and
        ! dissipative terms for a passively transported quantity.
        IF ( direction .EQ. PATH_DIR_X ) THEN
            u_normal_left=qp_adv_left_use(idx_u)
            u_normal_right=qp_adv_right_use(idx_u)
        ELSE
            u_normal_left=qp_adv_left_use(idx_v)
            u_normal_right=qp_adv_right_use(idx_v)

        END IF
        ! eval_inertial_flux returns a zero endpoint flux below this same
        ! threshold; represent that endpoint consistently in the factorized
        ! expression even if its conservative trace is nonzero by roundoff.
        IF ( qp_adv_left_use(1) .LE. EPSILON(1.0_wp) ) u_normal_left=0.0_wp
        IF ( qp_adv_right_use(1) .LE. EPSILON(1.0_wp) ) u_normal_right=0.0_wp

        CALL eval_pccu_advective_flux(a_minus,a_plus,u_normal_left,u_normal_right,                                  &
            q_adv_left_use,q_adv_right_use,flux_left,flux_right,interface_flux)
        IF (trace_front) THEN
            WRITE(*,'(A,4(ES24.16,1X))') 'FRONT_KERNEL input qL,qR,uL,uR=', &
                q_adv_left_use(1),q_adv_right_use(1),u_normal_left,u_normal_right
            WRITE(*,'(A,4(ES24.16,1X))') 'FRONT_KERNEL speeds,endpoint_fluxL,fluxR=', &
                a_minus,a_plus,flux_left(1),flux_right(1)
            WRITE(*,'(A,ES24.16)') 'FRONT_KERNEL after_pccu=',interface_flux(1)
            WRITE(*,'(A,3(L1,1X))') 'FRONT_KERNEL dry_flags enable,left,right=', &
                use_dry_donor,dry_donor_left,dry_donor_right
        END IF

        CALL limit_component_mass_flux(interface_flux)
        IF (trace_front) WRITE(*,'(A,ES24.16)') 'FRONT_KERNEL after_component_limiter=',interface_flux(1)

        ! Positivity guard: when the total mass flux points out of a dry donor,
        ! suppress the complete inertial flux vector, not only the mass equation.
        ! Keep the nonconservative hydrostatic path contribution below untouched.
        IF ( use_dry_donor ) THEN
            IF ( (interface_flux(1) .GT. 0.0_wp) .AND.                                                                 &
                (dry_donor_left .OR.                                                                                   &
                (qp_adv_left_use(1) .LE. dry_thickness_tolerance)) ) THEN
                IF (trace_front) WRITE(*,'(A)') 'FRONT_KERNEL suppressed_by_dry_left'
                interface_flux=0.0_wp
            ELSEIF ( (interface_flux(1) .LT. 0.0_wp) .AND.                                                             &
                (dry_donor_right .OR.                                                                                  &
                (qp_adv_right_use(1) .LE. dry_thickness_tolerance)) ) THEN
                IF (trace_front) WRITE(*,'(A)') 'FRONT_KERNEL suppressed_by_dry_right'
                interface_flux=0.0_wp

            END IF

        END IF

        IF (trace_front) WRITE(*,'(A,ES24.16)') 'FRONT_KERNEL after_dry_guard=',interface_flux(1)

        ! Suppress only the roundoff-level mass diffusion of a hydrostatic rest
        ! state. Zero normal velocity alone is not enough: on a flat bed, for
        ! example, a wet cell next to dry cells has unequal free-surface levels
        ! and must be allowed to drain. At an excavation crest, HP reconstruction
        ! makes both states dry (so no uphill flux is introduced), while equal-eta
        ! resting traces remain protected from a spurious central-upwind flux.
        eta_scale = MAX(1.0_wp,ABS(eta_left),ABS(eta_right))
        eta_roundoff = 8.0_wp*EPSILON(1.0_wp)*eta_scale
        IF ( ( qp_adv_left(normal_momentum_index) .EQ. 0.0_wp ) .AND.                                                  &
            ( qp_adv_right(normal_momentum_index) .EQ. 0.0_wp ) .AND.                                                  &
            ( ABS(eta_left-eta_right) .LE. eta_roundoff ) ) THEN
            IF (trace_front) WRITE(*,'(A)') 'FRONT_KERNEL suppressed_by_rest_guard'
            interface_flux(1) = 0.0_wp
            interface_flux(4:n_vars) = 0.0_wp

        END IF

        IF (trace_front) THEN
            WRITE(*,'(A,3(ES24.16,1X))') 'FRONT_KERNEL after_rest_guard,eta_diff,tolerance=', &
                interface_flux(1),eta_left-eta_right,eta_roundoff
        END IF

        qp_path_left=qp_adv_left_use
        qp_path_right=qp_adv_right_use
        eta_path_left=eta_left
        eta_path_right=eta_right
        IF (PRESENT(bed_face_left)) THEN
            IF (source_trace_left) eta_path_left=bed_face_left+h_eta_left
        ELSE
            IF ( source_trace_left ) eta_path_left=bed_face+h_eta_left
        END IF
        IF (PRESENT(bed_face_right)) THEN
            IF (source_trace_right) eta_path_right=bed_face_right+h_eta_right
        ELSE
            IF ( source_trace_right ) eta_path_right=bed_face+h_eta_right
        END IF

        IF ( use_dry_donor ) THEN

            IF ( qp_adv_left_use(1) .LE. dry_thickness_tolerance ) THEN
                qp_path_left(1)=0.0_wp
                IF (PRESENT(bed_face_left)) THEN
                    eta_path_left=bed_face_left
                ELSE
                    eta_path_left=bed_face
                END IF

            END IF

            IF ( qp_adv_right_use(1) .LE. dry_thickness_tolerance ) THEN
                qp_path_right(1)=0.0_wp
                IF (PRESENT(bed_face_right)) THEN
                    eta_path_right=bed_face_right
                ELSE
                    eta_path_right=bed_face
                END IF

            END IF

        END IF

        ! The characteristic-speed pass already evaluated rho_m*g' from the
        ! identical finalized temperature and component fractions at this face.
        ! Only the dry-donor path may change endpoint depth to zero; reproduce
        ! exactly the hydrostatic kernel's zero-depth branch in that case.
        gamma_left=gamma_left_cached
        gamma_right=gamma_right_cached
        IF (qp_path_left(1)<=EPSILON(1.0_wp)) gamma_left=0.0_wp
        IF (qp_path_right(1)<=EPSILON(1.0_wp)) gamma_right=0.0_wp
        CALL regularize_dry_gamma_pair(qp_path_left(1),gamma_left, qp_path_right(1),gamma_right)
        CALL eval_hydrostatic_path(direction,qp_path_left(1),gamma_left,eta_path_left,                                 &
            grav_coeff,qp_path_right(1),gamma_right,eta_path_right,grav_coeff,                                         &
            path_contribution)

        CALL eval_oriented_pccu_pair(interface_flux,path_contribution,a_minus, a_plus,value_left,value_right)
        IF (trace_front) THEN
            WRITE(*,'(A,3(ES24.16,1X))') 'FRONT_KERNEL oriented flux,path,left=', &
                interface_flux(1),path_contribution(1),value_left(1)
            WRITE(*,'(A,ES24.16)') 'FRONT_KERNEL oriented_right=',value_right(1)
        END IF

        ! A dry adjacent column receives no share of the nonconservative path;
        ! the wet side retains its own oriented contribution.
        IF ( use_dry_donor ) THEN
            IF ( qp_adv_left_use(1) .LE. dry_thickness_tolerance ) value_left=interface_flux
            IF ( qp_adv_right_use(1) .LE. dry_thickness_tolerance ) value_right=interface_flux

        END IF
        IF (trace_front) THEN
            WRITE(*,'(A,2(ES24.16,1X))') 'FRONT_KERNEL final left,right=',value_left(1),value_right(1)
        END IF

    END SUBROUTINE eval_pccu_interface

    !******************************************************************************
    !> \brief Use the wet endpoint pressure coefficient on a wet/dry path.
    !>
    !> Gamma is a material coefficient, not a thickness variable.  A dry final
    !> trace carries no composition and eval_hydrostatic_coefficient therefore
    !> returns Gamma=0.  Along a wet/dry shoreline this artificial jump would
    !> create a spurious -0.5*h^2*dGamma contribution in the hydrostatic path.
    !> Use the wet-side limiting value at the dry endpoint instead.
    !>
    !> \param[in] h_left Depth at the negative/left path endpoint [m].
    !> \param[in,out] gamma_left Left coefficient, replaced by the right one if only the left trace
    !> is
    !>                           dry.
    !> \param[in] h_right Depth at the positive/right path endpoint [m].
    !> \param[in,out] gamma_right Right coefficient, replaced by the left one if only the right
    !> trace
    !>                            is dry.
    !******************************************************************************

    !> @brief Extend the wet-side hydrostatic coefficient across a dry interface.
    PURE SUBROUTINE regularize_dry_gamma_pair(h_left,gamma_left,h_right,gamma_right)
        REAL(wp), INTENT(IN) :: h_left,h_right
        REAL(wp), INTENT(INOUT) :: gamma_left,gamma_right
        IF ( ( h_left .LE. dry_thickness_tolerance ) .AND.                                                             &
            ( h_right .GT. dry_thickness_tolerance ) ) THEN
            gamma_left = gamma_right
        ELSEIF ( ( h_right .LE. dry_thickness_tolerance ) .AND.                                                        &
            ( h_left .GT. dry_thickness_tolerance ) ) THEN
            gamma_right = gamma_left

        END IF
    END SUBROUTINE regularize_dry_gamma_pair

    !******************************************************************************
    !> \brief Identify lateral source traces on the domain boundary.
    !> \param[in] direction Cartesian orientation of the face.
    !> \param[in] face_i X interface index.
    !> \param[in] face_j Y interface index.
    !> \param[out] source_left True when the negative-side trace is a source cell.
    !> \param[out] source_right True when the positive-side trace is a source cell.
    !******************************************************************************
    SUBROUTINE lateral_source_trace_sides(direction,face_i,face_j, source_left,source_right)

    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

    IMPLICIT NONE

        INTEGER, INTENT(IN) :: direction,face_i,face_j
        LOGICAL, INTENT(OUT) :: source_left,source_right

        source_left=.FALSE.
        source_right=.FALSE.
        IF ( .NOT. lateral_source_flag ) RETURN

        SELECT CASE (direction)
        CASE (PATH_DIR_X)
            IF ( face_i .EQ. 1 ) source_left=source_cell(1,face_j) .EQ. 2
            IF ( face_i .EQ. comp_interfaces_x )                                                                       &
                source_right=source_cell(comp_cells_x,face_j) .EQ. 2
        CASE (PATH_DIR_Y)
            IF ( face_j .EQ. 1 ) source_left=source_cell(face_i,1) .EQ. 2
            IF ( face_j .EQ. comp_interfaces_y )                                                                       &
                source_right=source_cell(face_i,comp_cells_y) .EQ. 2

        END SELECT

    END SUBROUTINE lateral_source_trace_sides

    !******************************************************************************
    !> \brief Integrate the hydrostatic path between the final traces of each active cell.
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
    !> \param[in] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
    !> \param[in] qp_expl Cell-centered physical state at the explicit stage.
    !> \param[in] use_dry_donor Enable the dry-donor face treatment.
    !> \param[in] solve_cells Number of entries in the active-cell index lists.
    !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
    !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
    !******************************************************************************

    !> @brief Integrate within-cell nonconservative hydrostatic paths.
    !> @param[in] lower_recon Lower-layer reconstruction.
    SUBROUTINE eval_cell_hydrostatic_paths( this, recon, qp_expl, use_dry_donor, solve_cells, j_cent, k_cent, &
        lower_recon )

    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

    IMPLICIT NONE

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        REAL(wp), INTENT(IN) :: qp_expl(n_vars,comp_cells_x,comp_cells_y)
        LOGICAL, INTENT(IN) :: use_dry_donor
        INTEGER, INTENT(IN) :: solve_cells
        INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
        CLASS(reconstruction_workspace_type), INTENT(IN), OPTIONAL :: lower_recon

        INTEGER :: j, k, l
        REAL(wp) :: qp_path_minus(n_vars), qp_path_plus(n_vars)

        !$OMP PARALLEL DO private(l,j,k,qp_path_minus,qp_path_plus)
        DO l = 1, solve_cells

            j = j_cent(l)
            k = k_cent(l)

            IF ( use_dry_donor .AND. (qp_expl(1,j,k) .LE. dry_thickness_tolerance) ) THEN
                this%P_cell_x(:,j,k)=0.0_wp
                this%P_cell_y(:,j,k)=0.0_wp
                CYCLE

            END IF

            IF ( comp_cells_x .GT. 1 ) THEN
                qp_path_minus=recon%qp_adv_interfaceR(:,j,k)
                qp_path_plus=recon%qp_adv_interfaceL(:,j+1,k)
                IF (PRESENT(lower_recon)) THEN
                    qp_path_minus(1)=MAX(recon%eta_cellW(j,k)-lower_recon%eta_cellW(j,k),0.0_wp)
                ELSE
                qp_path_minus(1)=MAX(recon%eta_cellW(j,k)-B_face_x(j,k),0.0_wp)
                END IF
                IF (PRESENT(lower_recon)) THEN
                    qp_path_plus(1)=MAX(recon%eta_cellE(j,k)-lower_recon%eta_cellE(j,k),0.0_wp)
                ELSE
                qp_path_plus(1)=MAX(recon%eta_cellE(j,k)-B_face_x(j+1,k),0.0_wp)
                END IF
                IF (qp_path_minus(1).LE.dry_thickness_tolerance) qp_path_minus(2:3)=0.0_wp
                IF (qp_path_plus(1).LE.dry_thickness_tolerance) qp_path_plus(2:3)=0.0_wp
                CALL eval_cell_hydrostatic_path( qp_path_minus, recon%eta_cellW(j,k),grav_coeff_stag_x(j,k),           &
                    qp_path_plus,recon%eta_cellE(j,k),                                                                 &
                    grav_coeff_stag_x(j+1,k),PATH_DIR_X,this%P_cell_x(:,j,k),                                      &
                    this%gamma_interface_xR(j,k),recon%h_selected_interfaceR(j,k),                                   &
                    this%gamma_interface_xL(j+1,k),recon%h_selected_interfaceL(j+1,k),                            &
                    this%clamp_buoyant_pressure )

            END IF

            IF ( comp_cells_y .GT. 1 ) THEN
                qp_path_minus=recon%qp_adv_interfaceT(:,j,k)
                qp_path_plus=recon%qp_adv_interfaceB(:,j,k+1)
                IF (PRESENT(lower_recon)) THEN
                    qp_path_minus(1)=MAX(recon%eta_cellS(j,k)-lower_recon%eta_cellS(j,k),0.0_wp)
                ELSE
                qp_path_minus(1)=MAX(recon%eta_cellS(j,k)-B_face_y(j,k),0.0_wp)
                END IF
                IF (PRESENT(lower_recon)) THEN
                    qp_path_plus(1)=MAX(recon%eta_cellN(j,k)-lower_recon%eta_cellN(j,k),0.0_wp)
                ELSE
                qp_path_plus(1)=MAX(recon%eta_cellN(j,k)-B_face_y(j,k+1),0.0_wp)
                END IF
                IF (qp_path_minus(1).LE.dry_thickness_tolerance) qp_path_minus(2:3)=0.0_wp
                IF (qp_path_plus(1).LE.dry_thickness_tolerance) qp_path_plus(2:3)=0.0_wp
                CALL eval_cell_hydrostatic_path( qp_path_minus, recon%eta_cellS(j,k),grav_coeff_stag_y(j,k),           &
                    qp_path_plus,recon%eta_cellN(j,k),                                                                 &
                    grav_coeff_stag_y(j,k+1),PATH_DIR_Y,this%P_cell_y(:,j,k),                                      &
                    this%gamma_interface_yT(j,k),recon%h_selected_interfaceT(j,k),                                   &
                    this%gamma_interface_yB(j,k+1),recon%h_selected_interfaceB(j,k+1),                            &
                    this%clamp_buoyant_pressure )

            END IF

        END DO
        !$OMP END PARALLEL DO

    END SUBROUTINE eval_cell_hydrostatic_paths

    !> \brief Compute the hydrostatic path through one cell in either direction.
    !>
    !> Endpoint ordering is negative-to-positive: west-to-east in x and
    !> south-to-north in y. The cell-specific trace and gravity mapping stays in
    !> the caller, alongside its direction-specific output array.
    !> \param[in] qp_minus Physical state at the negative path endpoint.
    !> \param[in] eta_minus Free-surface elevation at the negative path endpoint.
    !> \param[in] grav_minus Large-slope gravity factor at the negative endpoint.
    !> \param[in] qp_plus Physical state at the positive path endpoint.
    !> \param[in] eta_plus Free-surface elevation at the positive path endpoint.
    !> \param[in] grav_plus Large-slope gravity factor at the positive endpoint.
    !> \param[in] direction Cartesian path orientation.
    !> \param[out] path Hydrostatic path contribution across the cell.
    !> @param[in] cached_gamma_minus Face pressure coefficient from the CFL pass.
    !> @param[in] selected_h_minus Final face depth used for its cached coefficient.
    !> @param[in] cached_gamma_plus Corresponding positive-side pressure coefficient.
    !> @param[in] selected_h_plus Final face depth of the positive side.
    !> @param[in] clamp_buoyant_pressure Suppress negative pressure at buoyant L2 endpoints.
    !> Optional resolved pressure coefficients.
    !> @param[out] gamma_minus_out Optional negative-endpoint hydrostatic coefficient.
    !> @param[out] gamma_plus_out Optional positive-endpoint hydrostatic coefficient.
    SUBROUTINE eval_cell_hydrostatic_path(qp_minus,eta_minus,grav_minus, qp_plus,eta_plus,grav_plus,direction,path, &
        cached_gamma_minus,selected_h_minus,cached_gamma_plus,selected_h_plus,clamp_buoyant_pressure, &
        gamma_minus_out,gamma_plus_out)

    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    USE pccu_2d, ONLY : eval_hydrostatic_path

    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: qp_minus(:),eta_minus,grav_minus
        REAL(wp), INTENT(IN) :: qp_plus(:),eta_plus,grav_plus
        REAL(wp), INTENT(IN) :: cached_gamma_minus,selected_h_minus,cached_gamma_plus,selected_h_plus
        LOGICAL, INTENT(IN) :: clamp_buoyant_pressure
        INTEGER, INTENT(IN) :: direction
        REAL(wp), INTENT(OUT) :: path(:)
        REAL(wp), INTENT(OUT), OPTIONAL :: gamma_minus_out,gamma_plus_out

        REAL(wp) :: gamma_minus,gamma_plus
        REAL(wp) :: reduced_gravity_minus,reduced_gravity_plus

        ! Face speed coefficients were evaluated on the same thermodynamic
        ! primitive traces. Their selected depths may differ from the path
        ! depths, so reuse only when both endpoints are genuinely wet.
        ! Otherwise retain the original density evaluation and shoreline
        ! regularization. No tolerance-based comparison of compositions.
        IF (qp_minus(1)>EPSILON(1.0_wp) .AND. selected_h_minus>EPSILON(1.0_wp)) THEN
            gamma_minus=cached_gamma_minus
        ELSE
            CALL eval_hydrostatic_coefficient(qp_minus,reduced_gravity_minus,gamma_minus)
        END IF
        IF (qp_plus(1)>EPSILON(1.0_wp) .AND. selected_h_plus>EPSILON(1.0_wp)) THEN
            gamma_plus=cached_gamma_plus
        ELSE
            CALL eval_hydrostatic_coefficient(qp_plus,reduced_gravity_plus,gamma_plus)
        END IF
        IF (clamp_buoyant_pressure) THEN
            gamma_minus=MAX(gamma_minus,0.0_wp)
            gamma_plus=MAX(gamma_plus,0.0_wp)
        END IF
        CALL regularize_dry_gamma_pair(qp_minus(1),gamma_minus,qp_plus(1), gamma_plus)
        IF (PRESENT(gamma_minus_out)) gamma_minus_out=gamma_minus
        IF (PRESENT(gamma_plus_out)) gamma_plus_out=gamma_plus
        CALL eval_hydrostatic_path(direction,qp_minus(1),gamma_minus,eta_minus,                                        &
            grav_minus,qp_plus(1),gamma_plus,eta_plus,grav_plus,path)

    END SUBROUTINE eval_cell_hydrostatic_path

    !> \brief Release the hyperbolic workspace arrays.
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.

    SUBROUTINE finalize_hyperbolic( this )

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this

        DEALLOCATE( this%a_interface_xNeg )
        DEALLOCATE( this%a_interface_xPos )
        DEALLOCATE( this%a_interface_yNeg )
        DEALLOCATE( this%a_interface_yPos )
        DEALLOCATE( this%gamma_interface_xL )
        DEALLOCATE( this%gamma_interface_xR )
        DEALLOCATE( this%gamma_interface_yB )
        DEALLOCATE( this%gamma_interface_yT )
        DEALLOCATE( this%G_interface_xL )
        DEALLOCATE( this%G_interface_xR )
        DEALLOCATE( this%G_interface_yB )
        DEALLOCATE( this%G_interface_yT )
        DEALLOCATE( this%P_cell_x )
        DEALLOCATE( this%P_cell_y )

    END SUBROUTINE finalize_hyperbolic

    !******************************************************************************
    !> \brief Evaluate facewise characteristic bounds from the reconstructed endpoint states.
    !
    !> This subroutine evaluates the largest characteristic speed at the
    !> cells interfaces from the reconstructed states.
    !> @author
    !> Mattia de' Michieli Vitturi
    !> \date 2019/11/11
    !>
    !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
    !> \param[in] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
    !> \param[in] solve_interfaces_x Number of active x-normal faces.
    !> \param[in] j_stag_x X-face indices paired with k_stag_x.
    !> \param[in] k_stag_x Y-cell indices of active x-normal faces.
    !> \param[in] solve_interfaces_y Number of active y-normal faces.
    !> \param[in] j_stag_y X-cell indices of active y-normal faces.
    !> \param[in] k_stag_y Y-face indices paired with j_stag_y.
    !******************************************************************************

    !> @brief Compute characteristic wave-speed bounds at reconstructed faces.
    SUBROUTINE eval_speeds( this, recon, solve_interfaces_x, j_stag_x,                                                 &
            k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

        ! External procedures
    USE equation_terms_2d, ONLY : eval_local_speeds_direction

    IMPLICIT NONE

        CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
        CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
        INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
        INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
        INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

        REAL(wp) :: min_r, max_r
        REAL(wp) :: h_eta_minus, h_eta_plus
        REAL(wp) :: qp_eta_minus(n_vars), qp_eta_plus(n_vars)

        INTEGER :: j,k,l

        !$OMP PARALLEL

        IF ( comp_cells_x .GT. 1 ) THEN

            !$OMP DO private(j,k,min_r,max_r,h_eta_minus,h_eta_plus, &
            !$OMP& qp_eta_minus,qp_eta_plus)

            x_interfaces_loop:DO l = 1,solve_interfaces_x

                j = j_stag_x(l)
                k = k_stag_x(l)
                ! Selected-depth endpoint states are the only wave-speed inputs.
                h_eta_minus=recon%h_selected_interfaceL(j,k)
                h_eta_plus=recon%h_selected_interfaceR(j,k)
                qp_eta_minus=recon%qp_adv_interfaceL(:,j,k)
                qp_eta_plus=recon%qp_adv_interfaceR(:,j,k)
                qp_eta_minus(1)=h_eta_minus
                qp_eta_plus(1)=h_eta_plus
                IF (h_eta_minus .LE. dry_thickness_tolerance)                                                          &
                    qp_eta_minus(idx_u:idx_v)=0.0_wp
                IF (h_eta_plus .LE. dry_thickness_tolerance)                                                           &
                    qp_eta_plus(idx_u:idx_v)=0.0_wp


                CALL eval_interface_speed_bounds(qp_eta_minus,qp_eta_plus,grav_coeff_stag_x(j,k),idx_u,min_r,max_r, &
                    this%gamma_interface_xL(j,k),this%gamma_interface_xR(j,k),this%clamp_buoyant_pressure)

                this%a_interface_xNeg(j,k) = min_r
                this%a_interface_xPos(j,k) = max_r

            END DO x_interfaces_loop

            !$OMP END DO NOWAIT

        END IF

        IF ( comp_cells_y .GT. 1 ) THEN

            !$OMP DO private(j,k,min_r,max_r,h_eta_minus,h_eta_plus, &
            !$OMP& qp_eta_minus,qp_eta_plus)

            y_interfaces_loop:DO l = 1,solve_interfaces_y

                j = j_stag_y(l)
                k = k_stag_y(l)
                ! Selected-depth endpoint states are the only wave-speed inputs.
                h_eta_minus=recon%h_selected_interfaceB(j,k)
                h_eta_plus=recon%h_selected_interfaceT(j,k)
                qp_eta_minus=recon%qp_adv_interfaceB(:,j,k)
                qp_eta_plus=recon%qp_adv_interfaceT(:,j,k)
                qp_eta_minus(1)=h_eta_minus
                qp_eta_plus(1)=h_eta_plus
                IF (h_eta_minus .LE. dry_thickness_tolerance)                                                          &
                    qp_eta_minus(idx_u:idx_v)=0.0_wp
                IF (h_eta_plus .LE. dry_thickness_tolerance)                                                           &
                    qp_eta_plus(idx_u:idx_v)=0.0_wp


                CALL eval_interface_speed_bounds(qp_eta_minus,qp_eta_plus,grav_coeff_stag_y(j,k),idx_v,min_r,max_r, &
                    this%gamma_interface_yB(j,k),this%gamma_interface_yT(j,k),this%clamp_buoyant_pressure)

                this%a_interface_yNeg(j,k) = min_r
                this%a_interface_yPos(j,k) = max_r

            END DO y_interfaces_loop

            !$OMP END DO

        END IF

        !$OMP END PARALLEL

        RETURN

    END SUBROUTINE eval_speeds

    !> \brief Reduce the characteristic bounds from both endpoints of one face.
    !>
    !> The caller supplies the face-normal velocity component and keeps its
    !> x/y-specific face list and output storage separate.
    !> \param[in] qp_minus Physical state at the negative-side face endpoint.
    !> \param[in] qp_plus Physical state at the positive-side face endpoint.
    !> \param[in] grav_coeff Large-slope gravity factor at the face.
    !> \param[in] normal_velocity_index Index of the face-normal velocity.
    !> \param[out] speed_minus Minimum wave speed, bounded above by zero.
    !> \param[out] speed_plus Maximum wave speed, bounded below by zero.
    !> Resolve the configuration string outside all active-interface loops.


    !> @brief Combine endpoint characteristic bounds for one oriented interface.
    !> @param[out] gamma_minus Hydrostatic coefficient of the negative trace.
    !> @param[out] gamma_plus Hydrostatic coefficient of the positive trace.
    !> @param[in] clamp_buoyant_pressure Suppress negative L2 hydrostatic coefficients.
    SUBROUTINE eval_interface_speed_bounds(qp_minus,qp_plus,grav_coeff, normal_velocity_index,speed_minus,speed_plus,  &
        gamma_minus,gamma_plus,clamp_buoyant_pressure)

    USE equation_terms_2d, ONLY : eval_local_speeds_direction

    IMPLICIT NONE

        REAL(wp), INTENT(IN) :: qp_minus(:),qp_plus(:),grav_coeff
        INTEGER, INTENT(IN) :: normal_velocity_index
        REAL(wp), INTENT(OUT) :: speed_minus,speed_plus
        REAL(wp), INTENT(OUT) :: gamma_minus,gamma_plus
        LOGICAL, INTENT(IN) :: clamp_buoyant_pressure
        REAL(wp) :: minus_min,minus_max,plus_min,plus_max

        CALL eval_local_speeds_direction(qp_minus,grav_coeff,normal_velocity_index,minus_min,minus_max,gamma_minus)
        CALL eval_local_speeds_direction(qp_plus,grav_coeff,normal_velocity_index,plus_min,plus_max,gamma_plus)

        ! The celerity already vanishes for negative g'; the cached pressure
        ! coefficient must obey the same L2 clamp at both face and cell paths.
        IF (clamp_buoyant_pressure) THEN
            gamma_minus=MAX(gamma_minus,0.0_wp)
            gamma_plus=MAX(gamma_plus,0.0_wp)
        END IF

        speed_minus = MIN(minus_min,plus_min,0.0_wp)
        speed_plus = MAX(minus_max,plus_max,0.0_wp)

    END SUBROUTINE eval_interface_speed_bounds

END MODULE hyperbolic_2d
