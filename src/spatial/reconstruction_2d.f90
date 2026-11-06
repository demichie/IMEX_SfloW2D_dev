!********************************************************************************
!> \brief Reconstruction of cell states at computational interfaces
!>
!> This module owns the reconstructed conservative/physical interface states
!> together with free-surface traces and reusable HP scratch storage.
!>
!> West/east/south/north arrays are cell-owned traces; left/right/bottom/top arrays are
!> face-oriented traces. HP evaluation retains separate candidate and commit phases to avoid reading
!> partially updated neighbours.
!********************************************************************************

MODULE reconstruction_2d

    USE parameters_2d, ONLY : idx_h
    USE parameters_2d, ONLY : wp, n_vars, kt_trace_reconstruction, dry_thickness_tolerance
    USE parameters_2d, ONLY : lateral_source_flag, radial_source_flag
    USE parameters_2d, ONLY : bc, bcW, bcE, bcS, bcN
    USE parameters_2d, ONLY : idx_u, idx_v
    USE parameters_2d, ONLY : limiter, reconstr_coeff

    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
    USE geometry_2d, ONLY : B_cent, B_face_x, B_face_y
    USE geometry_2d, ONLY : source_cell
    USE geometry_2d, ONLY : one_by_dx, one_by_dy, dx, dy
    USE geometry_2d, ONLY : limit
    USE geometry_2d, ONLY : x_comp, y_comp, dx2, dy2

    USE constitutive_parameters_2d, ONLY : grav

    IMPLICIT NONE

    PRIVATE

    ! Smooth transition from eta_x to h_x+B_x over this directional Froude band.
    REAL(wp), PARAMETER :: fr_blend_lower = 0.8_wp
    REAL(wp), PARAMETER :: fr_blend_upper = 1.2_wp

    !> \brief Persistent cell and face traces for the selected depth reconstruction.
    TYPE, PUBLIC :: reconstruction_workspace_type


        ! Direct h_x-based states for the advective part of the PCCU flux.
        ! The existing interface states retain the eta-based reconstruction and
        ! are used for central-upwind dissipation and hydrostatic path terms.
        REAL(wp), ALLOCATABLE :: qp_adv_interfaceL(:,:,:)
        REAL(wp), ALLOCATABLE :: h_selected_interfaceL(:,:)
        REAL(wp), ALLOCATABLE :: qp_adv_interfaceR(:,:,:)
        REAL(wp), ALLOCATABLE :: h_selected_interfaceR(:,:)
        REAL(wp), ALLOCATABLE :: qp_adv_interfaceB(:,:,:)
        REAL(wp), ALLOCATABLE :: h_selected_interfaceB(:,:)
        REAL(wp), ALLOCATABLE :: qp_adv_interfaceT(:,:,:)
        REAL(wp), ALLOCATABLE :: h_selected_interfaceT(:,:)

    ! Direct h_x-based cell-side candidates retained for the advective flux
    ! while the eta-based HP states are used by dissipation and path terms.
        REAL(wp), ALLOCATABLE :: qp_cellW(:,:,:)
        REAL(wp), ALLOCATABLE :: qp_cellE(:,:,:)
        REAL(wp), ALLOCATABLE :: qp_cellS(:,:,:)
        REAL(wp), ALLOCATABLE :: qp_cellN(:,:,:)

        ! Final free-surface traces.  Both cell-side and face-oriented forms are
        ! stored because the following PCCU path discretization needs both views.
        REAL(wp), ALLOCATABLE :: eta_cellW(:,:), eta_cellE(:,:)
        REAL(wp), ALLOCATABLE :: eta_cellS(:,:), eta_cellN(:,:)
        REAL(wp), ALLOCATABLE :: eta_interfaceL(:,:), eta_interfaceR(:,:)
        REAL(wp), ALLOCATABLE :: eta_interfaceB(:,:), eta_interfaceT(:,:)
        REAL(wp), ALLOCATABLE :: h_minmod_interfaceL(:,:), h_minmod_interfaceR(:,:)
        REAL(wp), ALLOCATABLE :: h_minmod_interfaceB(:,:), h_minmod_interfaceT(:,:)
        REAL(wp), ALLOCATABLE :: h_froude_interfaceL(:,:), h_froude_interfaceR(:,:)
        REAL(wp), ALLOCATABLE :: h_froude_interfaceB(:,:), h_froude_interfaceT(:,:)

            LOGICAL, ALLOCATABLE :: hp_eta_mask(:,:)
            INTEGER, ALLOCATABLE :: hp_eta_j(:), hp_eta_k(:)
            INTEGER :: hp_eta_cells
            LOGICAL :: hp_workset_prepared = .FALSE.
            ! Per direction: h-/h+, normal momentum-/+, then tangential momentum-/+.

    CONTAINS
        PROCEDURE :: initialize => initialize_reconstruction
        PROCEDURE :: finalize => finalize_reconstruction
        PROCEDURE :: reconstruct => reconstruction
        PROCEDURE :: prepare_hp_workset
        PROCEDURE :: reconstruct_kt_eta => reconstruct_kt_eta_traces
    END TYPE reconstruction_workspace_type

CONTAINS

    !> \brief Allocate interface traces, HP candidates and per-thread scratch storage.
    !>
    !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.

    SUBROUTINE initialize_reconstruction( this )

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this



        ALLOCATE( this%qp_adv_interfaceL( n_vars, comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_selected_interfaceL( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%qp_adv_interfaceR( n_vars, comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_selected_interfaceR( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%qp_adv_interfaceB( n_vars, comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_selected_interfaceB( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%qp_adv_interfaceT( n_vars, comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_selected_interfaceT( comp_cells_x, comp_interfaces_y ) )

        ALLOCATE( this%qp_cellW( n_vars, comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%qp_cellE( n_vars, comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%qp_cellS( n_vars, comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%qp_cellN( n_vars, comp_cells_x, comp_cells_y ) )

        ALLOCATE( this%eta_cellW( comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%eta_cellE( comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%eta_cellS( comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%eta_cellN( comp_cells_x, comp_cells_y ) )
        ALLOCATE( this%eta_interfaceL( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%eta_interfaceR( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%eta_interfaceB( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%eta_interfaceT( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_minmod_interfaceL( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_minmod_interfaceR( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_minmod_interfaceB( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_minmod_interfaceT( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_froude_interfaceL( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_froude_interfaceR( comp_interfaces_x, comp_cells_y ) )
        ALLOCATE( this%h_froude_interfaceB( comp_cells_x, comp_interfaces_y ) )
        ALLOCATE( this%h_froude_interfaceT( comp_cells_x, comp_interfaces_y ) )

        ALLOCATE( this%hp_eta_mask(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%hp_eta_j(comp_cells_x*comp_cells_y) )
        ALLOCATE( this%hp_eta_k(comp_cells_x*comp_cells_y) )
        this%hp_eta_mask = .FALSE.
        this%hp_eta_cells = 0
        this%hp_workset_prepared = .FALSE.


    END SUBROUTINE initialize_reconstruction

    !> \brief Reconstruct physical traces, apply boundary conditions and finalize HP states.
    !>
    !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
    !> \param[in] qp_expl Cell-centered physical states at the explicit stage, including the
    !> required
    !>                    stencil halo.
    !> \param[in] t Current simulation or stage time [s].
    !> \param[in] solve_cells Number of entries in the active-cell index lists.
    !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
    !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
    !> \param[in] use_kt_direct Use direct KT traces instead of the HP dry-state projection.
    !> \param[in] source_boundary_enabled Optional override for boundary-source traces.


    SUBROUTINE reconstruction( this, qp_expl, t, solve_cells, j_cent, k_cent, use_kt_direct, source_boundary_enabled )

        ! External procedures
    USE mass_boundary_injection_2d, ONLY : eval_boundary_injection_explicit
    USE parameters_2d, ONLY : limiter, theta
    USE geometry_2d, ONLY : minmod

    USE geometry_2d, ONLY : x_comp , x_stag , y_comp , y_stag , dx2 , dy2

    USE geometry_2d, ONLY : sourceW , sourceE , sourceN , sourceS
    USE geometry_2d, ONLY : sourceW_vect_x , sourceW_vect_y
    USE geometry_2d, ONLY : sourceE_vect_x , sourceE_vect_y
    USE geometry_2d, ONLY : sourceN_vect_x , sourceN_vect_y
    USE geometry_2d, ONLY : sourceS_vect_x , sourceS_vect_y

    USE parameters_2d, ONLY : reconstr_coeff

    IMPLICIT NONE

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: qp_expl(:,:,:)
        REAL(wp), INTENT(IN) :: t
        INTEGER, INTENT(IN) :: solve_cells
        INTEGER, INTENT(IN) :: j_cent(:)
        INTEGER, INTENT(IN) :: k_cent(:)
        LOGICAL, INTENT(IN) :: use_kt_direct
        LOGICAL, INTENT(IN), OPTIONAL :: source_boundary_enabled

        REAL(wp) :: qrecW(n_vars) !< recons var at the west edge of the cells
        REAL(wp) :: qrecE(n_vars) !< recons var at the east edge of the cells
        REAL(wp) :: qrecS(n_vars) !< recons var at the south edge of the cells
        REAL(wp) :: qrecN(n_vars) !< recons var at the north edge of the cells

        REAL(wp) :: source_bdry(n_vars)
        REAL(wp) :: qrec_prime_x(n_vars)      !< limited x-slopes for boundary cells
        REAL(wp) :: qrec_prime_y(n_vars)      !< limited y-slopes for boundary cells
        ! The five primitive stencil vectors are contiguous in their component index.
        ! They are private per OpenMP worker and never modify qp_expl.
        REAL(wp) :: center_state(n_vars), west_state(n_vars), east_state(n_vars)
        REAL(wp) :: south_state(n_vars), north_state(n_vars)

        TYPE(bc) :: bcW_effective, bcE_effective, bcS_effective, bcN_effective

        REAL(wp) :: qrec_stencil(3) !< recons variables stencil for the limiter
        REAL(wp) :: x_stencil(3)    !< grid stencil for the limiter
        REAL(wp) :: y_stencil(3)    !< grid stencil for the limiter

        INTEGER :: l,j,k            !< loop counters (cells)
        INTEGER :: i                !< loop counter (variables)

        REAL(wp) :: dq

        LOGICAL :: regular_interior
        LOGICAL :: use_boundary_source
        LOGICAL :: uniform_generalized_minmod

        ! Keep the same constant candidates for inactive neighbours, but initialize
        ! only cells that can contribute to an active HP trace.
        use_boundary_source=lateral_source_flag
        IF (PRESENT(source_boundary_enabled))                                                                          &
            use_boundary_source=source_boundary_enabled
        ! The configured limiter is constant across cells for this stage.
        ! A uniform generalized-minmod configuration can use a component-wise
        ! array kernel without choosing the limiter on every component.
        uniform_generalized_minmod=ALL(limiter(1:n_vars)==3)
        IF (.NOT. this%hp_workset_prepared)                                                                            &
            CALL build_hp_workset(this, solve_cells, j_cent, k_cent)

        !$OMP PARALLEL DO PRIVATE(l,j,k,i,qrecW,qrecE,qrecS,qrecN,x_stencil,y_stencil,                         &
        !$OMP & qrec_stencil,qrec_prime_x,qrec_prime_y,bcW_effective,bcE_effective,bcS_effective,bcN_effective, &
        !$OMP & source_bdry,dq,regular_interior,center_state,west_state,east_state,south_state,north_state)

        DO l = 1, this%hp_eta_cells
            j = this%hp_eta_j(l)
            k = this%hp_eta_k(l)
            qrecW = qp_expl(:,j,k)
            IF (qp_expl(1,j,k) .LE. dry_thickness_tolerance)                                                           &
                qrecW(idx_u:idx_v) = 0.0_wp
            qrecE = qrecW
            qrecS = qrecW
            qrecN = qrecW

            IF (l .LE. solve_cells) THEN
            x_stencil(2) = x_comp(j)
            y_stencil(2) = y_comp(k)

            ! Default source-side ghost state. For the conservative radial source,
            ! the ring is a wall and mass/momentum are injected through eval_expl_terms.
            ! The existing lateral boundary source keeps its Dirichlet treatment.
            source_bdry(1:n_vars) = qp_expl(1:n_vars,j,k)

            IF ( use_boundary_source .AND. lateral_source_flag .AND.                                                   &
                ( source_cell(j,k) .EQ. 2 ) ) THEN

                IF ( sourceE(j,k) ) THEN

                    CALL eval_boundary_injection_explicit( t, sourceE_vect_x(j,k), sourceE_vect_y(j,k), source_bdry )

                ELSEIF ( sourceW(j,k) ) THEN

                    CALL eval_boundary_injection_explicit( t, sourceW_vect_x(j,k), sourceW_vect_y(j,k), source_bdry )

                ELSEIF ( sourceS(j,k) ) THEN

                    CALL eval_boundary_injection_explicit( t, sourceS_vect_x(j,k), sourceS_vect_y(j,k), source_bdry )

                ELSEIF ( sourceN(j,k) ) THEN

                    CALL eval_boundary_injection_explicit( t, sourceN_vect_x(j,k), sourceN_vect_y(j,k), source_bdry )

                END IF

            END IF

            regular_interior = ( j .GT. 1 ) .AND. ( j .LT. comp_cells_x ) .AND.                                        &
                ( k .GT. 1 ) .AND. ( k .LT. comp_cells_y ) .AND.                                                       &
                ( source_cell(j,k) .NE. 2 )

            IF ( regular_interior ) THEN

                ! Interior fast path: fetch each neighbour's complete primitive state
                ! once. In Fortran the component index is contiguous, so these
                ! short vector loads avoid five indexed function calls per component.
                ! No boundary conditions or source-wall substitutions occur here.
                center_state = qp_expl(:,j,k)
                west_state   = qp_expl(:,j-1,k)
                east_state   = qp_expl(:,j+1,k)
                south_state  = qp_expl(:,j,k-1)
                north_state  = qp_expl(:,j,k+1)

                ! Reproduce reconstruction_cell_value's dry-neighbour rule:
                ! a dry cell cannot contribute an auxiliary velocity to a limiter.
                ! Other physical variables, including h, are left unchanged.
                IF (center_state(idx_h) <= dry_thickness_tolerance) center_state(idx_u:idx_v)=0.0_wp
                IF (west_state(idx_h)   <= dry_thickness_tolerance) west_state(idx_u:idx_v)=0.0_wp
                IF (east_state(idx_h)   <= dry_thickness_tolerance) east_state(idx_u:idx_v)=0.0_wp
                IF (south_state(idx_h)  <= dry_thickness_tolerance) south_state(idx_u:idx_v)=0.0_wp
                IF (north_state(idx_h)  <= dry_thickness_tolerance) north_state(idx_u:idx_v)=0.0_wp

                x_stencil = [x_comp(j-1),x_comp(j),x_comp(j+1)]
                y_stencil = [y_comp(k-1),y_comp(k),y_comp(k+1)]

                IF (uniform_generalized_minmod) THEN
                    ! Specialized interior kernel for limiter=3 on all components.
                    ! These array expressions expose the five-component primitive
                    ! direction to the compiler without a per-component SELECT CASE.
                    ! Denominators and nested minmod follow geometry_2d::limit.
                    ! Using the original coordinate differences also preserves
                    ! the established nonuniform-grid arithmetic where applicable.
                    qrec_prime_x=minmod(                                                                                 &
                        (east_state-west_state)/(x_stencil(3)-x_stencil(1)),                                          &
                        theta*minmod(                                                                                  &
                        (east_state-center_state)/(x_stencil(3)-x_stencil(2)),                                        &
                        (center_state-west_state)/(x_stencil(2)-x_stencil(1))))
                    qrec_prime_y=minmod(                                                                                 &
                        (north_state-south_state)/(y_stencil(3)-y_stencil(1)),                                        &
                        theta*minmod(                                                                                  &
                        (north_state-center_state)/(y_stencil(3)-y_stencil(2)),                                       &
                        (center_state-south_state)/(y_stencil(2)-y_stencil(1))))

                    ! Opposite faces share a limited slope; this matches the
                    ! previous symmetric reconstruction for every component.
                    qrecW=center_state-reconstr_coeff*dx2*qrec_prime_x
                    qrecE=center_state+reconstr_coeff*dx2*qrec_prime_x
                    qrecS=center_state-reconstr_coeff*dy2*qrec_prime_y
                    qrecN=center_state+reconstr_coeff*dy2*qrec_prime_y

                ELSE
                    ! General fallback preserves every configured limiter,
                    ! including mixtures of methods across primitive variables.
                    DO CONCURRENT (i=1:n_vars)
                        CALL limit([west_state(i),center_state(i),east_state(i)],                                      &
                            x_stencil,limiter(i),qrec_prime_x(i))
                        CALL limit([south_state(i),center_state(i),north_state(i)],                                    &
                            y_stencil,limiter(i),qrec_prime_y(i))
                        qrecW(i)=center_state(i)-reconstr_coeff*dx2*qrec_prime_x(i)
                        qrecE(i)=center_state(i)+reconstr_coeff*dx2*qrec_prime_x(i)
                        qrecS(i)=center_state(i)-reconstr_coeff*dy2*qrec_prime_y(i)
                        qrecN(i)=center_state(i)+reconstr_coeff*dy2*qrec_prime_y(i)
                    END DO
                END IF
            ELSE

            vars_loop:DO i=1,n_vars

                qrec_stencil(2) = reconstruction_cell_value(qp_expl,i,j,k)
                bcW_effective=bcW(i)
                bcE_effective=bcE(i)
                bcS_effective=bcS(i)
                bcN_effective=bcN(i)
                ! The input BCs on components 2/3 prescribe hu/hv, not velocities.
                ! Convert their trace values to velocities for the primitive limiter.
                IF (i==idx_u .OR. i==idx_v) THEN
                    IF (j==1 .AND. comp_cells_x>1) CALL discharge_bc_to_velocity(                                     &
                        bcW(i),qp_expl(1,j,k),qp_expl(i,j,k),                                                          &
                        qp_expl(1,j+1,k),qp_expl(i,j+1,k),qrecW(1),dx2,one_by_dx,-1,                                  &
                        i==idx_u,bcW_effective)
                    IF (j==comp_cells_x .AND. comp_cells_x>1) CALL discharge_bc_to_velocity(                          &
                        bcE(i),qp_expl(1,j,k),qp_expl(i,j,k),                                                          &
                        qp_expl(1,j-1,k),qp_expl(i,j-1,k),qrecE(1),dx2,one_by_dx,1,                                   &
                        i==idx_u,bcE_effective)
                    IF (k==1 .AND. comp_cells_y>1) CALL discharge_bc_to_velocity(                                     &
                        bcS(i),qp_expl(1,j,k),qp_expl(i,j,k),                                                          &
                        qp_expl(1,j,k+1),qp_expl(i,j,k+1),qrecS(1),dy2,one_by_dy,-1,                                  &
                        i==idx_v,bcS_effective)
                    IF (k==comp_cells_y .AND. comp_cells_y>1) CALL discharge_bc_to_velocity(                          &
                        bcN(i),qp_expl(1,j,k),qp_expl(i,j,k),                                                          &
                        qp_expl(1,j,k-1),qp_expl(i,j,k-1),qrecN(1),dy2,one_by_dy,1,                                   &
                        i==idx_v,bcN_effective)
                END IF

                ! x direction
                check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

                    ! west boundary
                    check_x_boundary:IF ( j .EQ. 1 ) THEN

                        x_stencil(1) = x_stag(1)
                        x_stencil(3) = x_comp(j+1)

                        IF ( source_cell(j,k).EQ.2 ) THEN

                            CALL boundary_reconstruction_slope(                                                        &
                                qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j+1,k),                           &
                                x_stag(1), x_comp(j), x_comp(j+1), 0, source_bdry(i),                                  &
                                one_by_dx,                                                                             &
                                limiter(i), .TRUE., qrec_prime_x(i))

                        ELSE

                            CALL boundary_reconstruction_slope(                                                        &
                                qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j+1,k),                           &
                                x_stag(1), x_comp(j), x_comp(j+1), bcW_effective%flag,                                        &
                                bcW_effective%value, one_by_dx, limiter(i), .TRUE.,                                           &
                                qrec_prime_x(i))

                        END IF

                        !east boundary
                    ELSEIF ( j .EQ. comp_cells_x ) THEN

                        x_stencil(3) = x_stag(comp_interfaces_x)
                        x_stencil(1) = x_comp(j-1)

                        IF ( source_cell(j,k).EQ.2 ) THEN

                            CALL boundary_reconstruction_slope(                                                        &
                                qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j-1,k),                           &
                                x_stag(comp_interfaces_x), x_comp(j), x_comp(j-1), 0,                                  &
                                source_bdry(i), one_by_dx, limiter(i), .FALSE.,                                        &
                                qrec_prime_x(i))

                        ELSE

                            CALL boundary_reconstruction_slope(                                                        &
                                qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j-1,k),                           &
                                x_stag(comp_interfaces_x), x_comp(j), x_comp(j-1),                                     &
                                bcE_effective%flag, bcE_effective%value, one_by_dx, limiter(i), .FALSE.,                           &
                                qrec_prime_x(i))

                        END IF

                    ELSE

                        ! internal x cells

                        x_stencil(1) = x_comp(j-1)
                        x_stencil(3) = x_comp(j+1)

                        qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j-1,k)
                        qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j+1,k)

                        ! correction for radial source inlet x-interfaces values
                        ! used for the linear reconstruction
                        IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                            IF ( sourceE(j,k) ) THEN

                                x_stencil(3) = x_stag(j+1)
                                qrec_stencil(3) = source_bdry(i)

                            ELSEIF ( sourceW(j,k) ) THEN

                                x_stencil(1) = x_stag(j)
                                qrec_stencil(1) = source_bdry(i)

                            END IF

                        END IF

                        CALL limit( qrec_stencil , x_stencil , limiter(i) , qrec_prime_x(i) )

                    ENDIF check_x_boundary

                    dq = reconstr_coeff* dx2 * qrec_prime_x(i)

                    qrecW(i) = qrec_stencil(2) - dq
                    qrecE(i) = qrec_stencil(2) + dq


                    IF ( j .EQ. 1 ) THEN

                        IF ( source_cell(j,k).EQ.2 ) THEN

                            qrecW(i) = source_bdry(i)

                        ELSE

                            ! Dirichelet boundary condition at the west of the domain
                            IF ( bcW_effective%flag .EQ. 0 ) THEN

                                qrecW(i) = bcW_effective%value

                            ELSE

                                IF ( i .EQ. 2 ) qrecW(i) = MIN( qrecW(i) , 0.0_wp )

                            END IF

                        END IF

                    END IF

                    IF ( j .EQ. comp_cells_x ) THEN

                        IF ( source_cell(j,k).EQ.2 ) THEN

                            qrecE(i) = source_bdry(i)

                        ELSE

                            ! Dirichelet boundary condition at the east of the domain
                            IF ( bcE_effective%flag .EQ. 0 ) THEN

                                qrecE(i) = bcE_effective%value

                            ELSE

                                IF ( i .EQ. 2 ) qrecE(i) = MAX( qrecE(i) , 0.0_wp )

                            END IF

                        END IF

                    END IF

                END IF check_comp_cells_x

                ! y-direction
                check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

                    ! South boundary
                    check_y_boundary:IF ( k .EQ. 1 ) THEN

                        y_stencil(1) = y_stag(1)
                        y_stencil(3) = y_comp(k+1)

                        CALL boundary_reconstruction_slope(                                                            &
                            qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j,k+1),                               &
                            y_stag(1), y_comp(k), y_comp(k+1), bcS_effective%flag,                                            &
                            bcS_effective%value,                                                                              &
                            one_by_dy, limiter(i), .TRUE., qrec_prime_y(i))

                        ! North boundary
                    ELSEIF ( k .EQ. comp_cells_y ) THEN

                            y_stencil(1) = y_comp(k-1)
                            y_stencil(3) = y_stag(comp_interfaces_y)

                            CALL boundary_reconstruction_slope(                                                        &
                            qrec_stencil(2), reconstruction_cell_value(qp_expl,i,j,k-1),                               &
                            y_stag(comp_interfaces_y), y_comp(k), y_comp(k-1),                                         &
                            bcN_effective%flag, bcN_effective%value, one_by_dy, limiter(i), .FALSE.,                               &
                            qrec_prime_y(i))

                    ELSE

                        ! Internal y cells

                        y_stencil(1) = y_comp(k-1)
                        y_stencil(3) = y_comp(k+1)

                        qrec_stencil(1) = reconstruction_cell_value(qp_expl,i,j,k-1)
                        qrec_stencil(3) = reconstruction_cell_value(qp_expl,i,j,k+1)

                        ! correction for radial source inlet y-interfaces
                        ! used for the linear reconstruction
                        IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                            IF ( sourceS(j,k) ) THEN

                                y_stencil(1) = y_stag(k)
                                qrec_stencil(1) = source_bdry(i)

                            ELSEIF ( sourceN(j,k) ) THEN

                                y_stencil(3) = y_stag(k+1)
                                qrec_stencil(3) = source_bdry(i)

                            END IF

                        END IF

                        CALL limit( qrec_stencil , y_stencil , limiter(i) , qrec_prime_y(i) )

                    ENDIF check_y_boundary

                    dq = reconstr_coeff * dy2 * qrec_prime_y(i)

                    qrecS(i) = qrec_stencil(2) - dq
                    qrecN(i) = qrec_stencil(2) + dq

                    IF ( k .EQ. 1 ) THEN

                        IF ( source_cell(j,k) .EQ. 2 ) THEN

                            qrecS(i) = source_bdry(i)

                        ELSE

                            ! Dirichlet boundary condition at the south of the domain
                            IF ( bcS_effective%flag .EQ. 0 ) THEN
                                qrecS(i) = bcS_effective%value
                            ELSE
                                IF ( i .EQ. 3 ) qrecS(i) = MIN( qrecS(i) , 0.0_wp )

                            END IF

                        END IF

                    END IF

                    IF ( k .EQ. comp_cells_y ) THEN

                        IF ( source_cell(j,k) .EQ. 2 ) THEN

                            qrecN(i) = source_bdry(i)

                        ELSE

                            ! Dirichlet boundary condition at the north of the domain
                            IF ( bcN_effective%flag .EQ. 0 ) THEN
                                qrecN(i) = bcN_effective%value
                            ELSE
                                IF ( i .EQ. 3 ) qrecN(i) = MAX( qrecN(i) , 0.0_wp )

                            END IF

                        END IF

                    END IF

                ENDIF check_comp_cells_y

            ENDDO vars_loop


            ! Radial inlet states are prescribed on the source-facing traces.
            IF (radial_source_flag .AND. source_cell(j,k)==2) THEN
                IF (sourceE(j,k)) qrecE=source_bdry
                IF (sourceW(j,k)) qrecW=source_bdry
                IF (sourceS(j,k)) qrecS=source_bdry
                IF (sourceN(j,k)) qrecN=source_bdry
            END IF


            END IF

            END IF ! active cell

            this%qp_cellW(:,j,k) = qrecW
            this%qp_cellE(:,j,k) = qrecE
            this%qp_cellS(:,j,k) = qrecS
            this%qp_cellN(:,j,k) = qrecN
            CALL reconstruct_kt_eta_cell(this,qp_expl,j,k)




            ! The transport interface arrays retain the finalized primitive traces.
            this%qp_adv_interfaceR(:,j,k) = this%qp_cellW(:,j,k)
            this%qp_adv_interfaceL(:,j+1,k) = this%qp_cellE(:,j,k)
            this%qp_adv_interfaceT(:,j,k) = this%qp_cellS(:,j,k)
            this%qp_adv_interfaceB(:,j,k+1) = this%qp_cellN(:,j,k)


        END DO

        !$OMP END PARALLEL DO

        CALL apply_hp_reconstruction( this, qp_expl, solve_cells, j_cent, k_cent, use_kt_direct )

        RETURN

    END SUBROUTINE reconstruction

    !> \brief Build a deduplicated eta workset containing active cells and immediate neighbours.
    !>
    !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
    !> \param[in] solve_cells Number of entries in the active-cell index lists.
    !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
    !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
    !>
    !> \note Writes only the parent workspace eta mask/list; immediate neighbours provide the
    !>       additional continuity stencil needed by active cells.

    SUBROUTINE build_hp_workset(this, solve_cells, j_cent, k_cent)

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: solve_cells
        INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
        INTEGER :: l, j, k

        DO l = 1, this%hp_eta_cells
                this%hp_eta_mask(this%hp_eta_j(l),this%hp_eta_k(l)) = .FALSE.

        END DO

        this%hp_eta_cells = 0

        ! Active cells first, followed by unique stencil halo neighbours.
        DO l = 1, solve_cells
                j = j_cent(l)
                k = k_cent(l)
                CALL append_eta_cell(j,k)

        END DO

        DO l = 1, solve_cells
                j = j_cent(l)
                k = k_cent(l)
                IF ( j .GT. 1 ) CALL append_eta_cell(j-1,k)
                IF ( j .LT. comp_cells_x ) CALL append_eta_cell(j+1,k)
                IF ( k .GT. 1 ) CALL append_eta_cell(j,k-1)
                IF ( k .LT. comp_cells_y ) CALL append_eta_cell(j,k+1)

        END DO

    CONTAINS

        !> \brief Append one cell to the parent eta workset only if it is not already present.
        !>
        !> \param[in] j_cell X index of the cell appended to the parent eta workset.
        !> \param[in] k_cell Y index of the cell appended to the parent eta workset.
        !>
        !> \note Updates the host-associated this workspace; the active-cell workset is not
        !> modified.

    SUBROUTINE append_eta_cell(j_cell,k_cell)

            INTEGER, INTENT(IN) :: j_cell, k_cell

            IF ( this%hp_eta_mask(j_cell,k_cell) ) RETURN

            this%hp_eta_mask(j_cell,k_cell) = .TRUE.
            this%hp_eta_cells = this%hp_eta_cells + 1
            this%hp_eta_j(this%hp_eta_cells) = j_cell
            this%hp_eta_k(this%hp_eta_cells) = k_cell

    END SUBROUTINE append_eta_cell

    END SUBROUTINE build_hp_workset

    !******************************************************************************
    !> \brief Reconstruct the KP free surface for KT dissipation and path terms.
    !> \param[in,out] this Workspace containing the active HP-eta cell list and traces.
    !> \param[in] qp_expl Cell-centered physical states at the explicit stage.
    !******************************************************************************
    !> @brief Reconstruct free-surface elevation and transport-depth candidates for one cell.
    !> @details Reconstructs local slopes using configured limiter and topography;
    !>          preserves KT and Froude transport-depth candidates for face selection.
    !> @param[in,out] this Reconstruction workspace containing reconstructed traces.
    !> @param[in] qp_expl Cell-centered primitive states including stencil halo.
    !> @param[in] j X index of the workset cell.
    !> @param[in] k Y index of the workset cell.
    SUBROUTINE reconstruct_kt_eta_cell(this,qp_expl,j,k)
        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: qp_expl(:,:,:)
        INTEGER, INTENT(IN) :: j,k
        REAL(wp) :: sx(3),sy(3),hx(3),hy(3),xx(3),yy(3)
        REAL(wp) :: slope_x,slope_y,slope_hx,slope_hy,eta_center
        REAL(wp) :: eta_w,eta_e,eta_s,eta_n,slope_eta2,slope_minmod
        REAL(wp) :: eta_w_mm,eta_e_mm,eta_s_mm,eta_n_mm
        REAL(wp) :: slope_eta2_x,slope_eta2_y
        REAL(wp) :: slope_blend_x,slope_blend_y,eta_blend_w,eta_blend_e
        REAL(wp) :: eta_blend_s,eta_blend_n,fr_weight_x,fr_weight_y
            eta_center=qp_expl(1,j,k)+B_cent(j,k)
            slope_x=0.0_wp
            slope_y=0.0_wp
            slope_hx=0.0_wp
            slope_hy=0.0_wp

            IF ((j.GT.1).AND.(j.LT.comp_cells_x)) THEN
                sx=[qp_expl(1,j-1,k)+B_cent(j-1,k),eta_center,                                                         &
                    qp_expl(1,j+1,k)+B_cent(j+1,k)]
                hx=[qp_expl(1,j-1,k),qp_expl(1,j,k),qp_expl(1,j+1,k)]
                xx=[x_comp(j-1),x_comp(j),x_comp(j+1)]
                CALL limit(sx,xx,limiter(1),slope_x)
                ! h_x is needed only by minmod and Froude-blended transport
                ! depths; the eta/path trace above is required by every mode.
                IF (kt_trace_reconstruction=='minmod' .OR.                                                           &
                    kt_trace_reconstruction=='froude_blend') CALL limit(hx,xx,limiter(1),slope_hx)

            END IF

            IF ((k.GT.1).AND.(k.LT.comp_cells_y)) THEN
                sy=[qp_expl(1,j,k-1)+B_cent(j,k-1),eta_center,                                                         &
                    qp_expl(1,j,k+1)+B_cent(j,k+1)]
                hy=[qp_expl(1,j,k-1),qp_expl(1,j,k),qp_expl(1,j,k+1)]
                yy=[y_comp(k-1),y_comp(k),y_comp(k+1)]
                CALL limit(sy,yy,limiter(1),slope_y)
                IF (kt_trace_reconstruction=='minmod' .OR.                                                           &
                    kt_trace_reconstruction=='froude_blend') CALL limit(hy,yy,limiter(1),slope_hy)

            END IF

            eta_w=eta_center-reconstr_coeff*dx2*slope_x
            eta_e=eta_center+reconstr_coeff*dx2*slope_x
            eta_s=eta_center-reconstr_coeff*dy2*slope_y
            eta_n=eta_center+reconstr_coeff*dy2*slope_y

            ! K&P positivity correction for the eta_x traces used by the
            ! hydrostatic/path discretization.
            IF ( eta_e .LT. B_face_x(j+1,k) ) THEN
                eta_e=B_face_x(j+1,k)
                eta_w=2.0_wp*eta_center-eta_e
            ELSEIF ( eta_w .LT. B_face_x(j,k) ) THEN
                eta_w=B_face_x(j,k)
                eta_e=2.0_wp*eta_center-eta_w

            END IF

            IF ( eta_n .LT. B_face_y(j,k+1) ) THEN
                eta_n=B_face_y(j,k+1)
                eta_s=2.0_wp*eta_center-eta_n
            ELSEIF ( eta_s .LT. B_face_y(j,k) ) THEN
                eta_s=B_face_y(j,k)
                eta_n=2.0_wp*eta_center-eta_s

            END IF

            this%eta_cellW(j,k)=eta_w
            this%eta_cellE(j,k)=eta_e
            this%eta_cellS(j,k)=eta_s
            this%eta_cellN(j,k)=eta_n

            IF (kt_trace_reconstruction .EQ. 'minmod') THEN
                ! Alternative transport traces: minmod(eta_x, h_x+B_x), followed
                ! by the face-bed positivity correction. Hydrostatic eta stays above.
                slope_eta2=reconstr_coeff*slope_hx +                                                                   &
                    (B_face_x(j+1,k)-B_face_x(j,k))*one_by_dx
                slope_minmod=minmod_slope_pair(reconstr_coeff*slope_x,slope_eta2)
                eta_w_mm=eta_center-dx2*slope_minmod
                eta_e_mm=eta_center+dx2*slope_minmod

                IF ( eta_e_mm .LT. B_face_x(j+1,k) ) THEN
                    eta_e_mm=B_face_x(j+1,k)
                    eta_w_mm=2.0_wp*eta_center-eta_e_mm
                ELSEIF ( eta_w_mm .LT. B_face_x(j,k) ) THEN
                    eta_w_mm=B_face_x(j,k)
                    eta_e_mm=2.0_wp*eta_center-eta_w_mm

                END IF

                this%h_minmod_interfaceR(j,k)=                                                                         &
                    MAX(eta_w_mm-B_face_x(j,k),0.0_wp)
                this%h_minmod_interfaceL(j+1,k)=                                                                       &
                    MAX(eta_e_mm-B_face_x(j+1,k),0.0_wp)

                slope_eta2=reconstr_coeff*slope_hy +                                                                   &
                    (B_face_y(j,k+1)-B_face_y(j,k))*one_by_dy
                slope_minmod=minmod_slope_pair(reconstr_coeff*slope_y,slope_eta2)
                eta_s_mm=eta_center-dy2*slope_minmod
                eta_n_mm=eta_center+dy2*slope_minmod

                IF ( eta_n_mm .LT. B_face_y(j,k+1) ) THEN
                    eta_n_mm=B_face_y(j,k+1)
                    eta_s_mm=2.0_wp*eta_center-eta_n_mm
                ELSEIF ( eta_s_mm .LT. B_face_y(j,k) ) THEN
                    eta_s_mm=B_face_y(j,k)
                    eta_n_mm=2.0_wp*eta_center-eta_s_mm

                END IF

                this%h_minmod_interfaceT(j,k)=                                                                         &
                    MAX(eta_s_mm-B_face_y(j,k),0.0_wp)
                this%h_minmod_interfaceB(j,k+1)=                                                                       &
                    MAX(eta_n_mm-B_face_y(j,k+1),0.0_wp)

            END IF

IF (kt_trace_reconstruction .EQ. 'froude_blend') THEN
            fr_weight_x=froude_h_blend_weight(qp_expl(1,j,k),                                                          &
                qp_expl(idx_u,j,k),grav)
            fr_weight_y=froude_h_blend_weight(qp_expl(1,j,k),                                                          &
                qp_expl(idx_v,j,k),grav)
            slope_eta2_x=reconstr_coeff*slope_hx +                                                                     &
                (B_face_x(j+1,k)-B_face_x(j,k))*one_by_dx
            slope_eta2_y=reconstr_coeff*slope_hy +                                                                     &
                (B_face_y(j,k+1)-B_face_y(j,k))*one_by_dy
            slope_blend_x=(1.0_wp-fr_weight_x)*reconstr_coeff*slope_x +                                                &
                fr_weight_x*slope_eta2_x
            slope_blend_y=(1.0_wp-fr_weight_y)*reconstr_coeff*slope_y +                                                &
                fr_weight_y*slope_eta2_y

            eta_blend_w=eta_center-dx2*slope_blend_x
            eta_blend_e=eta_center+dx2*slope_blend_x

            IF ( eta_blend_e .LT. B_face_x(j+1,k) ) THEN
                eta_blend_e=B_face_x(j+1,k)
                eta_blend_w=2.0_wp*eta_center-eta_blend_e
            ELSEIF ( eta_blend_w .LT. B_face_x(j,k) ) THEN
                eta_blend_w=B_face_x(j,k)
                eta_blend_e=2.0_wp*eta_center-eta_blend_w

            END IF

            this%h_froude_interfaceR(j,k)=                                                                             &
                MAX(eta_blend_w-B_face_x(j,k),0.0_wp)
            this%h_froude_interfaceL(j+1,k)=                                                                           &
                MAX(eta_blend_e-B_face_x(j+1,k),0.0_wp)

            eta_blend_s=eta_center-dy2*slope_blend_y
            eta_blend_n=eta_center+dy2*slope_blend_y

            IF ( eta_blend_n .LT. B_face_y(j,k+1) ) THEN
                eta_blend_n=B_face_y(j,k+1)
                eta_blend_s=2.0_wp*eta_center-eta_blend_n
            ELSEIF ( eta_blend_s .LT. B_face_y(j,k) ) THEN
                eta_blend_s=B_face_y(j,k)
                eta_blend_n=2.0_wp*eta_center-eta_blend_s

            END IF

            this%h_froude_interfaceT(j,k)=                                                                             &
                MAX(eta_blend_s-B_face_y(j,k),0.0_wp)
            this%h_froude_interfaceB(j,k+1)=                                                                           &
                MAX(eta_blend_n-B_face_y(j,k+1),0.0_wp)

            END IF

    END SUBROUTINE reconstruct_kt_eta_cell

    !******************************************************************************
    !> \brief Compute a smooth blend weight from the local Froude number.
    !> \param[in] h Flow depth.
    !> \param[in] u_normal Normal velocity.
    !> \param[in] gravity Effective gravitational acceleration.
    !> \return Blend weight between zero and one.
    !******************************************************************************
    PURE FUNCTION froude_h_blend_weight(h,u_normal,gravity) RESULT(weight)
        REAL(wp), INTENT(IN) :: h,u_normal,gravity
        REAL(wp) :: weight,fr,scaled

        weight=0.0_wp
        IF ( (h .LE. dry_thickness_tolerance) .OR. (gravity .LE. 0.0_wp) ) RETURN

        fr=ABS(u_normal)/SQRT(gravity*h)
        scaled=MIN(1.0_wp,MAX(0.0_wp,                                                                                  &
            (fr-fr_blend_lower)/(fr_blend_upper-fr_blend_lower)))
        weight=scaled*scaled*(3.0_wp-2.0_wp*scaled)
    END FUNCTION froude_h_blend_weight

    !******************************************************************************
    !> \brief Apply the two-argument minmod limiter to a pair of slopes.
    !> \param[in] slope_a First slope candidate.
    !> \param[in] slope_b Second slope candidate.
    !> \return Same-sign slope with the smaller magnitude, or zero.
    !******************************************************************************
    PURE FUNCTION minmod_slope_pair(slope_a,slope_b) RESULT(slope)
        REAL(wp), INTENT(IN) :: slope_a,slope_b
        REAL(wp) :: slope
        slope=0.0_wp
        IF ( slope_a*slope_b .GT. 0.0_wp )                                                                             &
            slope=SIGN(MIN(ABS(slope_a),ABS(slope_b)),slope_a)
    END FUNCTION minmod_slope_pair

    !> \brief Apply one side's boundary condition to the linear reconstruction slope.
    !> \details Flag 0 limits a Dirichlet boundary stencil, flag 1 prescribes the
    !>          Neumann slope, and flag 2 uses the adjacent-cell secant.
    !> \param[in] center_value Value in the boundary-adjacent cell.
    !> \param[in] adjacent_value Value in the next cell inward from the boundary.
    !> \param[in] boundary_coordinate Coordinate of the physical boundary.
    !> \param[in] center_coordinate Coordinate of the boundary-adjacent cell center.
    !> \param[in] adjacent_coordinate Coordinate of the next cell inward.
    !> \param[in] boundary_flag Boundary mode (Dirichlet, Neumann or secant).
    !> \param[in] boundary_value Dirichlet value or Neumann derivative.
    !> \param[in] inverse_spacing Inverse uniform grid spacing for secant mode.
    !> \param[in] limiter_id Limiter selector.
    !> \param[in] lower_boundary True for west/south, false for east/north.
    !> \param[out] slope Limited or prescribed derivative.
    SUBROUTINE boundary_reconstruction_slope(center_value,adjacent_value,                                              &
            boundary_coordinate,center_coordinate,adjacent_coordinate,                                                 &
            boundary_flag,boundary_value,inverse_spacing,limiter_id,                                                   &
            lower_boundary,slope)
        REAL(wp), INTENT(IN) :: center_value,adjacent_value
        REAL(wp), INTENT(IN) :: boundary_coordinate,center_coordinate
        REAL(wp), INTENT(IN) :: adjacent_coordinate
        REAL(wp), INTENT(IN) :: boundary_value,inverse_spacing
        INTEGER, INTENT(IN) :: boundary_flag,limiter_id
        LOGICAL, INTENT(IN) :: lower_boundary
        REAL(wp), INTENT(OUT) :: slope
        REAL(wp) :: stencil_values(3),stencil_coordinates(3)

        IF (boundary_flag .EQ. 0) THEN
            stencil_values(2)=center_value
            stencil_coordinates(2)=center_coordinate

            IF (lower_boundary) THEN
                stencil_values(1)=boundary_value
                stencil_values(3)=adjacent_value
                stencil_coordinates(1)=boundary_coordinate
                stencil_coordinates(3)=adjacent_coordinate
            ELSE
                stencil_values(1)=adjacent_value
                stencil_values(3)=boundary_value
                stencil_coordinates(1)=adjacent_coordinate
                stencil_coordinates(3)=boundary_coordinate

            END IF

            CALL limit(stencil_values,stencil_coordinates,limiter_id,slope)
        ELSEIF (boundary_flag .EQ. 1) THEN
            slope=boundary_value
        ELSEIF (boundary_flag .EQ. 2) THEN

            IF (lower_boundary) THEN
                slope=(adjacent_value-center_value)*inverse_spacing
            ELSE
                slope=(center_value-adjacent_value)*inverse_spacing

            END IF

        END IF
    END SUBROUTINE boundary_reconstruction_slope

    !> \brief Read one centre value for a direct reconstruction stencil.
    !> \details KT uses the cell-centred velocities already obtained by
    !>          qc_to_qp, which applies the v2 mass-based desingularization.
    !>          Auxiliary u/v are set to zero below the HP dry-depth threshold
    !>          so a dry cell cannot supply a velocity to a neighbouring limiter.
    !>          Depth, volumetric momenta, temperature and composition remain
    !>          unchanged; neither the input physical nor conservative state is
    !>          modified. The function is allocation-free and thread-safe.
    !> \param[in] qp_center Read-only physical centre fields, including halo.
    !> \param[in] component Primitive component index; idx_u/idx_v are velocities.
    !> \param[in] j X index of the centre stencil entry.
    !> \param[in] k Y index of the centre stencil entry.
    !> \return Original centre value, or zero for a dry velocity.
    PURE FUNCTION reconstruction_cell_value(qp_center,component,j,k) RESULT(value)
        REAL(wp), INTENT(IN) :: qp_center(:,:,:)
        INTEGER, INTENT(IN) :: component,j,k
        REAL(wp) :: value

        value = qp_center(component,j,k)

        IF ((component .EQ. idx_u) .OR. (component .EQ. idx_v)) THEN
            IF (qp_center(1,j,k) .LE. dry_thickness_tolerance) value = 0.0_wp

        END IF
    END FUNCTION reconstruction_cell_value

    !******************************************************************************
    !> \brief Evaluate hydrostatic candidates, blend them and commit final cell/face traces.
    !>
    !> The first reconstruction pass supplies direct h_x-based candidates for
    !> advective fluxes. The HP phases build eta-based traces for dissipative
    !> jumps and hydrostatic paths, then apply the existing bounded momentum
    !> correction on the active target cells.
    !>
    !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.
    !> \param[in] qp_center Unmodified cell-centered physical states used by all HP candidate
    !> phases.
    !> \param[in] solve_cells Number of entries in the active-cell index lists.
    !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
    !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
    !> \param[in] use_kt_direct Use direct KT traces instead of the HP dry-state projection.
    !******************************************************************************

    SUBROUTINE apply_hp_reconstruction( this, qp_center, solve_cells, j_cent, k_cent, use_kt_direct )

    USE geometry_2d, ONLY : sourceW, sourceE, sourceS, sourceN

    IMPLICIT NONE

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: qp_center(:,:,:)
        INTEGER, INTENT(IN) :: solve_cells
        INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
        LOGICAL, INTENT(IN) :: use_kt_direct

        INTEGER :: j, k, l

        ! Final primitive traces are already populated by reconstruct().
        ! This region handles only boundaries that require completed cell traces.
        !$OMP PARALLEL PRIVATE(l,j,k)

        ! Preserve the solver's existing external ghost convention after replacing
        ! the cell-owned traces.
        !$OMP DO
        DO k = 1, comp_cells_y
            this%qp_adv_interfaceL(:,1,k) = this%qp_adv_interfaceR(:,1,k)
            this%qp_adv_interfaceR(:,comp_interfaces_x,k) =                                                            &
                this%qp_adv_interfaceL(:,comp_interfaces_x,k)

        END DO
        !$OMP END DO

        !$OMP DO
        DO j = 1, comp_cells_x
            this%qp_adv_interfaceB(:,j,1) = this%qp_adv_interfaceT(:,j,1)
            this%qp_adv_interfaceT(:,j,comp_interfaces_y) =                                                            &
                this%qp_adv_interfaceB(:,j,comp_interfaces_y)

        END DO
        !$OMP END DO

        ! Rebuild the internal reflecting side of radial-source cells from the
        ! reconstructed transport state. Reflect primitive momentum and normal velocity.
        IF ( radial_source_flag ) THEN
            !$OMP DO
            DO l = 1, solve_cells
                j = j_cent(l)
                k = k_cent(l)
                IF ( source_cell(j,k) .NE. 2 ) CYCLE

                IF ( sourceE(j,k) ) THEN
                    this%qp_adv_interfaceR(:,j+1,k) =                                                                  &
                        this%qp_adv_interfaceL(:,j+1,k)
                    this%qp_adv_interfaceR(2,j+1,k) =                                                                  &
                        -this%qp_adv_interfaceR(2,j+1,k)
                    this%eta_interfaceR(j+1,k) = this%eta_interfaceL(j+1,k)
                ELSEIF ( sourceW(j,k) ) THEN
                    this%qp_adv_interfaceL(:,j,k) = this%qp_adv_interfaceR(:,j,k)
                    this%qp_adv_interfaceL(2,j,k) = -this%qp_adv_interfaceL(2,j,k)
                    this%eta_interfaceL(j,k) = this%eta_interfaceR(j,k)

                END IF

                IF ( sourceN(j,k) ) THEN
                    this%qp_adv_interfaceT(:,j,k+1) =                                                                  &
                        this%qp_adv_interfaceB(:,j,k+1)
                    this%qp_adv_interfaceT(3,j,k+1) =                                                                  &
                        -this%qp_adv_interfaceT(3,j,k+1)
                    this%eta_interfaceT(j,k+1) = this%eta_interfaceB(j,k+1)
                ELSEIF ( sourceS(j,k) ) THEN
                    this%qp_adv_interfaceB(:,j,k) = this%qp_adv_interfaceT(:,j,k)
                    this%qp_adv_interfaceB(3,j,k) = -this%qp_adv_interfaceB(3,j,k)
                    this%eta_interfaceB(j,k) = this%eta_interfaceT(j,k)

                END IF

            END DO
            !$OMP END DO
        END IF

        !$OMP END PARALLEL

    END SUBROUTINE apply_hp_reconstruction

    !> \brief Release all reconstructed-state and HP work arrays.
    !>
    !> \param[in,out] this Persistent reconstructed traces, HP candidates and per-thread scratch.

    SUBROUTINE finalize_reconstruction( this )

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this



        DEALLOCATE( this%qp_adv_interfaceL )
        DEALLOCATE( this%h_selected_interfaceL )
        DEALLOCATE( this%qp_adv_interfaceR )
        DEALLOCATE( this%h_selected_interfaceR )
        DEALLOCATE( this%qp_adv_interfaceB )
        DEALLOCATE( this%h_selected_interfaceB )
        DEALLOCATE( this%qp_adv_interfaceT )
        DEALLOCATE( this%h_selected_interfaceT )

        DEALLOCATE( this%qp_cellW )
        DEALLOCATE( this%qp_cellE )
        DEALLOCATE( this%qp_cellS )
        DEALLOCATE( this%qp_cellN )

        DEALLOCATE( this%eta_cellW )
        DEALLOCATE( this%eta_cellE )
        DEALLOCATE( this%eta_cellS )
        DEALLOCATE( this%eta_cellN )
        DEALLOCATE( this%eta_interfaceL )
        DEALLOCATE( this%eta_interfaceR )
        DEALLOCATE( this%eta_interfaceB )
        DEALLOCATE( this%eta_interfaceT )
        DEALLOCATE( this%h_minmod_interfaceL )
        DEALLOCATE( this%h_minmod_interfaceR )
        DEALLOCATE( this%h_minmod_interfaceB )
        DEALLOCATE( this%h_minmod_interfaceT )
        DEALLOCATE( this%h_froude_interfaceL )
        DEALLOCATE( this%h_froude_interfaceR )
        DEALLOCATE( this%h_froude_interfaceB )
        DEALLOCATE( this%h_froude_interfaceT )

        DEALLOCATE( this%hp_eta_mask )
        DEALLOCATE( this%hp_eta_j )
        DEALLOCATE( this%hp_eta_k )


    END SUBROUTINE finalize_reconstruction

    !> @brief Store the HP workset built once by the active-domain refresh.
    !> @param[in,out] this Reconstruction workspace whose HP list is updated.
    !> @param[in] hp_cells Number of HP workset entries.
    !> @param[in] hp_j X-cell indices of the HP workset.
    !> @param[in] hp_k Y-cell indices of the HP workset.
    SUBROUTINE prepare_hp_workset(this, hp_cells, hp_j, hp_k)
        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: hp_cells
        INTEGER, INTENT(IN) :: hp_j(:), hp_k(:)
        this%hp_eta_cells = hp_cells
        this%hp_eta_j(1:hp_cells) = hp_j(1:hp_cells)
        this%hp_eta_k(1:hp_cells) = hp_k(1:hp_cells)
        this%hp_workset_prepared = .TRUE.
    END SUBROUTINE prepare_hp_workset

    !> @brief Reconstruct KT free-surface and transport-depth traces on the HP workset.
    !> @details Applies the configured trace reconstruction to prepared HP cells.
    !> @param[in,out] this Reconstruction workspace and output face traces.
    !> @param[in] qp_expl Cell-centered primitive states with stencil halo.
    SUBROUTINE reconstruct_kt_eta_traces(this,qp_expl)
    USE geometry_2d, ONLY : sourceW,sourceE,sourceS,sourceN

        CLASS(reconstruction_workspace_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: qp_expl(:,:,:)
        INTEGER :: l,j,k,kt_mode

        ! Decode the configuration once, outside the reconstruction face loop.
        SELECT CASE (TRIM(kt_trace_reconstruction))
        CASE ('h_x')
            kt_mode=1
        CASE ('minmod')
            kt_mode=3
        CASE ('froude_blend')
            kt_mode=4
        CASE DEFAULT
            kt_mode=2

        END SELECT

        IF (kt_trace_reconstruction .EQ. 'froude_blend') THEN

            DO k=1,comp_cells_y
                this%h_froude_interfaceL(1,k)=this%h_froude_interfaceR(1,k)
                this%h_froude_interfaceR(comp_interfaces_x,k)=                                                         &
                    this%h_froude_interfaceL(comp_interfaces_x,k)

            END DO

            DO j=1,comp_cells_x
                this%h_froude_interfaceB(j,1)=this%h_froude_interfaceT(j,1)
                this%h_froude_interfaceT(j,comp_interfaces_y)=                                                         &
                    this%h_froude_interfaceB(j,comp_interfaces_y)

            END DO

        END IF

        !$OMP PARALLEL DO PRIVATE(l,j,k)
        DO l=1,this%hp_eta_cells
            j=this%hp_eta_j(l)
            k=this%hp_eta_k(l)
            this%eta_interfaceR(j,k)=this%eta_cellW(j,k)
            this%eta_interfaceL(j+1,k)=this%eta_cellE(j,k)
            this%eta_interfaceT(j,k)=this%eta_cellS(j,k)
            this%eta_interfaceB(j,k+1)=this%eta_cellN(j,k)

            IF (j.EQ.1) THEN
                this%eta_interfaceL(1,k)=this%eta_interfaceR(1,k)
                IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                           &
                    this%h_minmod_interfaceL(1,k)=this%h_minmod_interfaceR(1,k)

            END IF

            IF (j.EQ.comp_cells_x) THEN
                this%eta_interfaceR(comp_interfaces_x,k)=                                                              &
                    this%eta_interfaceL(comp_interfaces_x,k)
                IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                           &
                    this%h_minmod_interfaceR(comp_interfaces_x,k)=                                                     &
                    this%h_minmod_interfaceL(comp_interfaces_x,k)

            END IF

            IF (k.EQ.1) THEN
                this%eta_interfaceB(j,1)=this%eta_interfaceT(j,1)
                IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                           &
                    this%h_minmod_interfaceB(j,1)=this%h_minmod_interfaceT(j,1)

            END IF

            IF (k.EQ.comp_cells_y) THEN
                this%eta_interfaceT(j,comp_interfaces_y)=                                                              &
                    this%eta_interfaceB(j,comp_interfaces_y)
                IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                           &
                    this%h_minmod_interfaceT(j,comp_interfaces_y)=                                                     &
                    this%h_minmod_interfaceB(j,comp_interfaces_y)

            END IF

        END DO
        !$OMP END PARALLEL DO

        DO k=1,comp_cells_y
            this%qp_adv_interfaceL(:,1,k)=this%qp_adv_interfaceR(:,1,k)
            this%qp_adv_interfaceR(:,comp_interfaces_x,k)=                                                             &
                this%qp_adv_interfaceL(:,comp_interfaces_x,k)

        END DO

        DO j=1,comp_cells_x
            this%qp_adv_interfaceB(:,j,1)=this%qp_adv_interfaceT(:,j,1)
            this%qp_adv_interfaceT(:,j,comp_interfaces_y)=                                                             &
                this%qp_adv_interfaceB(:,j,comp_interfaces_y)

        END DO
        ! reconstruct() has already reflected the primitive states
        ! at the internal radial-source wall. The eta pass above overwrites its
        ! scalar ghost traces: restore equal depths on both sides, including the
        ! selected transport reconstruction, before evaluating speeds and fluxes.
        ! Injection is volumetric; no mass may cross this excluded inner boundary.
        IF ( radial_source_flag ) THEN

            DO l=1,this%hp_eta_cells
                j=this%hp_eta_j(l)
                k=this%hp_eta_k(l)
                IF ( source_cell(j,k) .NE. 2 ) CYCLE

                IF ( sourceE(j,k) ) THEN
                    this%eta_interfaceR(j+1,k)=this%eta_interfaceL(j+1,k)
                    IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                       &
                        this%h_minmod_interfaceR(j+1,k)=this%h_minmod_interfaceL(j+1,k)
                    IF ( kt_trace_reconstruction .EQ. 'froude_blend' )                                                 &
                        this%h_froude_interfaceR(j+1,k)=this%h_froude_interfaceL(j+1,k)
                ELSEIF ( sourceW(j,k) ) THEN
                    this%eta_interfaceL(j,k)=this%eta_interfaceR(j,k)
                    IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                       &
                        this%h_minmod_interfaceL(j,k)=this%h_minmod_interfaceR(j,k)
                    IF ( kt_trace_reconstruction .EQ. 'froude_blend' )                                                 &
                        this%h_froude_interfaceL(j,k)=this%h_froude_interfaceR(j,k)

                END IF

                IF ( sourceN(j,k) ) THEN
                    this%eta_interfaceT(j,k+1)=this%eta_interfaceB(j,k+1)
                    IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                       &
                        this%h_minmod_interfaceT(j,k+1)=this%h_minmod_interfaceB(j,k+1)
                    IF ( kt_trace_reconstruction .EQ. 'froude_blend' )                                                 &
                        this%h_froude_interfaceT(j,k+1)=this%h_froude_interfaceB(j,k+1)
                ELSEIF ( sourceS(j,k) ) THEN
                    this%eta_interfaceB(j,k)=this%eta_interfaceT(j,k)
                    IF ( kt_trace_reconstruction .EQ. 'minmod' )                                                       &
                        this%h_minmod_interfaceB(j,k)=this%h_minmod_interfaceT(j,k)
                    IF ( kt_trace_reconstruction .EQ. 'froude_blend' )                                                 &
                        this%h_froude_interfaceB(j,k)=this%h_froude_interfaceT(j,k)

                END IF

            END DO

        END IF

        ! Select final depths after eta and radial boundary corrections.
        ! Keep the direct h_x candidate available for characteristic speed bounds.
        !$OMP PARALLEL DO PRIVATE(l,j,k)
        DO l=1,this%hp_eta_cells
            j=this%hp_eta_j(l)
            k=this%hp_eta_k(l)
            CALL select_final_transport_trace(this%qp_adv_interfaceR(:,j,k), this%eta_interfaceR(j,k),B_face_x(j,k),   &
                this%h_minmod_interfaceR(j,k),this%h_froude_interfaceR(j,k),                                           &
                lateral_source_flag .AND. j==1 .AND. source_cell(j,k)==2,                                              &
                this%h_selected_interfaceR(j,k),'KT reconstruction west',kt_mode)
            CALL select_final_transport_trace(this%qp_adv_interfaceL(:,j+1,k),                                         &
                this%eta_interfaceL(j+1,k),B_face_x(j+1,k),                                                            &
                this%h_minmod_interfaceL(j+1,k),this%h_froude_interfaceL(j+1,k),                                       &
                lateral_source_flag .AND. j==comp_cells_x .AND. source_cell(j,k)==2,                                   &
                this%h_selected_interfaceL(j+1,k),'KT reconstruction east',kt_mode)
            CALL select_final_transport_trace(this%qp_adv_interfaceT(:,j,k), this%eta_interfaceT(j,k),B_face_y(j,k),   &
                this%h_minmod_interfaceT(j,k),this%h_froude_interfaceT(j,k),                                           &
                lateral_source_flag .AND. k==1 .AND. source_cell(j,k)==2,                                              &
                this%h_selected_interfaceT(j,k),'KT reconstruction south',kt_mode)
            CALL select_final_transport_trace(this%qp_adv_interfaceB(:,j,k+1),                                         &
                this%eta_interfaceB(j,k+1),B_face_y(j,k+1),                                                            &
                this%h_minmod_interfaceB(j,k+1),this%h_froude_interfaceB(j,k+1),                                       &
                lateral_source_flag .AND. k==comp_cells_y .AND. source_cell(j,k)==2,                                   &
                this%h_selected_interfaceB(j,k+1),'KT reconstruction north',kt_mode)

        END DO
        !$OMP END PARALLEL DO

        ! Apply prescribed boundary discharge to the final HP transport depths.
        ! This is intentionally after h_selected is finalized and before ghost copies.
        DO k=1,comp_cells_y
            IF (comp_cells_x>1) THEN
                IF (source_cell(1,k)/=2) THEN
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceR(:,1,k),                                   &
                        this%h_selected_interfaceR(1,k),qp_expl(:,1,k),                                               &
                        qp_expl(:,2,k),bcW(2),dx2,one_by_dx,-1,.TRUE.,idx_u)
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceR(:,1,k),                                   &
                        this%h_selected_interfaceR(1,k),qp_expl(:,1,k),                                               &
                        qp_expl(:,2,k),bcW(3),dx2,one_by_dx,-1,.FALSE.,idx_v)
                END IF
                IF (source_cell(comp_cells_x,k)/=2) THEN
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceL(:,comp_interfaces_x,k),                  &
                        this%h_selected_interfaceL(comp_interfaces_x,k),qp_expl(:,comp_cells_x,k),                    &
                        qp_expl(:,comp_cells_x-1,k),bcE(2),dx2,one_by_dx,1,.TRUE.,idx_u)
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceL(:,comp_interfaces_x,k),                  &
                        this%h_selected_interfaceL(comp_interfaces_x,k),qp_expl(:,comp_cells_x,k),                    &
                        qp_expl(:,comp_cells_x-1,k),bcE(3),dx2,one_by_dx,1,.FALSE.,idx_v)
                END IF
            END IF
        END DO
        DO j=1,comp_cells_x
            IF (comp_cells_y>1) THEN
                IF (source_cell(j,1)/=2) THEN
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceT(:,j,1),                                   &
                        this%h_selected_interfaceT(j,1),qp_expl(:,j,1),                                               &
                        qp_expl(:,j,2),bcS(2),dy2,one_by_dy,-1,.FALSE.,idx_u)
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceT(:,j,1),                                   &
                        this%h_selected_interfaceT(j,1),qp_expl(:,j,1),                                               &
                        qp_expl(:,j,2),bcS(3),dy2,one_by_dy,-1,.TRUE.,idx_v)
                END IF
                IF (source_cell(j,comp_cells_y)/=2) THEN
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceB(:,j,comp_interfaces_y),                  &
                        this%h_selected_interfaceB(j,comp_interfaces_y),qp_expl(:,j,comp_cells_y),                    &
                        qp_expl(:,j,comp_cells_y-1),bcN(2),dy2,one_by_dy,1,.FALSE.,idx_u)
                    CALL enforce_boundary_discharge(this%qp_adv_interfaceB(:,j,comp_interfaces_y),                  &
                        this%h_selected_interfaceB(j,comp_interfaces_y),qp_expl(:,j,comp_cells_y),                    &
                        qp_expl(:,j,comp_cells_y-1),bcN(3),dy2,one_by_dy,1,.TRUE.,idx_v)
                END IF
            END IF
        END DO
        ! Keep both sides of an external interface identical after the correction.
        DO k=1,comp_cells_y
            this%qp_adv_interfaceL(:,1,k)=this%qp_adv_interfaceR(:,1,k)
            this%qp_adv_interfaceR(:,comp_interfaces_x,k)=this%qp_adv_interfaceL(:,comp_interfaces_x,k)
        END DO
        DO j=1,comp_cells_x
            this%qp_adv_interfaceB(:,j,1)=this%qp_adv_interfaceT(:,j,1)
            this%qp_adv_interfaceT(:,j,comp_interfaces_y)=this%qp_adv_interfaceB(:,j,comp_interfaces_y)
        END DO

        ! External ghost states inherit the selected depths of their owners.
        DO k=1,comp_cells_y
            this%h_selected_interfaceL(1,k)=this%h_selected_interfaceR(1,k)
            this%h_selected_interfaceR(comp_interfaces_x,k)=                                                           &
                this%h_selected_interfaceL(comp_interfaces_x,k)

        END DO

        DO j=1,comp_cells_x
            this%h_selected_interfaceB(j,1)=this%h_selected_interfaceT(j,1)
            this%h_selected_interfaceT(j,comp_interfaces_y)=                                                           &
                this%h_selected_interfaceB(j,comp_interfaces_y)

        END DO
    END SUBROUTINE reconstruct_kt_eta_traces




        !> \brief Compute positivity-constrained free-surface traces for one parent-workset cell.
        !>
        !> \param[in] jc X index of the current cell in the parent reconstruction.
        !> \param[in] kc Y index of the current cell in the parent reconstruction.
        !>
        !> \note Reads the host qp_center and shared face bed; writes only this cell eta candidates
        !> in
        !>       the parent workspace.



        !> @brief Convert direct velocity traces to transport momenta using nonnegative depth.
        !> @details A dry trace has zero momentum and velocity; wet traces use
        !>          the same selected depth for both momentum components.
        !> @param[in,out] trace Primitive face trace modified in place.



    !> Convert an imposed discharge (Dirichlet or Neumann) into a primitive
    !> velocity boundary value, using the reconstructed first-pass depth.
    SUBROUTINE discharge_bc_to_velocity(discharge_bc,hc,uc,hn,un,hface,half,inv_dx,side,normal,velocity_bc)
        TYPE(bc), INTENT(IN) :: discharge_bc
        REAL(wp), INTENT(IN) :: hc,uc,hn,un,hface,half,inv_dx
        INTEGER, INTENT(IN) :: side
        LOGICAL, INTENT(IN) :: normal
        TYPE(bc), INTENT(OUT) :: velocity_bc
        REAL(wp) :: discharge
        CALL boundary_discharge_value(discharge_bc,hc*uc,hn*un,half,inv_dx,side,normal,discharge)
        velocity_bc%flag=0
        velocity_bc%value=0.0_wp
        IF (hface>dry_thickness_tolerance) THEN
            velocity_bc%value=discharge/hface
        ELSE
            ! The first-pass h trace is only a provisional limiter candidate.
            ! It can be dry even when the subsequent HP selection is wet.
            ! Use a zero provisional velocity rather than dividing by a dry
            ! depth; the physical hu/hv constraint is applied after HP.
            velocity_bc%value=0.0_wp
        END IF
    END SUBROUTINE discharge_bc_to_velocity

    !> Evaluate the boundary discharge in the original physical input units.
    SUBROUTINE boundary_discharge_value(discharge_bc,qc,qn,half,inv_dx,side,normal,discharge)
        TYPE(bc), INTENT(IN) :: discharge_bc
        REAL(wp), INTENT(IN) :: qc,qn,half,inv_dx
        INTEGER, INTENT(IN) :: side
        LOGICAL, INTENT(IN) :: normal
        REAL(wp), INTENT(OUT) :: discharge
        REAL(wp) :: gradient
        SELECT CASE(discharge_bc%flag)
        CASE(0)
            discharge=discharge_bc%value
        CASE(1)
            discharge=qc+side*reconstr_coeff*half*discharge_bc%value
        CASE(2)
            IF (side<0) THEN
                gradient=(qn-qc)*inv_dx
            ELSE
                gradient=(qc-qn)*inv_dx
            END IF
            discharge=qc+side*reconstr_coeff*half*gradient
        CASE DEFAULT
            ERROR STOP 'Invalid boundary discharge condition'
        END SELECT
        IF (normal .AND. discharge_bc%flag/=0) THEN
            IF (side<0) discharge=MIN(discharge,0.0_wp)
            IF (side>0) discharge=MAX(discharge,0.0_wp)
        END IF
    END SUBROUTINE boundary_discharge_value

    !> Enforce prescribed hu/hv after HP depth selection, preserving the
    !> original physical meaning of HU_BC/HV_BC (including Neumann gradients).
    SUBROUTINE enforce_boundary_discharge(trace,hselected,qcenter,qneighbor,discharge_bc,half,inv_dx,side,normal,idx)
        REAL(wp), INTENT(INOUT) :: trace(n_vars)
        REAL(wp), INTENT(IN) :: hselected,qcenter(n_vars),qneighbor(n_vars),half,inv_dx
        TYPE(bc), INTENT(IN) :: discharge_bc
        INTEGER, INTENT(IN) :: side,idx
        LOGICAL, INTENT(IN) :: normal
        REAL(wp) :: discharge
        CALL boundary_discharge_value(discharge_bc,qcenter(idx_h)*qcenter(idx),                                        &
            qneighbor(idx_h)*qneighbor(idx),half,inv_dx,side,normal,discharge)
        IF (hselected>dry_thickness_tolerance) THEN
            ! A wet HP face retains the prescribed discharge, including the
            ! original Neumann gradient in hu/hv rather than in u/v.
            trace(idx)=discharge/hselected
        ELSE
            ! An empty face cannot transport mass or momentum. For Neumann
            ! conditions the extrapolated discharge is not an imposed inflow:
            ! its admissible dry-face value is zero.
            ! A nonzero Dirichlet discharge is a genuine incompatible inlet
            ! prescription, so reject it instead of silently dropping mass.
            IF (discharge_bc%flag==0 .AND. ABS(discharge)>dry_thickness_tolerance)                                      &
                ERROR STOP 'Nonzero Dirichlet discharge at a dry HP boundary'
            trace(idx)=0.0_wp
        END IF
    END SUBROUTINE enforce_boundary_discharge

    !> @brief Select the admissible transport depth and finalize face fractions.
    !> @details Couples the depth choice to the component-fraction closure rather
    !>          than projecting those fields independently.
    !> @param[in,out] qp Face primitive state whose fractions and depth are finalized.
    !> @param[in] eta Reconstructed free-surface elevation [m].
    !> @param[in] bed Bed elevation at the face [m].
    !> @param[in] h_minmod Minmod transport-depth candidate [m].
    !> @param[in] h_froude Froude-limited transport-depth candidate [m].
    !> @param[in] source_trace Whether this face represents a source boundary.
    !> @param[out] h_selected Accepted nonnegative transport depth [m].
    !> @param[in] context Context string used in physical closure diagnostics.
    !> @param[in] kt_mode Selected KT trace reconstruction mode.
    SUBROUTINE select_final_transport_trace(qp,eta,bed,h_minmod,h_froude, source_trace,h_selected,context,kt_mode)
        USE state_composition, ONLY : normalize_face_composition
        USE diagnostics_2d, ONLY : fatal_error
        USE parameters_2d, ONLY : n_solid,n_add_gas,idx_solid_first,idx_solid_last, &
            idx_add_gas_first,idx_add_gas_last
        REAL(wp), INTENT(INOUT) :: qp(n_vars)
        REAL(wp), INTENT(IN) :: eta,bed,h_minmod,h_froude
        LOGICAL, INTENT(IN) :: source_trace
        REAL(wp), INTENT(OUT) :: h_selected
        CHARACTER(LEN=*), INTENT(IN) :: context
        INTEGER, INTENT(IN) :: kt_mode
        REAL(wp) :: xs(n_solid),xg(n_add_gas)

        SELECT CASE (kt_mode)
        CASE (1)
            h_selected=MAX(qp(idx_h),0.0_wp)
        CASE (3)
            h_selected=MAX(h_minmod,0.0_wp)
        CASE (4)
            h_selected=MAX(h_froude,0.0_wp)
        CASE DEFAULT
            h_selected=MAX(eta-bed,0.0_wp)

        END SELECT

        IF (source_trace) h_selected=MAX(qp(idx_h),0.0_wp)
        IF (h_selected .GT. dry_thickness_tolerance) THEN
            xs=qp(idx_solid_first:idx_solid_last)
            xg=qp(idx_add_gas_first:idx_add_gas_last)
            IF (.NOT.normalize_face_composition(xs,xg)) &
                CALL fatal_error(TRIM(context)//': Material negative face mass fraction')
            qp(idx_solid_first:idx_solid_last)=xs
            qp(idx_add_gas_first:idx_add_gas_last)=xg
        END IF
    END SUBROUTINE select_final_transport_trace

END MODULE reconstruction_2d
