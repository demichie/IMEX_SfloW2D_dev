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

  !> \brief Persistent wave speeds, oriented face contributions and cell paths.
  !> \details Face arrays distinguish the contribution seen by either adjacent
  !>          cell; hydrostatic pressure belongs to the path terms, not the
  !>          inertial endpoint flux.
  TYPE, PUBLIC :: hyperbolic_workspace_type
     ! One common characteristic bound per face and sign, not per equation.
     REAL(wp), ALLOCATABLE :: a_interface_xNeg(:,:)
     REAL(wp), ALLOCATABLE :: a_interface_xPos(:,:)
     REAL(wp), ALLOCATABLE :: a_interface_yNeg(:,:)
     REAL(wp), ALLOCATABLE :: a_interface_yPos(:,:)
     REAL(wp), ALLOCATABLE :: G_interface_xL(:,:,:)
     REAL(wp), ALLOCATABLE :: G_interface_xR(:,:,:)
     REAL(wp), ALLOCATABLE :: G_interface_yB(:,:,:)
     REAL(wp), ALLOCATABLE :: G_interface_yT(:,:,:)
     REAL(wp), ALLOCATABLE :: P_cell_x(:,:,:)
     REAL(wp), ALLOCATABLE :: P_cell_y(:,:,:)
   CONTAINS
     PROCEDURE :: initialize => initialize_hyperbolic
     PROCEDURE :: finalize => finalize_hyperbolic
     PROCEDURE :: evaluate_terms => eval_hyperbolic_terms
     PROCEDURE :: evaluate_speeds => eval_speeds
  END TYPE hyperbolic_workspace_type

CONTAINS

  !> \brief Allocate and zero wave-speed and oriented KT work arrays.
  !>
  !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.

  SUBROUTINE initialize_hyperbolic( this )

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this

    ALLOCATE( this%a_interface_xNeg(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%a_interface_xPos(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%a_interface_yNeg(comp_cells_x,comp_interfaces_y) )
    ALLOCATE( this%a_interface_yPos(comp_cells_x,comp_interfaces_y) )

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
    this%G_interface_xL = 0.0_wp
    this%G_interface_xR = 0.0_wp
    this%G_interface_yB = 0.0_wp
    this%G_interface_yT = 0.0_wp
    this%P_cell_x = 0.0_wp
    this%P_cell_y = 0.0_wp

  END SUBROUTINE initialize_hyperbolic

  !> \brief Release the hyperbolic workspace arrays.
  !>
  !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.

  SUBROUTINE finalize_hyperbolic( this )

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this

    DEALLOCATE( this%a_interface_xNeg )
    DEALLOCATE( this%a_interface_xPos )
    DEALLOCATE( this%a_interface_yNeg )
    DEALLOCATE( this%a_interface_yPos )
    DEALLOCATE( this%G_interface_xL )
    DEALLOCATE( this%G_interface_xR )
    DEALLOCATE( this%G_interface_yB )
    DEALLOCATE( this%G_interface_yT )
    DEALLOCATE( this%P_cell_x )
    DEALLOCATE( this%P_cell_y )

  END SUBROUTINE finalize_hyperbolic

  !> \brief Assemble the complete HP-PCCU spatial term on active cells.
  !>
  !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
  !> \param[in,out] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
  !> \param[in] qp_expl Cell-centered physical states at the explicit stage, including the required
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

  SUBROUTINE eval_hyperbolic_terms( this, recon, qp_expl,                     &
       divFlux_iRK, t, solve_cells, j_cent, k_cent, solve_interfaces_x,       &
       j_stag_x, k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y,           &
       source_boundary_enabled )

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(INOUT) :: recon
    REAL(wp), INTENT(IN) :: qp_expl(n_vars+2,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: divFlux_iRK(n_eqns,comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)
    LOGICAL, INTENT(IN), OPTIONAL :: source_boundary_enabled

    INTEGER :: l , i, j, k      !< loop counters

    ! Linear reconstruction of the physical variables at the interfaces
    CALL recon%reconstruct( qp_expl, t, solve_cells, j_cent, k_cent, .TRUE., &
         source_boundary_enabled=source_boundary_enabled )

    CALL recon%reconstruct_kt_eta(qp_expl)

    ! Evaluation of the maximum local speeds at the interfaces
    CALL this%evaluate_speeds( recon, solve_interfaces_x, j_stag_x,          &
         k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    ! The production spatial operator is the single HP-PCCU baseline.
    CALL eval_flux_PCCU( this, recon, qp_expl,                              &
         .TRUE., t,                                                         &
         solve_interfaces_x, j_stag_x,                                       &
         k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    CALL eval_cell_hydrostatic_paths( this, recon, qp_expl,                  &
         .TRUE., solve_cells, j_cent, k_cent )

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

             divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                       &
                  ( this%G_interface_xL(i,j+1,k)                             &
                  - this%G_interface_xR(i,j,k)                              &
                  - this%P_cell_x(i,j,k) ) * one_by_dx

          END IF

          IF ( comp_cells_y .GT. 1 ) THEN

             divFlux_iRK(i,j,k) = divFlux_iRK(i,j,k) +                       &
                  ( this%G_interface_yB(i,j,k+1)                             &
                  - this%G_interface_yT(i,j,k)                              &
                  - this%P_cell_y(i,j,k) ) * one_by_dy

          END IF

       END DO

       ! Only a completely stationary stencil may discard hydrostatic
       ! cancellation noise. Never classify a nonzero physical velocity as
       ! rest: even tiny motion must retain the unmodified spatial operator.
       ! A pressure-scaled roundoff residual otherwise creates a tiny discharge
       ! which disables the exact-rest scalar-flux guard at the next stage.
       IF ( exactly_resting_stencil(recon,qp_expl(:,j,k),j,k) ) THEN
          IF ( comp_cells_x .GT. 1 ) CALL cancel_resting_roundoff(             &
               recon,j,k,1,divFlux_iRK(2,j,k) )
          IF ( comp_cells_y .GT. 1 ) CALL cancel_resting_roundoff(             &
               recon,j,k,2,divFlux_iRK(3,j,k) )
       END IF

    END DO cells_loop

    !$OMP END PARALLEL DO

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
       resting = ALL(recon%qp_adv_interfaceL(2:3,j:j+1,k) .EQ. 0.0_wp) .AND. &
            ALL(recon%qp_adv_interfaceR(2:3,j:j+1,k) .EQ. 0.0_wp)
       IF ( .NOT. resting ) RETURN
    END IF
    IF ( comp_cells_y .GT. 1 ) THEN
       resting = ALL(recon%qp_adv_interfaceB(2:3,j,k:k+1) .EQ. 0.0_wp) .AND. &
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
  SUBROUTINE cancel_resting_roundoff(recon,j,k,axis,residual)
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    INTEGER, INTENT(IN) :: j,k,axis
    REAL(wp), INTENT(INOUT) :: residual
    REAL(wp) :: scale, inverse_spacing

    IF (axis .EQ. 1) THEN
       scale = MAX(                                                      &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceL(:,j,k),       &
            recon%eta_interfaceL(j,k),grav_coeff_stag_x(j,k)),             &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceR(:,j,k),       &
            recon%eta_interfaceR(j,k),grav_coeff_stag_x(j,k)),             &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceL(:,j+1,k),     &
            recon%eta_interfaceL(j+1,k),grav_coeff_stag_x(j+1,k)),          &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceR(:,j+1,k),     &
            recon%eta_interfaceR(j+1,k),grav_coeff_stag_x(j+1,k)) )
       inverse_spacing = one_by_dx
    ELSE
       scale = MAX(                                                      &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceB(:,j,k),       &
            recon%eta_interfaceB(j,k),grav_coeff_stag_y(j,k)),             &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceT(:,j,k),       &
            recon%eta_interfaceT(j,k),grav_coeff_stag_y(j,k)),             &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceB(:,j,k+1),     &
            recon%eta_interfaceB(j,k+1),grav_coeff_stag_y(j,k+1)),          &
            hydrostatic_arithmetic_scale(recon%qp_adv_interfaceT(:,j,k+1),     &
            recon%eta_interfaceT(j,k+1),grav_coeff_stag_y(j,k+1)) )
       inverse_spacing = one_by_dy
    END IF
    IF ( ABS(residual) .LE. 64.0_wp*EPSILON(1.0_wp)*scale*inverse_spacing ) &
         residual = 0.0_wp
  END SUBROUTINE cancel_resting_roundoff

  !> \brief Estimate the pressure scale underlying one hydrostatic endpoint calculation.
  !> \param[in] qp_face Reconstructed physical state, with depth in entry one.
  !> \param[in] eta Face free surface [m], included for subtraction roundoff.
  !> \param[in] gravity_factor Dimensionless large-slope factor at the face.
  !> \return Absolute path arithmetic scale [kg s^-2]; zero for a dry endpoint.
  FUNCTION hydrostatic_arithmetic_scale(qp_face,eta,gravity_factor) RESULT(scale)
    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    REAL(wp), INTENT(IN) :: qp_face(:),eta,gravity_factor
    REAL(wp) :: scale,gamma,reduced_gravity,h

    h = MAX(0.0_wp,qp_face(1))
    scale = 0.0_wp
    IF (h .EQ. 0.0_wp) RETURN
    CALL eval_hydrostatic_coefficient(qp_face,reduced_gravity,gamma)
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

  SUBROUTINE eval_flux_PCCU( this, recon, qp_center, use_dry_donor, t,       &
       solve_interfaces_x,                                                   &
       j_stag_x,                                                             &
       k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    REAL(wp), INTENT(IN) :: qp_center(n_vars+2,comp_cells_x,comp_cells_y)
    LOGICAL, INTENT(IN) :: use_dry_donor
    REAL(wp), INTENT(IN) :: t
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: a_minus, a_plus
    LOGICAL :: dry_donor_left,dry_donor_right
    INTEGER :: j, k, l

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,a_minus,a_plus,dry_donor_left,dry_donor_right)
       DO l = 1, solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          a_minus = this%a_interface_xNeg(j,k)
          a_plus = this%a_interface_xPos(j,k)
          dry_donor_left = recon%qp_adv_interfaceL(1,j,k) .LE. dry_thickness_tolerance
          dry_donor_right = recon%qp_adv_interfaceR(1,j,k) .LE. dry_thickness_tolerance
          IF ( j .GT. 1 ) dry_donor_left =                                  &
               qp_center(1,j-1,k) .LE. dry_thickness_tolerance
          IF ( j .LE. comp_cells_x ) dry_donor_right =                      &
               qp_center(1,j,k) .LE. dry_thickness_tolerance
          CALL eval_pccu_interface( recon%eta_interfaceL(j,k),              &
               recon%eta_interfaceR(j,k), B_face_x(j,k),                     &
                recon%h_selected_interfaceL(j,k),                              &
                recon%h_selected_interfaceR(j,k),                              &
               recon%qp_adv_interfaceL(:,j,k),                                &
               recon%qp_adv_interfaceR(:,j,k), grav_coeff_stag_x(j,k),       &
               a_minus, a_plus, PATH_DIR_X, 2, j, k, use_dry_donor,         &
               dry_donor_left,dry_donor_right,t,                            &
               this%G_interface_xL(:,j,k), this%G_interface_xR(:,j,k) )

          ! Apply the boundary's one-way condition to the final numerical flux.
          ! V2 limited the reconstructed normal momentum to outward-only values
          ! (west/south <= 0; east/north >= 0). HP reconstruction follows the
          ! primitive boundary trace and can otherwise reintroduce inflow, so
          ! remove incoming mass-carried transport here. Keep normal momentum
          ! flux, which contains the wall pressure contribution. Source cells
          ! are exempt, as in the legacy reconstruction.
          IF ( ( j .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),  &
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                    &
               ( bcW(2)%flag .NE. 0 ) .AND.                                  &
               ( this%G_interface_xL(1,j,k) .GT. 0.0_wp ) ) THEN
             this%G_interface_xL(1,j,k) = 0.0_wp
             this%G_interface_xR(1,j,k) = 0.0_wp
             this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
             this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( j .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),&
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                     &
               ( bcW(2)%flag .EQ. 0 ) .AND.                                  &
               ( bcW(2)%value .EQ. 0.0_wp ) ) THEN
             this%G_interface_xL(1,j,k) = 0.0_wp
             this%G_interface_xR(1,j,k) = 0.0_wp
             this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
             this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( j .EQ. comp_interfaces_x ) .AND.                        &
               ( source_cell(MAX(1,MIN(j-1,SIZE(source_cell,1))),             &
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                    &
               ( bcE(2)%flag .NE. 0 ) .AND.                                  &
               ( this%G_interface_xL(1,j,k) .LT. 0.0_wp ) ) THEN
             this%G_interface_xL(1,j,k) = 0.0_wp
             this%G_interface_xR(1,j,k) = 0.0_wp
             this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
             this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( j .EQ. comp_interfaces_x ) .AND.                        &
               ( source_cell(MAX(1,MIN(j-1,SIZE(source_cell,1))),             &
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                    &
               ( bcE(2)%flag .EQ. 0 ) .AND.                                  &
               ( bcE(2)%value .EQ. 0.0_wp ) ) THEN
             this%G_interface_xL(1,j,k) = 0.0_wp
             this%G_interface_xR(1,j,k) = 0.0_wp
             this%G_interface_xL(3:n_eqns,j,k) = 0.0_wp
             this%G_interface_xR(3:n_eqns,j,k) = 0.0_wp
          END IF

       END DO
       !$OMP END PARALLEL DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,a_minus,a_plus,dry_donor_left,dry_donor_right)
       DO l = 1, solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          a_minus = this%a_interface_yNeg(j,k)
          a_plus = this%a_interface_yPos(j,k)
          dry_donor_left = recon%qp_adv_interfaceB(1,j,k) .LE. dry_thickness_tolerance
          dry_donor_right = recon%qp_adv_interfaceT(1,j,k) .LE. dry_thickness_tolerance
          IF ( k .GT. 1 ) dry_donor_left =                                  &
               qp_center(1,j,k-1) .LE. dry_thickness_tolerance
          IF ( k .LE. comp_cells_y ) dry_donor_right =                      &
               qp_center(1,j,k) .LE. dry_thickness_tolerance
          CALL eval_pccu_interface( recon%eta_interfaceB(j,k),              &
               recon%eta_interfaceT(j,k), B_face_y(j,k),                     &
                recon%h_selected_interfaceB(j,k),                              &
                recon%h_selected_interfaceT(j,k),                              &
               recon%qp_adv_interfaceB(:,j,k),                                &
               recon%qp_adv_interfaceT(:,j,k), grav_coeff_stag_y(j,k),       &
               a_minus, a_plus, PATH_DIR_Y, 3, j, k, use_dry_donor,         &
               dry_donor_left,dry_donor_right,t,                            &
               this%G_interface_yB(:,j,k), this%G_interface_yT(:,j,k) )

          IF ( ( k .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),  &
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                    &
               ( bcS(3)%flag .NE. 0 ) .AND.                                  &
               ( this%G_interface_yB(1,j,k) .GT. 0.0_wp ) ) THEN
             this%G_interface_yB(1:2,j,k) = 0.0_wp
             this%G_interface_yT(1:2,j,k) = 0.0_wp
             this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
             this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( k .EQ. 1 ) .AND. ( source_cell(MIN(j,SIZE(source_cell,1)),&
               MIN(k,SIZE(source_cell,2))) .NE. 2 ) .AND.                     &
               ( bcS(3)%flag .EQ. 0 ) .AND.                                 &
               ( bcS(3)%value .EQ. 0.0_wp ) ) THEN
             this%G_interface_yB(1:2,j,k) = 0.0_wp
             this%G_interface_yT(1:2,j,k) = 0.0_wp
             this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
             this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( k .EQ. comp_interfaces_y ) .AND.                        &
               ( source_cell(MIN(j,SIZE(source_cell,1)),                      &
               MAX(1,MIN(k-1,SIZE(source_cell,2)))) .NE. 2 ) .AND.           &
               ( bcN(3)%flag .NE. 0 ) .AND.                                  &
               ( this%G_interface_yB(1,j,k) .LT. 0.0_wp ) ) THEN
             this%G_interface_yB(1:2,j,k) = 0.0_wp
             this%G_interface_yT(1:2,j,k) = 0.0_wp
             this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
             this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
          ELSEIF ( ( k .EQ. comp_interfaces_y ) .AND.                        &
               ( source_cell(MIN(j,SIZE(source_cell,1)),                      &
               MAX(1,MIN(k-1,SIZE(source_cell,2)))) .NE. 2 ) .AND.           &
               ( bcN(3)%flag .EQ. 0 ) .AND.                                  &
               ( bcN(3)%value .EQ. 0.0_wp ) ) THEN
             this%G_interface_yB(1:2,j,k) = 0.0_wp
             this%G_interface_yT(1:2,j,k) = 0.0_wp
             this%G_interface_yB(4:n_eqns,j,k) = 0.0_wp
             this%G_interface_yT(4:n_eqns,j,k) = 0.0_wp
          END IF

       END DO
       !$OMP END PARALLEL DO

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
  !> \param[in] h_minmod_left Minmod-reconstructed depth on the negative side.
  !> \param[in] h_minmod_right Minmod-reconstructed depth on the positive side.
  !> \param[in] h_froude_left Froude-blend depth on the negative side.
  !> \param[in] h_froude_right Froude-blend depth on the positive side.
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
  !> \param[out] value_right PCCU contribution assigned to the positive-side cell.
  SUBROUTINE eval_pccu_interface(eta_left,eta_right,bed_face,                &
       h_selected_left,h_selected_right,                                      &
       qp_adv_left,qp_adv_right,                                              &
       grav_coeff,a_minus,a_plus,direction,normal_momentum_index,face_i,     &
       face_j,use_dry_donor,dry_donor_left,dry_donor_right,t,                &
       value_left,value_right)

    USE equation_terms_2d, ONLY : eval_inertial_flux
    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    USE equation_terms_2d, ONLY : limit_component_mass_flux
    USE state_conversion_2d, ONLY : qp_to_qc
    USE pccu_2d, ONLY : eval_hydrostatic_path
    USE pccu_2d, ONLY : eval_oriented_pccu_pair
    USE pccu_2d, ONLY : pccu_speed_tolerance
    USE nonconservative_2d, ONLY : PATH_DIR_X

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: eta_left, eta_right, bed_face
    REAL(wp), INTENT(IN) :: h_selected_left,h_selected_right
    REAL(wp), INTENT(IN) :: qp_adv_left(:), qp_adv_right(:)
    REAL(wp), INTENT(IN) :: grav_coeff, a_minus, a_plus
    INTEGER, INTENT(IN) :: direction, normal_momentum_index
    INTEGER, INTENT(IN) :: face_i, face_j
    LOGICAL, INTENT(IN) :: use_dry_donor
    LOGICAL, INTENT(IN) :: dry_donor_left,dry_donor_right
    REAL(wp), INTENT(IN) :: t
    REAL(wp), INTENT(OUT) :: value_left(:), value_right(:)

    REAL(wp) :: q_adv_left_use(n_vars), q_adv_right_use(n_vars)
    REAL(wp) :: qp_adv_left_use(n_vars+2), qp_adv_right_use(n_vars+2)
    REAL(wp) :: flux_left(n_eqns), flux_right(n_eqns)
    REAL(wp) :: path_contribution(n_eqns), interface_flux(n_eqns)
    REAL(wp) :: qp_path_left(n_vars+2), qp_path_right(n_vars+2)
    REAL(wp) :: gamma_left, gamma_right
    REAL(wp) :: reduced_gravity_left, reduced_gravity_right
    REAL(wp) :: eta_scale, eta_roundoff
    REAL(wp) :: denominator,u_normal_left,u_normal_right
    REAL(wp) :: eta_path_left, eta_path_right
    REAL(wp) :: h_eta_left,h_eta_right
    ! Constant contexts avoid expensive formatted I/O in each interface.

    qp_adv_left_use=qp_adv_left
    qp_adv_right_use=qp_adv_right
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
       qp_adv_left_use(2)=h_eta_left*qp_adv_left(idx_u)
       qp_adv_left_use(3)=h_eta_left*qp_adv_left(idx_v)
    END IF
    IF ( h_eta_right .LE. dry_thickness_tolerance ) THEN
       qp_adv_right_use(1:3)=0.0_wp
       qp_adv_right_use(idx_u:idx_v)=0.0_wp
    ELSE
       qp_adv_right_use(2)=h_eta_right*qp_adv_right(idx_u)
       qp_adv_right_use(3)=h_eta_right*qp_adv_right(idx_v)
    END IF
    CALL qp_to_qc(qp_adv_left_use,q_adv_left_use,                        &
         'KT interface left advective trace')
    CALL qp_to_qc(qp_adv_right_use,q_adv_right_use,                      &
         'KT interface right advective trace')

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
    denominator=a_plus-a_minus
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
    IF ( denominator .GT. pccu_speed_tolerance ) THEN
       interface_flux=(a_plus*(u_normal_left-a_minus)*q_adv_left_use          &
            -a_minus*(u_normal_right-a_plus)*q_adv_right_use)/denominator
    ELSE
       interface_flux=0.5_wp*(flux_left+flux_right)
    END IF

    CALL limit_component_mass_flux(interface_flux)

    ! Positivity guard: when the total mass flux points out of a dry donor,
    ! suppress the complete inertial flux vector, not only the mass equation.
    ! Keep the nonconservative hydrostatic path contribution below untouched.
    IF ( use_dry_donor ) THEN
       IF ( (interface_flux(1) .GT. 0.0_wp) .AND.                         &
            (dry_donor_left .OR.                                         &
            (qp_adv_left_use(1) .LE. dry_thickness_tolerance)) ) THEN
          interface_flux=0.0_wp
       ELSEIF ( (interface_flux(1) .LT. 0.0_wp) .AND.                    &
            (dry_donor_right .OR.                                        &
            (qp_adv_right_use(1) .LE. dry_thickness_tolerance)) ) THEN
          interface_flux=0.0_wp
       END IF
    END IF

    ! Suppress only the roundoff-level mass diffusion of a hydrostatic rest
    ! state. Zero normal velocity alone is not enough: on a flat bed, for
    ! example, a wet cell next to dry cells has unequal free-surface levels
    ! and must be allowed to drain. At an excavation crest, HP reconstruction
    ! makes both states dry (so no uphill flux is introduced), while equal-eta
    ! resting traces remain protected from a spurious central-upwind flux.
    eta_scale = MAX(1.0_wp,ABS(eta_left),ABS(eta_right))
    eta_roundoff = 8.0_wp*EPSILON(1.0_wp)*eta_scale
    IF ( ( qp_adv_left(normal_momentum_index) .EQ. 0.0_wp ) .AND.           &
         ( qp_adv_right(normal_momentum_index) .EQ. 0.0_wp ) .AND.          &
         ( ABS(eta_left-eta_right) .LE. eta_roundoff ) ) THEN
       interface_flux(1) = 0.0_wp
       interface_flux(4:n_vars) = 0.0_wp
    END IF

    qp_path_left=qp_adv_left_use
    qp_path_right=qp_adv_right_use
    eta_path_left=eta_left
    eta_path_right=eta_right
    IF ( source_trace_left ) eta_path_left=bed_face+h_eta_left
    IF ( source_trace_right ) eta_path_right=bed_face+h_eta_right
    IF ( use_dry_donor ) THEN
       IF ( qp_adv_left_use(1) .LE. dry_thickness_tolerance ) THEN
          qp_path_left(1)=0.0_wp
          eta_path_left=bed_face
       END IF
       IF ( qp_adv_right_use(1) .LE. dry_thickness_tolerance ) THEN
          qp_path_right(1)=0.0_wp
          eta_path_right=bed_face
       END IF
    END IF
    CALL eval_hydrostatic_coefficient(qp_path_left,reduced_gravity_left,gamma_left)
    CALL eval_hydrostatic_coefficient(qp_path_right,reduced_gravity_right,gamma_right)
    CALL regularize_dry_gamma_pair(qp_path_left(1),gamma_left,               &
         qp_path_right(1),gamma_right)
    CALL eval_hydrostatic_path(direction,qp_path_left(1),gamma_left,eta_path_left,&
         grav_coeff,qp_path_right(1),gamma_right,eta_path_right,grav_coeff,  &
         path_contribution)

    CALL eval_oriented_pccu_pair(interface_flux,path_contribution,a_minus,   &
         a_plus,value_left,value_right)

    ! A dry adjacent column receives no share of the nonconservative path;
    ! the wet side retains its own oriented contribution.
    IF ( use_dry_donor ) THEN
       IF ( qp_adv_left_use(1) .LE. dry_thickness_tolerance ) value_left=interface_flux
       IF ( qp_adv_right_use(1) .LE. dry_thickness_tolerance ) value_right=interface_flux
    END IF

  END SUBROUTINE eval_pccu_interface

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

  SUBROUTINE eval_cell_hydrostatic_paths( this, recon, qp_expl, use_dry_donor,&
       solve_cells, j_cent, k_cent )

    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    REAL(wp), INTENT(IN) :: qp_expl(n_vars+2,comp_cells_x,comp_cells_y)
    LOGICAL, INTENT(IN) :: use_dry_donor
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)

    INTEGER :: j, k, l
    REAL(wp) :: qp_path_minus(n_vars+2), qp_path_plus(n_vars+2)

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
          qp_path_minus(1)=MAX(recon%eta_cellW(j,k)-B_face_x(j,k),0.0_wp)
          qp_path_plus(1)=MAX(recon%eta_cellE(j,k)-B_face_x(j+1,k),0.0_wp)
          IF (qp_path_minus(1).LE.dry_thickness_tolerance) qp_path_minus(2:3)=0.0_wp
          IF (qp_path_plus(1).LE.dry_thickness_tolerance) qp_path_plus(2:3)=0.0_wp
          CALL eval_cell_hydrostatic_path( qp_path_minus,                    &
               recon%eta_cellW(j,k),grav_coeff_stag_x(j,k),                 &
               qp_path_plus,recon%eta_cellE(j,k),                            &
               grav_coeff_stag_x(j+1,k),PATH_DIR_X,this%P_cell_x(:,j,k) )
       END IF

       IF ( comp_cells_y .GT. 1 ) THEN
          qp_path_minus=recon%qp_adv_interfaceT(:,j,k)
          qp_path_plus=recon%qp_adv_interfaceB(:,j,k+1)
          qp_path_minus(1)=MAX(recon%eta_cellS(j,k)-B_face_y(j,k),0.0_wp)
          qp_path_plus(1)=MAX(recon%eta_cellN(j,k)-B_face_y(j,k+1),0.0_wp)
          IF (qp_path_minus(1).LE.dry_thickness_tolerance) qp_path_minus(2:3)=0.0_wp
          IF (qp_path_plus(1).LE.dry_thickness_tolerance) qp_path_plus(2:3)=0.0_wp
          CALL eval_cell_hydrostatic_path( qp_path_minus,                    &
               recon%eta_cellS(j,k),grav_coeff_stag_y(j,k),                 &
               qp_path_plus,recon%eta_cellN(j,k),                            &
               grav_coeff_stag_y(j,k+1),PATH_DIR_Y,this%P_cell_y(:,j,k) )
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
  SUBROUTINE eval_cell_hydrostatic_path(qp_minus,eta_minus,grav_minus,       &
       qp_plus,eta_plus,grav_plus,direction,path)

    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    USE pccu_2d, ONLY : eval_hydrostatic_path

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qp_minus(:),eta_minus,grav_minus
    REAL(wp), INTENT(IN) :: qp_plus(:),eta_plus,grav_plus
    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(OUT) :: path(:)

    REAL(wp) :: gamma_minus,gamma_plus
    REAL(wp) :: reduced_gravity_minus,reduced_gravity_plus

    CALL eval_hydrostatic_coefficient(qp_minus,reduced_gravity_minus,       &
         gamma_minus)
    CALL eval_hydrostatic_coefficient(qp_plus,reduced_gravity_plus,gamma_plus)
    CALL regularize_dry_gamma_pair(qp_minus(1),gamma_minus,qp_plus(1),       &
         gamma_plus)
    CALL eval_hydrostatic_path(direction,qp_minus(1),gamma_minus,eta_minus, &
         grav_minus,qp_plus(1),gamma_plus,eta_plus,grav_plus,path)

  END SUBROUTINE eval_cell_hydrostatic_path

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
  !> \param[in,out] gamma_left Left coefficient, replaced by the right one if only the left trace is
  !>                           dry.
  !> \param[in] h_right Depth at the positive/right path endpoint [m].
  !> \param[in,out] gamma_right Right coefficient, replaced by the left one if only the right trace
  !>                            is dry.
  !******************************************************************************

  PURE SUBROUTINE regularize_dry_gamma_pair(h_left,gamma_left,h_right,gamma_right)
    REAL(wp), INTENT(IN) :: h_left,h_right
    REAL(wp), INTENT(INOUT) :: gamma_left,gamma_right
    IF ( ( h_left .LE. dry_thickness_tolerance ) .AND.                              &
         ( h_right .GT. dry_thickness_tolerance ) ) THEN
       gamma_left = gamma_right
    ELSEIF ( ( h_right .LE. dry_thickness_tolerance ) .AND.                         &
         ( h_left .GT. dry_thickness_tolerance ) ) THEN
       gamma_right = gamma_left
    END IF
  END SUBROUTINE regularize_dry_gamma_pair

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

  SUBROUTINE eval_speeds( this, recon, solve_interfaces_x, j_stag_x,         &
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
    REAL(wp) :: qp_eta_minus(n_vars+2), qp_eta_plus(n_vars+2)

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
          IF (h_eta_minus .LE. dry_thickness_tolerance) &
               qp_eta_minus(idx_u:idx_v)=0.0_wp
          IF (h_eta_plus .LE. dry_thickness_tolerance) &
               qp_eta_plus(idx_u:idx_v)=0.0_wp
          qp_eta_minus(2)=h_eta_minus*qp_eta_minus(idx_u)
          qp_eta_minus(3)=h_eta_minus*qp_eta_minus(idx_v)
          qp_eta_plus(2)=h_eta_plus*qp_eta_plus(idx_u)
          qp_eta_plus(3)=h_eta_plus*qp_eta_plus(idx_v)
          CALL eval_interface_speed_bounds(qp_eta_minus,qp_eta_plus, &
               grav_coeff_stag_x(j,k),idx_u,min_r,max_r)

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
          IF (h_eta_minus .LE. dry_thickness_tolerance) &
               qp_eta_minus(idx_u:idx_v)=0.0_wp
          IF (h_eta_plus .LE. dry_thickness_tolerance) &
               qp_eta_plus(idx_u:idx_v)=0.0_wp
          qp_eta_minus(2)=h_eta_minus*qp_eta_minus(idx_u)
          qp_eta_minus(3)=h_eta_minus*qp_eta_minus(idx_v)
          qp_eta_plus(2)=h_eta_plus*qp_eta_plus(idx_u)
          qp_eta_plus(3)=h_eta_plus*qp_eta_plus(idx_v)
          CALL eval_interface_speed_bounds(qp_eta_minus,qp_eta_plus, &
               grav_coeff_stag_y(j,k),idx_v,min_r,max_r)

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


  SUBROUTINE eval_interface_speed_bounds(qp_minus,qp_plus,grav_coeff,       &
       normal_velocity_index,speed_minus,speed_plus)

    USE equation_terms_2d, ONLY : eval_local_speeds_direction

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qp_minus(:),qp_plus(:),grav_coeff
    INTEGER, INTENT(IN) :: normal_velocity_index
    REAL(wp), INTENT(OUT) :: speed_minus,speed_plus

    REAL(wp) :: minus_min,minus_max,plus_min,plus_max

    CALL eval_local_speeds_direction(qp_minus,grav_coeff,                   &
         normal_velocity_index,minus_min,minus_max)
    CALL eval_local_speeds_direction(qp_plus,grav_coeff,                    &
         normal_velocity_index,plus_min,plus_max)

    speed_minus = MIN(minus_min,plus_min,0.0_wp)
    speed_plus = MAX(minus_max,plus_max,0.0_wp)

  END SUBROUTINE eval_interface_speed_bounds

  !******************************************************************************
  !> \brief Identify lateral source traces on the domain boundary.
  !> \param[in] direction Cartesian orientation of the face.
  !> \param[in] face_i X interface index.
  !> \param[in] face_j Y interface index.
  !> \param[out] source_left True when the negative-side trace is a source cell.
  !> \param[out] source_right True when the positive-side trace is a source cell.
  !******************************************************************************
  SUBROUTINE lateral_source_trace_sides(direction,face_i,face_j,           &
       source_left,source_right)

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
       IF ( face_i .EQ. comp_interfaces_x ) &
            source_right=source_cell(comp_cells_x,face_j) .EQ. 2
    CASE (PATH_DIR_Y)
       IF ( face_j .EQ. 1 ) source_left=source_cell(face_i,1) .EQ. 2
       IF ( face_j .EQ. comp_interfaces_y ) &
            source_right=source_cell(face_i,comp_cells_y) .EQ. 2
    END SELECT

  END SUBROUTINE lateral_source_trace_sides

END MODULE hyperbolic_2d
