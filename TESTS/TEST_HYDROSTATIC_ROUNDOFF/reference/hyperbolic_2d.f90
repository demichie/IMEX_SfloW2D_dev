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

  USE parameters_2d, ONLY : wp, n_eqns, n_vars
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : grav_coeff_stag_x, grav_coeff_stag_y
  USE geometry_2d, ONLY : one_by_dx, one_by_dy

  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hp_reconstruction_2d, ONLY : hp_dry_tolerance

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

  !> \brief Allocate and zero wave-speed and oriented HP-PCCU work arrays.
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

  SUBROUTINE eval_hyperbolic_terms( this, recon, qp_expl,                     &
       divFlux_iRK, t, solve_cells, j_cent, k_cent, solve_interfaces_x,       &
       j_stag_x, k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

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

    INTEGER :: l , i, j, k      !< loop counters

    !WRITE(*,*) 'SUBROUTINE eval_hyperbolic_terms'
    !WRITE(*,*) 'qp_expl(4,1,1)',qp_expl(4,1,1)
    !WRITE(*,*)
    
    ! Linear reconstruction of the physical variables at the interfaces
     CALL recon%reconstruct( qp_expl, t, solve_cells, &
         j_cent, k_cent )

    ! Evaluation of the maximum local speeds at the interfaces
    CALL this%evaluate_speeds( recon, solve_interfaces_x, j_stag_x,          &
         k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    ! The production spatial operator is the single HP-PCCU baseline.
    CALL eval_flux_PCCU( this, recon, solve_interfaces_x, j_stag_x,          &
         k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    CALL eval_cell_hydrostatic_paths( this, recon, solve_cells, j_cent, k_cent )

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

    END DO cells_loop

    !$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE eval_hyperbolic_terms

  !******************************************************************************
  !> \brief Build both oriented HP-PCCU values at every active face.
  !>
  !> The conservative endpoint flux is purely inertial.  Hydrostatic pressure
  !> and bed geometry enter only through the path contribution used to form the
  !> two oriented values at each face.
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

  SUBROUTINE eval_flux_PCCU( this, recon, solve_interfaces_x, j_stag_x,       &
       k_stag_x, solve_interfaces_y, j_stag_y, k_stag_y )

    USE equation_terms_2d, ONLY : eval_inertial_flux
    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    USE equation_terms_2d, ONLY : limit_component_mass_flux
    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y
    USE pccu_2d, ONLY : eval_hydrostatic_path
    USE pccu_2d, ONLY : eval_central_upwind_flux
    USE pccu_2d, ONLY : eval_oriented_pccu_pair

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: flux_left(n_eqns), flux_right(n_eqns)
    REAL(wp) :: path_contribution(n_eqns), interface_flux(n_eqns)
    REAL(wp) :: gamma_left, gamma_right
    REAL(wp) :: reduced_gravity_left, reduced_gravity_right
    REAL(wp) :: a_minus, a_plus
    INTEGER :: j, k, l

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,flux_left,flux_right,path_contribution,&
       !$OMP & interface_flux,                                               &
       !$OMP & gamma_left,gamma_right,reduced_gravity_left,                   &
       !$OMP & reduced_gravity_right,a_minus,a_plus)
       DO l = 1, solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_inertial_flux( recon%q_interfaceL(:,j,k),                &
               recon%qp_interfaceL(:,j,k), PATH_DIR_X, flux_left )
          CALL eval_inertial_flux( recon%q_interfaceR(:,j,k),                &
               recon%qp_interfaceR(:,j,k), PATH_DIR_X, flux_right )

          a_minus = this%a_interface_xNeg(j,k)
          a_plus = this%a_interface_xPos(j,k)
          CALL eval_central_upwind_flux( a_minus, a_plus, flux_left,         &
               flux_right, recon%q_interfaceL(:,j,k),                        &
               recon%q_interfaceR(:,j,k), interface_flux )

          CALL limit_component_mass_flux(interface_flux)

          ! Match the frozen baseline: a face at exact rest transports neither
          ! total mass nor thermal/composition scalars.
          IF ( ( recon%qp_interfaceL(2,j,k) .EQ. 0.0_wp ) .AND.             &
               ( recon%qp_interfaceR(2,j,k) .EQ. 0.0_wp ) ) THEN
             interface_flux(1) = 0.0_wp
             interface_flux(4:n_vars) = 0.0_wp
          END IF

          CALL eval_hydrostatic_coefficient( recon%qp_interfaceL(:,j,k),     &
               reduced_gravity_left, gamma_left )
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceR(:,j,k),     &
               reduced_gravity_right, gamma_right )
          CALL regularize_dry_gamma_pair( recon%qp_interfaceL(1,j,k), gamma_left, &
               recon%qp_interfaceR(1,j,k), gamma_right )
          CALL eval_hydrostatic_path( PATH_DIR_X,                            &
               recon%qp_interfaceL(1,j,k), gamma_left,                       &
               recon%eta_interfaceL(j,k), grav_coeff_stag_x(j,k),            &
               recon%qp_interfaceR(1,j,k), gamma_right,                      &
               recon%eta_interfaceR(j,k), grav_coeff_stag_x(j,k),            &
               path_contribution )

          CALL eval_oriented_pccu_pair( interface_flux,                       &
               path_contribution, a_minus, a_plus,                           &
               this%G_interface_xL(:,j,k), this%G_interface_xR(:,j,k) )

       END DO
       !$OMP END PARALLEL DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP PARALLEL DO private(l,j,k,flux_left,flux_right,path_contribution,&
       !$OMP & interface_flux,                                               &
       !$OMP & gamma_left,gamma_right,reduced_gravity_left,                   &
       !$OMP & reduced_gravity_right,a_minus,a_plus)
       DO l = 1, solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_inertial_flux( recon%q_interfaceB(:,j,k),                &
               recon%qp_interfaceB(:,j,k), PATH_DIR_Y, flux_left )
          CALL eval_inertial_flux( recon%q_interfaceT(:,j,k),                &
               recon%qp_interfaceT(:,j,k), PATH_DIR_Y, flux_right )

          a_minus = this%a_interface_yNeg(j,k)
          a_plus = this%a_interface_yPos(j,k)
          CALL eval_central_upwind_flux( a_minus, a_plus, flux_left,         &
               flux_right, recon%q_interfaceB(:,j,k),                        &
               recon%q_interfaceT(:,j,k), interface_flux )

          CALL limit_component_mass_flux(interface_flux)

          IF ( ( recon%qp_interfaceB(3,j,k) .EQ. 0.0_wp ) .AND.             &
               ( recon%qp_interfaceT(3,j,k) .EQ. 0.0_wp ) ) THEN
             interface_flux(1) = 0.0_wp
             interface_flux(4:n_vars) = 0.0_wp
          END IF

          CALL eval_hydrostatic_coefficient( recon%qp_interfaceB(:,j,k),     &
               reduced_gravity_left, gamma_left )
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceT(:,j,k),     &
               reduced_gravity_right, gamma_right )
          CALL regularize_dry_gamma_pair( recon%qp_interfaceB(1,j,k), gamma_left, &
               recon%qp_interfaceT(1,j,k), gamma_right )
          CALL eval_hydrostatic_path( PATH_DIR_Y,                            &
               recon%qp_interfaceB(1,j,k), gamma_left,                       &
               recon%eta_interfaceB(j,k), grav_coeff_stag_y(j,k),            &
               recon%qp_interfaceT(1,j,k), gamma_right,                      &
               recon%eta_interfaceT(j,k), grav_coeff_stag_y(j,k),            &
               path_contribution )

          CALL eval_oriented_pccu_pair( interface_flux,                       &
               path_contribution, a_minus, a_plus,                           &
               this%G_interface_yB(:,j,k), this%G_interface_yT(:,j,k) )

       END DO
       !$OMP END PARALLEL DO

    END IF

  END SUBROUTINE eval_flux_PCCU

  !******************************************************************************
  !> \brief Integrate the hydrostatic path between the final traces of each active cell.
  !>
  !> \param[in,out] this Persistent face-speed, oriented-flux and cell-path workspace.
  !> \param[in] recon Reconstructed cell-owned and face-oriented physical/conservative traces.
  !> \param[in] solve_cells Number of entries in the active-cell index lists.
  !> \param[in] j_cent X indices of active cells; only the first solve_cells entries are used.
  !> \param[in] k_cent Y indices of active cells; only the first solve_cells entries are used.
  !******************************************************************************

  SUBROUTINE eval_cell_hydrostatic_paths( this, recon, solve_cells, j_cent,   &
       k_cent )

    USE equation_terms_2d, ONLY : eval_hydrostatic_coefficient
    USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y
    USE pccu_2d, ONLY : eval_hydrostatic_path

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    INTEGER, INTENT(IN) :: solve_cells
    INTEGER, INTENT(IN) :: j_cent(:), k_cent(:)

    REAL(wp) :: gamma_minus, gamma_plus
    REAL(wp) :: reduced_gravity_minus, reduced_gravity_plus
    INTEGER :: j, k, l

    !$OMP PARALLEL DO private(l,j,k,gamma_minus,gamma_plus,                   &
    !$OMP & reduced_gravity_minus,reduced_gravity_plus)
    DO l = 1, solve_cells

       j = j_cent(l)
       k = k_cent(l)

       IF ( comp_cells_x .GT. 1 ) THEN
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceR(:,j,k),     &
               reduced_gravity_minus, gamma_minus )
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceL(:,j+1,k),   &
               reduced_gravity_plus, gamma_plus )
          CALL regularize_dry_gamma_pair( recon%qp_interfaceR(1,j,k), gamma_minus, &
               recon%qp_interfaceL(1,j+1,k), gamma_plus )
          CALL eval_hydrostatic_path( PATH_DIR_X,                            &
               recon%qp_interfaceR(1,j,k), gamma_minus,                      &
               recon%eta_cellW(j,k), grav_coeff_stag_x(j,k),                 &
               recon%qp_interfaceL(1,j+1,k), gamma_plus,                     &
               recon%eta_cellE(j,k), grav_coeff_stag_x(j+1,k),               &
               this%P_cell_x(:,j,k) )
       END IF

       IF ( comp_cells_y .GT. 1 ) THEN
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceT(:,j,k),     &
               reduced_gravity_minus, gamma_minus )
          CALL eval_hydrostatic_coefficient( recon%qp_interfaceB(:,j,k+1),   &
               reduced_gravity_plus, gamma_plus )
          CALL regularize_dry_gamma_pair( recon%qp_interfaceT(1,j,k), gamma_minus, &
               recon%qp_interfaceB(1,j,k+1), gamma_plus )
          CALL eval_hydrostatic_path( PATH_DIR_Y,                            &
               recon%qp_interfaceT(1,j,k), gamma_minus,                      &
               recon%eta_cellS(j,k), grav_coeff_stag_y(j,k),                 &
               recon%qp_interfaceB(1,j,k+1), gamma_plus,                     &
               recon%eta_cellN(j,k), grav_coeff_stag_y(j,k+1),               &
               this%P_cell_y(:,j,k) )
       END IF

    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE eval_cell_hydrostatic_paths

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
    IF ( ( h_left .LE. hp_dry_tolerance ) .AND.                              &
         ( h_right .GT. hp_dry_tolerance ) ) THEN
       gamma_left = gamma_right
    ELSEIF ( ( h_right .LE. hp_dry_tolerance ) .AND.                         &
         ( h_left .GT. hp_dry_tolerance ) ) THEN
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
    USE equation_terms_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y

    IMPLICIT NONE

    CLASS(hyperbolic_workspace_type), INTENT(INOUT) :: this
    CLASS(reconstruction_workspace_type), INTENT(IN) :: recon
    INTEGER, INTENT(IN) :: solve_interfaces_x, solve_interfaces_y
    INTEGER, INTENT(IN) :: j_stag_x(:), k_stag_x(:)
    INTEGER, INTENT(IN) :: j_stag_y(:), k_stag_y(:)

    REAL(wp) :: abslambdaL_min , abslambdaL_max
    REAL(wp) :: abslambdaR_min , abslambdaR_max
    REAL(wp) :: abslambdaB_min , abslambdaB_max
    REAL(wp) :: abslambdaT_min , abslambdaT_max
    REAL(wp) :: min_r , max_r

    INTEGER :: j,k,l

    !$OMP PARALLEL

    IF ( comp_cells_x .GT. 1 ) THEN

       !$OMP DO private(j , k , abslambdaL_min , abslambdaL_max ,               &
       !$OMP & abslambdaR_min , abslambdaR_max , min_r , max_r )

       x_interfaces_loop:DO l = 1,solve_interfaces_x

          j = j_stag_x(l)
          k = k_stag_x(l)

          CALL eval_local_speeds_x( recon%qp_interfaceL(:,j,k) ,             &
               grav_coeff_stag_x(j,k), abslambdaL_min , abslambdaL_max )

          CALL eval_local_speeds_x( recon%qp_interfaceR(:,j,k) ,             &
               grav_coeff_stag_x(j,k), abslambdaR_min , abslambdaR_max )

          min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0_wp)
          max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0_wp)

          this%a_interface_xNeg(j,k) = min_r
          this%a_interface_xPos(j,k) = max_r

       END DO x_interfaces_loop

       !$OMP END DO NOWAIT

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       !$OMP DO private(j , k , abslambdaB_min , abslambdaB_max ,               &
       !$OMP & abslambdaT_min , abslambdaT_max , min_r , max_r )

       y_interfaces_loop:DO l = 1,solve_interfaces_y

          j = j_stag_y(l)
          k = k_stag_y(l)

          CALL eval_local_speeds_y( recon%qp_interfaceB(:,j,k) ,             &
               grav_coeff_stag_y(j,k), abslambdaB_min , abslambdaB_max )
          
          CALL eval_local_speeds_y( recon%qp_interfaceT(:,j,k) ,             &
               grav_coeff_stag_y(j,k), abslambdaT_min , abslambdaT_max )

          min_r = MIN(abslambdaB_min , abslambdaT_min , 0.0_wp)
          max_r = MAX(abslambdaB_max , abslambdaT_max , 0.0_wp)

          this%a_interface_yNeg(j,k) = min_r
          this%a_interface_yPos(j,k) = max_r

       END DO y_interfaces_loop

       !$OMP END DO

    END IF

    !$OMP END PARALLEL

    RETURN
    
  END SUBROUTINE eval_speeds

END MODULE hyperbolic_2d
