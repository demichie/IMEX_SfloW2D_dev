!********************************************************************************
!> \brief Model-facing explicit spatial operator
!>
!> This module is the boundary between time integration and the numerical
!> discretization of spatial terms. It owns the KT backend and
!> exposes only complete spatial-term and CFL evaluations.  Time integration
!> therefore does not depend on flux, path-contribution, or wave-speed storage.
!********************************************************************************

MODULE spatial_operator_2d
    USE parameters_2d, ONLY : idx_totMassEqn, idx_h, idx_T, idx_u, idx_v, dry_thickness_tolerance

    USE parameters_2d, ONLY : wp, n_eqns, n_vars, max_dt, cfl, cfl_policy
    USE diagnostics_2d, ONLY : fatal_error
    USE constitutive_parameters_2d, ONLY : T_ambient
    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, comp_interfaces_x, comp_interfaces_y, dx, dy
    USE state_conversions, ONLY : conservative_to_primitive
    USE domain_2d, ONLY : domain_type
    USE reconstruction_2d, ONLY : reconstruction_workspace_type
    USE hyperbolic_2d, ONLY : hyperbolic_workspace_type, add_interlayer_terms
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite

    IMPLICIT NONE

    PRIVATE

    !> \brief Numerical backend owned by one model-facing spatial operator.
    !> \details Its private reconstruction/flux workspaces are reused by both
    !>          spatial-term evaluation and the CFL calculation.
    TYPE, PUBLIC :: spatial_operator_type
    PRIVATE
        TYPE(reconstruction_workspace_type) :: reconstruction
        TYPE(hyperbolic_workspace_type) :: hp_pccu
        LOGICAL :: profile_enabled=.FALSE.
        INTEGER :: profile_cfl_calls=0
        REAL(wp) :: profile_cfl_seconds=0.0_wp
    CONTAINS
        PROCEDURE, PUBLIC :: initialize => initialize_spatial_operator
        PROCEDURE, PUBLIC :: finalize => finalize_spatial_operator
        PROCEDURE, PUBLIC :: evaluate => evaluate_spatial_operator
        PROCEDURE, PUBLIC :: evaluate_coupled => evaluate_coupled_spatial_operator
        PROCEDURE, PUBLIC :: diagnose_cell => diagnose_spatial_cell
        PROCEDURE, PUBLIC :: front_mass_flux => spatial_front_mass_flux
        PROCEDURE, PUBLIC :: compute_timestep => compute_spatial_timestep
        PROCEDURE, PUBLIC :: compute_coupled_timestep => compute_coupled_spatial_timestep
        PROCEDURE, PUBLIC :: prepare_hp_workset => prepare_spatial_hp_workset
    END TYPE spatial_operator_type

CONTAINS

    !> @brief Report oriented HP-PCCU fluxes of a failed spatial cell.
    !> @param[in] this Workspace holding the last evaluated spatial fluxes.
    !> @param[in] j Horizontal cell index.
    !> @param[in] k Vertical cell index.
    !> @param[in] eq Conservative equation index.
    SUBROUTINE diagnose_spatial_cell(this,j,k,eq)
        CLASS(spatial_operator_type), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: j,k,eq
        CALL this%hp_pccu%diagnose_cell(this%reconstruction,j,k,eq)
    END SUBROUTINE diagnose_spatial_cell

    !> @brief Read the mass flux facing each cell at a selected x-interface.
    !> @details Only used for opt-in diagnostics; never recomputes numerical fluxes.
    SUBROUTINE spatial_front_mass_flux(this,j,k,flux_left,flux_right)
        CLASS(spatial_operator_type), INTENT(IN) :: this
        INTEGER, INTENT(IN) :: j,k
        REAL(wp), INTENT(OUT) :: flux_left,flux_right

        IF (j<1 .OR. j>comp_interfaces_x .OR. k<1 .OR. k>comp_cells_y) &
            CALL fatal_error('Front diagnostic face outside computational domain')
        flux_left=this%hp_pccu%G_interface_xL(idx_totMassEqn,j,k)
        flux_right=this%hp_pccu%G_interface_xR(idx_totMassEqn,j,k)
        WRITE(*,'(A,4(ES24.16,1X))') 'FRONT_TRACE selected_hL,hR,etaL,etaR=', &
            this%reconstruction%h_selected_interfaceL(j,k), &
            this%reconstruction%h_selected_interfaceR(j,k), &
            this%reconstruction%eta_interfaceL(j,k),this%reconstruction%eta_interfaceR(j,k)
        WRITE(*,'(A,4(ES24.16,1X))') 'FRONT_TRACE reconstructed_hL,hR,uL,uR=', &
            this%reconstruction%qp_adv_interfaceL(idx_h,j,k), &
            this%reconstruction%qp_adv_interfaceR(idx_h,j,k), &
            this%reconstruction%qp_adv_interfaceL(idx_u,j,k), &
            this%reconstruction%qp_adv_interfaceR(idx_u,j,k)
        WRITE(*,'(A,4(ES24.16,1X))') 'FRONT_TRACE h_froudeL,h_froudeR,a_minus,a_plus=', &
            this%reconstruction%h_froude_interfaceL(j,k), &
            this%reconstruction%h_froude_interfaceR(j,k), &
            this%hp_pccu%a_interface_xNeg(j,k),this%hp_pccu%a_interface_xPos(j,k)
    END SUBROUTINE spatial_front_mass_flux

    !> \brief Initialize the reconstruction and KT flux workspace.
    !>
    !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.
    !> @param[in] layer Optional material layer slot (1 or 2); omitted means legacy L1 behavior.

    SUBROUTINE initialize_spatial_operator( this, layer )

        CLASS(spatial_operator_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN), OPTIONAL :: layer
        INTEGER :: selected_layer

        BLOCK
            CHARACTER(LEN=8) :: value
            INTEGER :: status
            value=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_PROFILE_LIE',value,STATUS=status)
            this%profile_enabled=status==0 .AND. TRIM(value)=='1'
        END BLOCK
        selected_layer=1
        IF (PRESENT(layer)) selected_layer=layer
        IF (selected_layer<1 .OR. selected_layer>2) CALL fatal_error('Invalid spatial layer slot')

        CALL this%reconstruction%initialize(layer=selected_layer)
        CALL this%hp_pccu%initialize(clamp_buoyant_pressure=selected_layer==2)

    END SUBROUTINE initialize_spatial_operator

    !******************************************************************************
    !> \brief Evaluate the complete explicit spatial term on the active domain.
    !>
    !> The returned term includes conservative interface transport and the
    !> cell/interface path contributions assembled by HP-PCCU.  Its sign follows
    !> q_t + spatial_term = local_sources.
    !>
    !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.
    !> \param[in] qp Canonical physical states: h, hu, hv, T, component mass fractions, then
    !> appended
    !>               u and v.
    !> \param[out] spatial_term Complete explicit spatial term with sign
    !>                          q_t+spatial_term=local_sources.
    !> \param[in] time Current simulation time [s].
    !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
    !> \param[in] boundary_source_enabled Optional override for boundary-source traces.
    !>
    !> \note The sign convention is q_t+spatial_term=local_sources; hydrostatic pressure/topography
    !>       are already included through path terms.
    !******************************************************************************

    SUBROUTINE evaluate_spatial_operator( this, qp, spatial_term, time, domain, boundary_source_enabled,               &
        reuse_cfl_traces, lower_operator, h_lower )

        CLASS(spatial_operator_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: qp(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(OUT) :: spatial_term(n_eqns,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: time
        CLASS(domain_type), INTENT(IN) :: domain
        LOGICAL, INTENT(IN), OPTIONAL :: boundary_source_enabled,reuse_cfl_traces
        CLASS(spatial_operator_type), INTENT(IN), OPTIONAL :: lower_operator
        REAL(wp), INTENT(IN), OPTIONAL :: h_lower(:,:)

        IF (PRESENT(lower_operator) .NEQV. PRESENT(h_lower)) &
            ERROR STOP 'Upper-layer spatial evaluation needs L1 operator and thickness'
        IF (PRESENT(lower_operator)) THEN
            IF (this%reconstruction%layer/=2 .OR. lower_operator%reconstruction%layer/=1) &
                ERROR STOP 'Invalid L1/L2 spatial operator pairing'
        END IF

        IF (PRESENT(lower_operator)) THEN
            CALL this%hp_pccu%evaluate_terms( this%reconstruction, qp, spatial_term, &
                time, domain%solve_cells, domain%j_cent, domain%k_cent, &
                domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x, &
                domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y, &
                source_boundary_enabled=boundary_source_enabled, reuse_cfl_traces=reuse_cfl_traces, &
                lower_recon=lower_operator%reconstruction,h_lower=h_lower)
        ELSE
            CALL this%hp_pccu%evaluate_terms( this%reconstruction, qp, spatial_term,                                   &
                time, domain%solve_cells, domain%j_cent, domain%k_cent,                                                &
                domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,                                           &
                domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y,                                           &
                source_boundary_enabled=boundary_source_enabled, reuse_cfl_traces=reuse_cfl_traces )
        END IF

    END SUBROUTINE evaluate_spatial_operator

    !> @brief Evaluate both HP-PCCU residuals with shared total-depth face bounds.
    !> @details Requires up-to-date L1/L2 reconstruction and hydrostatic
    !> coefficients from CFL passes on the SAME union workset. Reuses both
    !> layer operators and adds the interlayer terms once; no second engine.
    !> This remains test-only until coupled eigenvalues and shear are validated.
    !> @param[in,out] this Prepared lower-layer spatial operator.
    !> @param[in,out] upper Prepared upper-layer spatial operator.
    !> @param[in] qp1,qp2 Current layer primitive fields.
    !> @param[out] term1,term2 Coupled spatial residuals, on active cells.
    !> @param[in] time Explicit stage time [s].
    !> @param[in] domain Shared OR cell, face and HP reconstruction workset.
    !> @param[out] dt Preliminary shared CFL timestep [s].
    SUBROUTINE evaluate_coupled_spatial_operator(this,upper,qp1,qp2,term1,term2,time,domain,dt)
        CLASS(spatial_operator_type), INTENT(INOUT) :: this,upper
        REAL(wp), INTENT(IN) :: qp1(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: qp2(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(OUT) :: term1(n_eqns,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(OUT) :: term2(n_eqns,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: time
        CLASS(domain_type), INTENT(IN) :: domain
        REAL(wp), INTENT(OUT) :: dt
        INTEGER :: l,j,k

        IF (this%reconstruction%layer/=1 .OR. upper%reconstruction%layer/=2) &
            ERROR STOP 'Coupled spatial operator requires ordered L1/L2'

        ! Hydrostatic gamma coefficients were cached by each isolated
        ! reconstruction pass; overwrite only the characteristic bounds.
        CALL this%hp_pccu%evaluate_total_depth_speeds(this%reconstruction,upper%reconstruction, &
            domain%solve_interfaces_x,domain%j_stag_x,domain%k_stag_x, &
            domain%solve_interfaces_y,domain%j_stag_y,domain%k_stag_y)

        ! Both backends use the same coupled speeds only on active faces.
        ! Do not copy unrelated inactive grid entries at every timestep.
        ! Keep hydrostatic Gamma and reconstructed geometry layer-local.
        DO l=1,domain%solve_interfaces_x
            j=domain%j_stag_x(l)
            k=domain%k_stag_x(l)
            upper%hp_pccu%a_interface_xNeg(j,k)=this%hp_pccu%a_interface_xNeg(j,k)
            upper%hp_pccu%a_interface_xPos(j,k)=this%hp_pccu%a_interface_xPos(j,k)
        END DO
        DO l=1,domain%solve_interfaces_y
            j=domain%j_stag_y(l)
            k=domain%k_stag_y(l)
            upper%hp_pccu%a_interface_yNeg(j,k)=this%hp_pccu%a_interface_yNeg(j,k)
            upper%hp_pccu%a_interface_yPos(j,k)=this%hp_pccu%a_interface_yPos(j,k)
        END DO

        CALL this%compute_coupled_timestep(domain, &
            this%hp_pccu%a_interface_xNeg,this%hp_pccu%a_interface_xPos, &
            this%hp_pccu%a_interface_yNeg,this%hp_pccu%a_interface_yPos,dt)

        ! Use the two previously prepared reconstructions, without rebuilding
        ! the faces and undoing the common characteristic bounds.
        CALL this%evaluate(qp1,term1,time,domain,reuse_cfl_traces=.TRUE.)
        CALL upper%evaluate(qp2,term2,time,domain,reuse_cfl_traces=.TRUE., &
            lower_operator=this,h_lower=qp1(1,:,:))

        CALL add_interlayer_terms(this%reconstruction,upper%reconstruction,upper%hp_pccu, &
            term1,term2,domain%solve_mask,domain%solve_cells,domain%j_cent,domain%k_cent,qp2, &
            domain%solve_interfaces_x,domain%j_stag_x,domain%k_stag_x, &
            domain%solve_interfaces_y,domain%j_stag_y,domain%k_stag_y, &
            this%hp_pccu%a_interface_xNeg,this%hp_pccu%a_interface_xPos, &
            this%hp_pccu%a_interface_yNeg,this%hp_pccu%a_interface_yPos)
    END SUBROUTINE evaluate_coupled_spatial_operator

    !******************************************************************************
    !> \brief Refresh the reconstruction states and compute the face-speed CFL timestep.
    !>
    !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.
    !> \param[in] q Cell-centered conservative states, indexed as (variable,x-cell,y-cell).
    !> \param[in,out] qp Physical state cache refreshed from q on the reconstruction workset.
    !> \param[in] time Current simulation time [s].
    !> \param[out] dt Returned CFL-limited timestep, bounded by max_dt [s].
    !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
    !> \param[in] boundary_source_enabled Optional override for boundary-source traces during CFL
    !> evaluation.
    !> @param[out] max_speed_x Optional maximum x-face characteristic speed for diagnostics [m/s].
    !> @param[out] max_speed_y Optional maximum y-face characteristic speed for diagnostics [m/s].
    !> @param[in] diagnose_max_x_face Print the CFL-limiting x-face and its input states when requested.
    !> @param[in] lower_operator Optional prepared L1 spatial operator; L2 reconstruction only.
    !> @param[in] h_lower Optional lower-layer thickness, paired with lower_operator.
    !> @note The returned CFL still uses isolated-layer characteristic bounds;
    !>       it is not a validated coupled two-layer stability bound.
    !******************************************************************************

    SUBROUTINE compute_spatial_timestep( this, q, qp, time, dt, domain, boundary_source_enabled, &
        lower_operator, h_lower, max_speed_x, max_speed_y, diagnose_max_x_face )

        CLASS(spatial_operator_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: q(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(INOUT) :: qp(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: time
        REAL(wp), INTENT(OUT) :: dt
        CLASS(domain_type), INTENT(IN) :: domain
        LOGICAL, INTENT(IN), OPTIONAL :: boundary_source_enabled
        CLASS(spatial_operator_type), INTENT(IN), OPTIONAL :: lower_operator
        REAL(wp), INTENT(IN), OPTIONAL :: h_lower(:,:)
        REAL(wp), INTENT(OUT), OPTIONAL :: max_speed_x,max_speed_y
        LOGICAL, INTENT(IN), OPTIONAL :: diagnose_max_x_face

        INTEGER :: j, k, l, face_j, face_k
        REAL(wp) :: face_speed
        REAL(wp) :: qp_face_l(n_vars),qp_face_r(n_vars)
        REAL(wp) :: max_a_x, max_a_y
        REAL(wp) :: dynamic_pressure
        INTEGER :: tick,tock,tick_rate
        LOGICAL :: conversion_ok,conversion_failed
        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tick,tick_rate)
            this%profile_cfl_calls=this%profile_cfl_calls+1
        END IF

        dt = max_dt
        IF (PRESENT(max_speed_x)) max_speed_x=0.0_wp
        IF (PRESENT(max_speed_y)) max_speed_y=0.0_wp

        ! The sentinel selects the prescribed max_dt without evaluating face speeds.
        IF ( cfl .EQ. -1.0_wp ) THEN
            IF (this%profile_enabled) THEN
                CALL SYSTEM_CLOCK(tock)
                this%profile_cfl_seconds=this%profile_cfl_seconds+REAL(tock-tick,wp)/REAL(tick_rate,wp)
            END IF
            RETURN
        END IF

        ! Refresh wet cells and the dry halo: the HP stencil reads beyond solve_cells.
        ! The end-of-loop barrier makes every cache entry available to reconstruction.
        conversion_failed=.FALSE.
        !$OMP PARALLEL DO PRIVATE(j,k,dynamic_pressure,conversion_ok) REDUCTION(.OR.:conversion_failed)
        DO l = 1, domain%reconstruction_cells
            j = domain%j_reconstruction(l)
            k = domain%k_reconstruction(l)

            IF ( q(idx_totMassEqn,j,k) .GT. 0.0_wp ) THEN
                CALL conservative_to_primitive(q(:,j,k),qp(:,j,k),accepted=conversion_ok, &
                    p_dyn=dynamic_pressure)
                IF (.NOT.conversion_ok) conversion_failed=.TRUE.
            ELSE
                qp(:,j,k) = 0.0_wp
                qp(idx_T,j,k) = T_ambient

            END IF

        END DO
        !$OMP END PARALLEL DO
        IF (conversion_failed) CALL fatal_error('CFL reconstruction conversion rejected')

        ! Preserve the existing single-layer reconstruction path unchanged.
        ! Upper-layer traces use the already prepared L1 HP surfaces as local
        ! cell-side bases; this does not yet implement coupled eigenvalues.
        IF (PRESENT(lower_operator) .NEQV. PRESENT(h_lower)) &
            ERROR STOP 'Upper-layer CFL requires lower operator and thickness'
        IF (PRESENT(lower_operator)) THEN
            IF (this%reconstruction%layer/=2 .OR. lower_operator%reconstruction%layer/=1) &
                ERROR STOP 'Invalid layer pairing in upper-layer CFL'
            CALL this%reconstruction%reconstruct(qp,time,domain%solve_cells, &
                domain%j_cent,domain%k_cent,.TRUE.,source_boundary_enabled=boundary_source_enabled, &
                lower_recon=lower_operator%reconstruction,h_lower=h_lower)
            CALL this%reconstruction%reconstruct_kt_eta(qp,lower_recon=lower_operator%reconstruction)
        ELSE
            ! CFL speeds must use the same selected advective traces as the spatial operator.
            CALL this%reconstruction%reconstruct( qp, time, domain%solve_cells,                                        &
                domain%j_cent, domain%k_cent,                                                                          &
                .TRUE., source_boundary_enabled=boundary_source_enabled )

            ! Match the KT face depths and desingularized velocities used by the
            ! spatial operator before evaluating the CFL wave speeds.
            CALL this%reconstruction%reconstruct_kt_eta( qp )
        END IF

        CALL this%hp_pccu%evaluate_speeds( this%reconstruction,                                                        &
            domain%solve_interfaces_x, domain%j_stag_x, domain%k_stag_x,                                               &
            domain%solve_interfaces_y, domain%j_stag_y, domain%k_stag_y )

        max_a_x = 0.0_wp
        max_a_y = 0.0_wp

        ! Each active cell contributes all of its faces; a global maximum gives
        ! one timestep valid for every participating OpenMP thread.
        !$OMP PARALLEL DO PRIVATE(j,k) REDUCTION(MAX:max_a_x,max_a_y)
        DO l = 1, domain%solve_cells
            j = domain%j_cent(l)
            k = domain%k_cent(l)

            CALL accumulate_directional_face_speeds(this%hp_pccu,j,k, max_a_x,max_a_y)

        END DO
        !$OMP END PARALLEL DO

        CALL limit_timestep_by_spatial_speeds(max_a_x,max_a_y,dt)
        IF (PRESENT(max_speed_x)) max_speed_x=max_a_x
        IF (PRESENT(max_speed_y)) max_speed_y=max_a_y
        ! Run this serial diagnostic only after all OpenMP speed reductions
        ! have completed. Reuse the exact selected traces and cached bounds.
        IF (PRESENT(diagnose_max_x_face)) THEN
            IF (diagnose_max_x_face .AND. max_a_x>100.0_wp) THEN
                face_j=0
                face_k=0
                DO l=1,domain%solve_interfaces_x
                    j=domain%j_stag_x(l)
                    k=domain%k_stag_x(l)
                    face_speed=MAX(this%hp_pccu%a_interface_xPos(j,k), &
                        -this%hp_pccu%a_interface_xNeg(j,k))
                    IF (face_speed==max_a_x) THEN
                        face_j=j
                        face_k=k
                        EXIT
                    END IF
                END DO
                IF (face_j>0) THEN
                    j=face_j
                    k=face_k
                    qp_face_l=this%reconstruction%qp_adv_interfaceL(:,j,k)
                    qp_face_r=this%reconstruction%qp_adv_interfaceR(:,j,k)
                    qp_face_l(1)=this%reconstruction%h_selected_interfaceL(j,k)
                    qp_face_r(1)=this%reconstruction%h_selected_interfaceR(j,k)
                    IF (qp_face_l(1)<=dry_thickness_tolerance) qp_face_l(idx_u:idx_v)=0.0_wp
                    IF (qp_face_r(1)<=dry_thickness_tolerance) qp_face_r(idx_u:idx_v)=0.0_wp
                    WRITE(*,'(A,2(A,I0),3(A,ES24.16))') 'CFL_X_FACE', &
                        ' j=',j,' k=',k, &
                        ' a_neg=',this%hp_pccu%a_interface_xNeg(j,k), &
                        ' a_pos=',this%hp_pccu%a_interface_xPos(j,k), &
                        ' max_a_x=',max_a_x
                    WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE left_speed_input=',qp_face_l
                    WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE right_speed_input=',qp_face_r
                    WRITE(*,'(A,4(A,ES24.16))') 'CFL_X_FACE geometry', &
                        ' selected_h_l=',qp_face_l(1),' selected_h_r=',qp_face_r(1), &
                        ' gamma_l=',this%hp_pccu%gamma_interface_xL(j,k), &
                        ' gamma_r=',this%hp_pccu%gamma_interface_xR(j,k)
                    IF (j>1) THEN
                        WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE left_center_q=',q(:,j-1,k)
                        WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE left_center_qp=',qp(:,j-1,k)
                    END IF
                    IF (j<=comp_cells_x) THEN
                        WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE right_center_q=',q(:,j,k)
                        WRITE(*,'(A,*(ES24.16,1X))') 'CFL_X_FACE right_center_qp=',qp(:,j,k)
                    END IF
                END IF
            END IF
        END IF

        IF (this%profile_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            this%profile_cfl_seconds=this%profile_cfl_seconds+REAL(tock-tick,wp)/REAL(tick_rate,wp)
        END IF
    END SUBROUTINE compute_spatial_timestep

    !> @brief Limit a timestep using externally validated coupled face-speed bounds.
    !> @details The caller owns the coupled eigensystem, hyperbolicity check
    !> and bounds. The method only reduces the same directional face rates
    !> and applies the unchanged single-layer CFL policy. In particular,
    !> isolated L1/L2 eigenvalues are not substituted for coupled bounds.
    !> @param[in] this Spatial operator carrying the public timestep API.
    !> @param[in] domain Shared OR cell/face workset for the coupled pair.
    !> @param[in] axneg,axpos Coupled nonpositive/nonnegative x-face bounds.
    !> @param[in] ayneg,aypos Coupled nonpositive/nonnegative y-face bounds.
    !> @param[out] dt Shared CFL-limited timestep [s].
    SUBROUTINE compute_coupled_spatial_timestep(this,domain,axneg,axpos,ayneg,aypos,dt)
        CLASS(spatial_operator_type), INTENT(IN) :: this
        CLASS(domain_type), INTENT(IN) :: domain
        REAL(wp), INTENT(IN) :: axneg(comp_interfaces_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: axpos(comp_interfaces_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: ayneg(comp_cells_x,comp_interfaces_y)
        REAL(wp), INTENT(IN) :: aypos(comp_cells_x,comp_interfaces_y)
        REAL(wp), INTENT(OUT) :: dt
        REAL(wp) :: max_a_x,max_a_y
        INTEGER :: l,j,k

        dt=max_dt
        IF (cfl==-1.0_wp) RETURN

        ! Validate each actual face once; inactive array entries are irrelevant.
        DO l=1,domain%solve_interfaces_x
            j=domain%j_stag_x(l)
            k=domain%k_stag_x(l)
            IF (.NOT.ieee_is_finite(axneg(j,k)) .OR. .NOT.ieee_is_finite(axpos(j,k))) &
                CALL fatal_error('Nonfinite coupled x-face wave bound')
            IF (axneg(j,k)>0.0_wp .OR. axpos(j,k)<0.0_wp) &
                CALL fatal_error('Coupled x-face wave bounds have invalid signs')
        END DO

        DO l=1,domain%solve_interfaces_y
            j=domain%j_stag_y(l)
            k=domain%k_stag_y(l)
            IF (.NOT.ieee_is_finite(ayneg(j,k)) .OR. .NOT.ieee_is_finite(aypos(j,k))) &
                CALL fatal_error('Nonfinite coupled y-face wave bound')
            IF (ayneg(j,k)>0.0_wp .OR. aypos(j,k)<0.0_wp) &
                CALL fatal_error('Coupled y-face wave bounds have invalid signs')
        END DO

        max_a_x=0.0_wp
        max_a_y=0.0_wp
        DO l=1,domain%solve_cells
            j=domain%j_cent(l)
            k=domain%k_cent(l)
            ! Match the four-face rate and MAX operand order of the L1 CFL.
            IF (comp_cells_x>1) THEN
                max_a_x=MAX(max_a_x,axpos(j,k),-axneg(j,k), &
                    axpos(j+1,k),-axneg(j+1,k))
            END IF
            IF (comp_cells_y>1) THEN
                max_a_y=MAX(max_a_y,aypos(j,k),-ayneg(j,k), &
                    aypos(j,k+1),-ayneg(j,k+1))
            END IF
        END DO

        CALL limit_timestep_by_spatial_speeds(max_a_x,max_a_y,dt)
    END SUBROUTINE compute_coupled_spatial_timestep

    !> @brief Apply the canonical CFL policy to maximum directional wave speeds.
    !> @details Keeps the legacy floating-point evaluation order and 1D branch.
    !> @param[in] max_a_x,max_a_y Maximum face speed magnitudes [m/s].
    !> @param[in,out] dt Existing maximum permitted timestep [s].
    SUBROUTINE limit_timestep_by_spatial_speeds(max_a_x,max_a_y,dt)
        REAL(wp), INTENT(IN) :: max_a_x,max_a_y
        REAL(wp), INTENT(INOUT) :: dt
        REAL(wp) :: cfl_spatial,propagation_rate

        cfl_spatial = cfl
        cfl_spatial = MIN(cfl_spatial,0.25_wp)

        SELECT CASE (TRIM(cfl_policy))
        CASE ('directional_min')
            IF ( max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dx/max_a_x)
            IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dy/max_a_y)
        CASE ('combined')
            ! Use the legacy division order in 1D: mathematical equivalence alone
            ! would allow a one-ulp difference from cfl/(lambda/dx).
            IF ( comp_cells_y .EQ. 1 ) THEN
                IF ( comp_cells_x .GT. 1 .AND. max_a_x .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dx/max_a_x)
            ELSEIF ( comp_cells_x .EQ. 1 ) THEN
                IF ( max_a_y .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial*dy/max_a_y)
            ELSE
                propagation_rate = max_a_x/dx + max_a_y/dy
                IF ( propagation_rate .GT. 0.0_wp ) dt = MIN(dt,cfl_spatial/propagation_rate)

            END IF
        CASE DEFAULT
            ! Also protect callers that configure the module without input I/O.
            CALL fatal_error('CFL_POLICY must be directional_min or combined')

        END SELECT
    END SUBROUTINE limit_timestep_by_spatial_speeds

    !> \brief Include all four reconstructed face speeds for one active cell.
    !> \details The MAX operand order is kept identical to the original CFL loop.
    !> \param[in] speeds Reconstructed face-speed arrays.
    !> \param[in] j X index of the active cell.
    !> \param[in] k Y index of the active cell.
    !> \param[in,out] max_a_x Maximum outward wave speed in x.
    !> \param[in,out] max_a_y Maximum outward wave speed in y.
    SUBROUTINE accumulate_directional_face_speeds(speeds,j,k,max_a_x,max_a_y)
        CLASS(hyperbolic_workspace_type), INTENT(IN) :: speeds
        INTEGER, INTENT(IN) :: j,k
        REAL(wp), INTENT(INOUT) :: max_a_x,max_a_y

        max_a_x = MAX(max_a_x,                                                                                         &
            speeds%a_interface_xPos(j,k),                                                                              &
            -speeds%a_interface_xNeg(j,k),                                                                             &
            speeds%a_interface_xPos(j+1,k),                                                                            &
            -speeds%a_interface_xNeg(j+1,k))

        max_a_y = MAX(max_a_y,                                                                                         &
            speeds%a_interface_yPos(j,k),                                                                              &
            -speeds%a_interface_yNeg(j,k),                                                                             &
            speeds%a_interface_yPos(j,k+1),                                                                            &
            -speeds%a_interface_yNeg(j,k+1))
    END SUBROUTINE accumulate_directional_face_speeds

    !> \brief Release the reconstruction and KT flux workspace.
    !>
    !> \param[in,out] this Spatial operator owning reconstruction and HP-PCCU backend storage.

    SUBROUTINE finalize_spatial_operator( this )

        CLASS(spatial_operator_type), INTENT(INOUT) :: this

        IF (this%profile_enabled) THEN
            WRITE(*,'(A,I0,A,ES15.7)') 'LIE_PROFILE cfl_calls=',this%profile_cfl_calls, &
                ' seconds=',this%profile_cfl_seconds
            WRITE(*,'(A,I0)') 'LIE_PROFILE spatial_calls=',this%hp_pccu%profile_calls
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE reconstruction_seconds=',this%hp_pccu%profile_reconstruction
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE flux_and_paths_seconds=',this%hp_pccu%profile_flux
            WRITE(*,'(A,ES15.7)') 'LIE_PROFILE assembly_seconds=',this%hp_pccu%profile_assembly
        END IF
        CALL this%hp_pccu%finalize
        CALL this%reconstruction%finalize

    END SUBROUTINE finalize_spatial_operator

    !> @brief Reuse the HP workset prepared by the active-domain refresh.
    !> @param[in,out] this Spatial operator and owned reconstruction workset.
    !> @param[in] domain Active-domain worklist supplying HP workset indices.
    SUBROUTINE prepare_spatial_hp_workset(this, domain)
        CLASS(spatial_operator_type), INTENT(INOUT) :: this
        CLASS(domain_type), INTENT(IN) :: domain
        CALL this%reconstruction%prepare_hp_workset(                                                                   &
            domain%hp_eta_cells, domain%hp_eta_j, domain%hp_eta_k)
    END SUBROUTINE prepare_spatial_hp_workset

END MODULE spatial_operator_2d
