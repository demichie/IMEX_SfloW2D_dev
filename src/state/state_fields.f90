!********************************************************************************
!> \brief Solver state ownership
!>
!> This module owns the conservative and physical state arrays together with
!> diagnostics that evolve during a simulation.
!>
!> Stores conservative/physical cell states and accumulated hazard/runout diagnostics. Configuration
!> and numerical scratch storage belong to other modules.
!********************************************************************************

MODULE state_fields
    USE, INTRINSIC :: iso_fortran_env, ONLY : int64
    USE parameters_2d, ONLY : idx_h, idx_T,                                                                            &
            idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last

    USE constitutive_parameters_2d, ONLY : T_ambient
    USE substrate_parameters_2d, ONLY : substrate

    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE model_layout_2d, ONLY : model_layout_type

    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE parameters_2d, ONLY : n_eqns,idx_add_gas_first,idx_add_gas_last
    USE parameters_2d, ONLY : wp,n_solid,lava_flag
    USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: initialize_dry_physical_cell
    PUBLIC :: physical_cell_state_finite, finite_real_value
    !> \brief Owned cell states and diagnostics for one simulation.
    !> \details Fluid arrays use (component,layer,x-cell,y-cell); diagnostic maps use
    !>          (x-cell,y-cell). The model layout fixes component counts and indices.
    TYPE, PUBLIC :: state_type

        !> Conservative variables: M=rho_m*h, M*u, M*v, M*cp*T, then transported masses.
        REAL(wp), ALLOCATABLE :: q(:,:,:,:)

        !> Primitive variables: h, u, v, T, then transported mass fractions and scalars.
        REAL(wp), ALLOCATABLE :: qp(:,:,:,:)

        !> Physical substrate inventories owned by the same simulation state as q and qp.
        REAL(wp), ALLOCATABLE :: deposit(:,:,:),deposit_tot(:,:)
        !> Cumulative eroded solid inventory by (x,y,class), plus total map (x,y).
        REAL(wp), ALLOCATABLE :: erosion(:,:,:),erosion_tot(:,:)
        !> Remaining erodible substrate inventory by (class,x,y).
        REAL(wp), ALLOCATABLE :: erodible(:,:,:)

        !> Map of positive thickness
        LOGICAL, ALLOCATABLE :: hpos(:,:)

        !> Map of positive thickness at previous output step
        LOGICAL, ALLOCATABLE :: hpos_old(:,:)

        !> Maximum over time of thickness
        REAL(wp), ALLOCATABLE :: hmax(:,:)

        !> Maximum over time of dynamic pressure
        REAL(wp), ALLOCATABLE :: pdynmax(:,:)

        !> Maximum over time of velocity magnitude
        REAL(wp), ALLOCATABLE :: mod_vel_max(:,:)

        !> Vulnerability diagnostic table indexed by threshold and horizontal cell.
        LOGICAL, ALLOCATABLE :: vuln_table(:,:,:)

    CONTAINS

        PROCEDURE :: initialize => initialize_state
        PROCEDURE :: finalize => finalize_state
        PROCEDURE :: refresh_substrate_totals => refresh_substrate_totals

    END TYPE state_type

CONTAINS

    !> @brief Allocate and initialize all simulation-owned state fields.
    !> @param[inout] this Authoritative simulation state owning grid arrays and substrate
    !> inventories.
    !> @param[in] model_layout Validated runtime counts and canonical variable indices.
    SUBROUTINE initialize_state(this, model_layout)

        CLASS(state_type), INTENT(INOUT) :: this
        TYPE(model_layout_type), INTENT(IN) :: model_layout

        IF (.NOT. model_layout%is_initialized()) THEN
            ERROR STOP 'State allocation requires an initialized model layout'

        END IF

        ALLOCATE( this%q(model_layout%n_vars,model_layout%n_layers,comp_cells_x,comp_cells_y) )

        ALLOCATE( this%hpos(comp_cells_x,comp_cells_y),                                                                &
            this%hpos_old(comp_cells_x,comp_cells_y) )

        ALLOCATE( this%qp(model_layout%n_vars,model_layout%n_layers,comp_cells_x,                     &
            comp_cells_y) )

        ALLOCATE(this%deposit(comp_cells_x,comp_cells_y,n_solid))
        ALLOCATE(this%erosion(comp_cells_x,comp_cells_y,n_solid))
        ALLOCATE(this%erodible(n_solid,comp_cells_x,comp_cells_y))
        ALLOCATE(this%deposit_tot(comp_cells_x,comp_cells_y))
        ALLOCATE(this%erosion_tot(comp_cells_x,comp_cells_y))
        this%deposit=0.0_wp
        this%erosion=0.0_wp
        this%erodible=0.0_wp
        this%deposit_tot=0.0_wp
        this%erosion_tot=0.0_wp
        this%q = 0.0_wp
        this%qp = 0.0_wp
        this%qp(model_layout%idx_T,1:model_layout%n_layers,1:comp_cells_x,1:comp_cells_y) = T_ambient

        ALLOCATE( this%hmax(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%pdynmax(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%mod_vel_max(comp_cells_x,comp_cells_y) )

        this%hpos = .FALSE.
        this%hpos_old = .FALSE.
        this%hmax = 0.0_wp
        this%pdynmax = 0.0_wp
        this%mod_vel_max = 0.0_wp

        ALLOCATE( this%vuln_table(n_thickness_levels*n_dyn_pres_levels,                                                &
            comp_cells_x,comp_cells_y) )


        this%vuln_table = .FALSE.

    END SUBROUTINE initialize_state

    !> @brief Deallocate simulation-owned state arrays and diagnostics.
    !> @param[inout] this Authoritative simulation state owning grid arrays and substrate
    !> inventories.
    SUBROUTINE finalize_state(this)

        CLASS(state_type), INTENT(INOUT) :: this

        DEALLOCATE( this%q, this%hpos, this%hpos_old )

        DEALLOCATE( this%hmax, this%pdynmax, this%mod_vel_max )

        DEALLOCATE( this%vuln_table )


        DEALLOCATE(this%deposit,this%deposit_tot,this%erosion,this%erosion_tot,                                        &
            this%erodible)
        DEALLOCATE( this%qp )

    END SUBROUTINE finalize_state

    !> @brief Initialize a primitive cell vector to the canonical dry state.
    !> @param[out] qp_cell Primitive cell vector with all variables for one cell.
    PURE SUBROUTINE initialize_dry_physical_cell(qp_cell)
        REAL(wp), INTENT(OUT) :: qp_cell(:)
        qp_cell = 0.0_wp
        qp_cell(idx_T) = T_ambient
    END SUBROUTINE initialize_dry_physical_cell

    !> @brief Refresh substrate totals without replacing the authoritative lava deposit.
    !> @details For lava deposit_tot is accumulated immobilized thickness; for
    !>          granular flows it is derived from class inventories and porosity.
    !> @param[inout] this Authoritative simulation state owning grid arrays and substrate
    !> inventories.
    !> @param[in] active_j Optional x indices of cells changed by the current Lie step.
    !> @param[in] active_k Optional y indices, paired with active_j.
    SUBROUTINE refresh_substrate_totals(this,active_j,active_k)
        CLASS(state_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN), OPTIONAL :: active_j(:),active_k(:)
        INTEGER :: j,k,cell,ncells
        ! Lava deposit_tot is authoritative; never reset it or use granular porosity.
        IF (lava_flag) THEN
            this%erosion_tot=0.0_wp
            RETURN
        END IF
        ! No granular solids: there is no solid bed layer to convert into thickness.
        IF (n_solid==0) THEN
            this%deposit_tot=0.0_wp
            this%erosion_tot=0.0_wp
            RETURN
        END IF
        ! Totals include the carrier-filled pore volume of the deposited/eroded layer.
        ! A zero-solid-volume substrate (porosity=1) cannot define this conversion.
        IF (.NOT.ieee_is_finite(substrate%erodible_porosity) .OR. &
            substrate%erodible_porosity<0.0_wp .OR. substrate%erodible_porosity>=1.0_wp) &
            ERROR STOP 'Invalid substrate porosity for aggregate thickness'
        ! Initial/restart refreshes cover the full grid; a Lie step only changes active cells.
        IF (PRESENT(active_j) .NEQV. PRESENT(active_k)) &
            ERROR STOP 'Substrate refresh requires paired active indices'
        ncells=comp_cells_x*comp_cells_y
        IF (PRESENT(active_j)) THEN
            IF (SIZE(active_j)/=SIZE(active_k)) ERROR STOP 'Mismatched substrate active indices'
            ncells=SIZE(active_j)
        END IF
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(cell,j,k) SCHEDULE(STATIC)
        DO cell=1,ncells
            IF (PRESENT(active_j)) THEN
                j=active_j(cell)
                k=active_k(cell)
            ELSE
                j=MOD(cell-1,comp_cells_x)+1
                k=(cell-1)/comp_cells_x+1
            END IF
            this%deposit_tot(j,k)=SUM(this%deposit(j,k,:))/(1.0_wp-substrate%erodible_porosity)
            this%erosion_tot(j,k)=SUM(this%erosion(j,k,:))/(1.0_wp-substrate%erodible_porosity)
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE refresh_substrate_totals

    !> @brief Test finiteness of all primitive components in one cell.
    !> @param[in] qp_cell Primitive cell vector with all variables for one cell.
    !> @return True if the requested state or composition check succeeds.
    PURE LOGICAL FUNCTION physical_cell_state_finite(qp_cell) RESULT(valid)
        REAL(wp), INTENT(IN) :: qp_cell(:)
        valid=ALL(finite_real_value(qp_cell))
    END FUNCTION physical_cell_state_finite

    !> @brief Detect finite binary64 values without altering IEEE runtime settings.
    !> @param[in] value One REAL scalar to test for finite binary64 representation.
    !> @return True if the requested state or composition check succeeds.
    PURE ELEMENTAL LOGICAL FUNCTION finite_real_value(value) RESULT(valid)
        REAL(wp), INTENT(IN) :: value
        INTEGER(int64), PARAMETER :: exponent_mask=INT(Z'7FF0000000000000',int64)
        INTEGER(int64) :: bits
        bits=TRANSFER(value,0_int64)
        valid=IAND(bits,exponent_mask)/=exponent_mask
    END FUNCTION finite_real_value

END MODULE state_fields
