!********************************************************************************
!> \brief Solver state ownership
!>
!> This module owns the conservative and physical state arrays together with
!> diagnostics that evolve during a simulation.
!>
!> Stores conservative/physical cell states and accumulated hazard/runout diagnostics. Configuration
!> and numerical scratch storage belong to other modules.
!********************************************************************************

MODULE state_2d
    USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T,                                                            &
            idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last

    USE constitutive_parameters_2d, ONLY : T_ambient,erodible_porosity

    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE model_layout_2d, ONLY : model_layout_type

    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE parameters_2d, ONLY : n_eqns,idx_add_gas_first,idx_add_gas_last
    USE parameters_2d, ONLY : wp,n_solid
    USE parameters_2d, ONLY : n_thickness_levels, n_dyn_pres_levels

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: initialize_dry_physical_cell
    PUBLIC :: conservative_mass_admissible
    PUBLIC :: physical_cell_state_finite
    PUBLIC :: substrate_snapshot_type

    !> Pre-step substrate inventory snapshot for transactional time integration.
    TYPE :: substrate_snapshot_type
        REAL(wp), ALLOCATABLE :: deposit(:,:,:),erosion(:,:,:),erodible(:,:,:)
        REAL(wp), ALLOCATABLE :: lava_deposit(:,:)
        REAL(wp), ALLOCATABLE :: deposit_tot(:,:),erosion_tot(:,:)
    END TYPE substrate_snapshot_type

    !> \brief Owned cell states and diagnostics for one simulation.
    !> \details State arrays use (component,x-cell,y-cell); diagnostic maps use
    !>          (x-cell,y-cell). The model layout fixes component counts and indices.
    TYPE, PUBLIC :: state_type

        !> Conservative variables: M=rho_m*h, M*u, M*v, M*cp*T, then transported masses.
        REAL(wp), ALLOCATABLE :: q(:,:,:)

        !> Physical variables: h, h*u, h*v, T, mass fractions, and appended u and v.
        REAL(wp), ALLOCATABLE :: qp(:,:,:)

        !> Physical substrate inventories owned by the same simulation state as q and qp.
        REAL(wp), ALLOCATABLE :: deposit(:,:,:),deposit_tot(:,:)
        REAL(wp), ALLOCATABLE :: erosion(:,:,:),erosion_tot(:,:)
        REAL(wp), ALLOCATABLE :: erodible(:,:,:),lava_deposit(:,:)

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
        PROCEDURE :: snapshot_substrate => snapshot_substrate
        PROCEDURE :: restore_substrate => restore_substrate

    END TYPE state_type

CONTAINS


    !> @brief Check that all primitive cell components are finite.
    !> @param[in] qp_cell One cell's physical variables, including velocities.
    !> @return True if every entry is finite; this procedure does not change state.
    PURE LOGICAL FUNCTION physical_cell_state_finite(qp_cell) RESULT(valid)
        REAL(wp), INTENT(IN) :: qp_cell(:)
        valid=ALL(ieee_is_finite(qp_cell))
    END FUNCTION physical_cell_state_finite


    !******************************************************************************
    !> \brief Allocate and reset cell states and accumulated diagnostics for one model layout.
    !>
    !> \param[in,out] this Allocated conservative/physical cell state and accumulated diagnostics.
    !> \param[in] model_layout Validated variable counts and canonical component indices.
    !******************************************************************************

    SUBROUTINE initialize_state(this, model_layout)

        CLASS(state_type), INTENT(INOUT) :: this
        TYPE(model_layout_type), INTENT(IN) :: model_layout

        IF (.NOT. model_layout%is_initialized()) THEN
            ERROR STOP 'State allocation requires an initialized model layout'

        END IF

        ALLOCATE( this%q(model_layout%n_variables,comp_cells_x,comp_cells_y) )

        ALLOCATE( this%hpos(comp_cells_x,comp_cells_y),                                                                &
            this%hpos_old(comp_cells_x,comp_cells_y) )

        ALLOCATE( this%qp(model_layout%n_physical_variables,comp_cells_x,                                              &
            comp_cells_y) )

        ALLOCATE(this%deposit(comp_cells_x,comp_cells_y,n_solid))
        ALLOCATE(this%erosion(comp_cells_x,comp_cells_y,n_solid))
        ALLOCATE(this%erodible(n_solid,comp_cells_x,comp_cells_y))
        ALLOCATE(this%lava_deposit(comp_cells_x,comp_cells_y))
        ALLOCATE(this%deposit_tot(comp_cells_x,comp_cells_y))
        ALLOCATE(this%erosion_tot(comp_cells_x,comp_cells_y))
        this%deposit=0.0_wp
        this%erosion=0.0_wp
        this%erodible=0.0_wp
        this%lava_deposit=0.0_wp
        this%deposit_tot=0.0_wp
        this%erosion_tot=0.0_wp
        this%q = 0.0_wp
        this%qp = 0.0_wp
        this%qp(model_layout%idx_energy,1:comp_cells_x,1:comp_cells_y) = T_ambient

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

    !> @brief Initialize one primitive cell, including halo, to the dry state.
    !> @details Sets all components to zero except ambient temperature.
    !> @param[out] qp_cell Local physical state vector to reset.
    PURE SUBROUTINE initialize_dry_physical_cell(qp_cell)
        REAL(wp), INTENT(OUT) :: qp_cell(:)
        qp_cell = 0.0_wp
        qp_cell(idx_T) = T_ambient
    END SUBROUTINE initialize_dry_physical_cell

    !******************************************************************************
    !> \brief Release the physical/conservative state and diagnostic fields.
    !>
    !> \param[in,out] this Allocated conservative/physical cell state and accumulated diagnostics.
    !******************************************************************************

    SUBROUTINE finalize_state(this)

        CLASS(state_type), INTENT(INOUT) :: this

        DEALLOCATE( this%q, this%hpos, this%hpos_old )

        DEALLOCATE( this%hmax, this%pdynmax, this%mod_vel_max )

        DEALLOCATE( this%vuln_table )


        DEALLOCATE(this%deposit,this%deposit_tot,this%erosion,this%erosion_tot,                                        &
            this%erodible,this%lava_deposit)
        DEALLOCATE( this%qp )

    END SUBROUTINE finalize_state

    !> @brief Test finite, nonnegative and internally consistent mass inventories.
    !> @details Checks total mass against summed solid and additional gas mass;
    !>          does not clip or otherwise modify a candidate conservative state.
    !> @param[in] q_cell One conservative cell vector of n_eqns entries.
    !> @return True if the candidate satisfies all mass admissibility guards.
    PURE LOGICAL FUNCTION conservative_mass_admissible(q_cell) RESULT(valid)
        REAL(wp), INTENT(IN) :: q_cell(n_eqns)
        valid=.FALSE.
        IF (.NOT.ALL(ieee_is_finite(q_cell))) RETURN
        IF (q_cell(idx_totMassEqn)<0.0_wp) RETURN
        IF (ANY(q_cell(idx_solidEqn_first:idx_solidEqn_last)<0.0_wp)) RETURN
        IF (ANY(q_cell(idx_add_gas_first:idx_add_gas_last)<0.0_wp)) RETURN
        ! All transported component masses share the same total-mixture budget.
        IF (SUM(q_cell(idx_solidEqn_first:idx_solidEqn_last)) +                                                        &
            SUM(q_cell(idx_add_gas_first:idx_add_gas_last)) > q_cell(idx_totMassEqn)) RETURN
        valid=.TRUE.
    END FUNCTION conservative_mass_admissible

    !> @brief Snapshot all owned substrate inventories and their aggregates.
    !> @param[in] this Authoritative state before the Lie source sweep.
    !> @param[out] snapshot Independent substrate snapshot for rollback.
    SUBROUTINE snapshot_substrate(this,snapshot)
        CLASS(state_type), INTENT(IN) :: this
        TYPE(substrate_snapshot_type), INTENT(OUT) :: snapshot
        snapshot%deposit=this%deposit
        snapshot%erosion=this%erosion
        snapshot%erodible=this%erodible
        snapshot%lava_deposit=this%lava_deposit
        snapshot%deposit_tot=this%deposit_tot
        snapshot%erosion_tot=this%erosion_tot
    END SUBROUTINE snapshot_substrate

    !> @brief Restore substrate inventories after a rejected Lie source sweep.
    !> @param[in,out] this State whose substrate fields are restored.
    !> @param[in] snapshot Previously captured substrate snapshot.
    SUBROUTINE restore_substrate(this,snapshot)
        CLASS(state_type), INTENT(INOUT) :: this
        TYPE(substrate_snapshot_type), INTENT(IN) :: snapshot
        this%deposit=snapshot%deposit
        this%erosion=snapshot%erosion
        this%erodible=snapshot%erodible
        this%lava_deposit=snapshot%lava_deposit
        this%deposit_tot=snapshot%deposit_tot
        this%erosion_tot=snapshot%erosion_tot
    END SUBROUTINE restore_substrate

    !> @brief Recompute deposit and erosion thickness totals from per-solid mass.
    !> @details Rejects invalid substrate porosity instead of silently clipping.
    !> @param[in,out] this State whose aggregate diagnostic fields are refreshed.
    SUBROUTINE refresh_substrate_totals(this)
        CLASS(state_type), INTENT(INOUT) :: this
        ! Totals include the carrier-filled pore volume of the deposited/eroded layer.
        ! A zero-solid-volume substrate (porosity=1) cannot define this conversion.
        IF (.NOT.ieee_is_finite(erodible_porosity) .OR.                                                                &
            erodible_porosity<0.0_wp .OR. erodible_porosity>=1.0_wp)                                                   &
            ERROR STOP 'Invalid substrate porosity for aggregate thickness'
        this%deposit_tot=SUM(this%deposit,DIM=3)/(1.0_wp-erodible_porosity)
        this%erosion_tot=SUM(this%erosion,DIM=3)/(1.0_wp-erodible_porosity)
    END SUBROUTINE refresh_substrate_totals

    !******************************************************************************
END MODULE state_2d
