!********************************************************************************
!> @brief Describe immutable structural component indices for one configured model.
!>
!> The layout snapshots variable counts and component indices before solver
!> workspaces are allocated. State arrays are rank four (component,layer,x,y). The production driver
!> controls experimental two-layer execution at the input boundary.
!********************************************************************************

MODULE model_layout_2d

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: model_layout_type
    PUBLIC :: supported_layer_count

    !> @brief Validated component indices shared by states and numerical workspaces.
    !> @details Zero optional indices mean that the corresponding transported
    !>          variable is absent. Appended physical velocities are not equations.
    TYPE :: model_layout_type

        !> Number of material layers; one or two structural layer slots.
        INTEGER :: n_layers = 0
        !> Number of conservative state components per cell.
        INTEGER :: n_vars = 0
        !> Number of conservative balance equations.
        INTEGER :: n_eqns = 0

        INTEGER :: n_solid = 0
        INTEGER :: n_add_gas = 0
        INTEGER :: n_stoch_vars = 0
        INTEGER :: n_pore_vars = 0
        INTEGER :: idx_totMassEqn = 0
        INTEGER :: idx_uEqn = 0
        INTEGER :: idx_vEqn = 0
        INTEGER :: idx_engyEqn = 0
        !> Temperature component of the primitive state, not a conservative energy equation.
        INTEGER :: idx_T = 0
        INTEGER :: idx_solid_first = 0
        INTEGER :: idx_solid_last = 0
        INTEGER :: idx_add_gas_first = 0
        INTEGER :: idx_add_gas_last = 0
        INTEGER :: idx_stoch = 0
        INTEGER :: idx_pore = 0
        INTEGER :: idx_u = 0
        INTEGER :: idx_v = 0

    CONTAINS

        PROCEDURE :: initialize => initialize_model_layout
        PROCEDURE :: finalize => finalize_model_layout
        PROCEDURE :: is_initialized => model_layout_is_initialized

    END TYPE model_layout_type

CONTAINS

    !******************************************************************************
    !> @brief Report whether the requested layer count is implemented.
    !>
    !> @param[in] n_layers Requested material-layer count; one or two layers.
    !> @return True for the structurally implemented one- or two-layer layouts.
    !******************************************************************************

    PURE LOGICAL FUNCTION supported_layer_count( n_layers )

        INTEGER, INTENT(IN) :: n_layers

        supported_layer_count = n_layers >= 1 .AND. n_layers <= 2

    END FUNCTION supported_layer_count


    !******************************************************************************
    !> @brief Validate variable counts and build the canonical component index layout.
    !>
    !> @param[in,out] this Canonical model variable-count and component-index descriptor.
    !> @param[in] n_layers Requested material-layer count; one or two layers.
    !> @param[in] n_vars Number of conservative variables in each cell state.
    !> @param[in] n_eqns Number of balance equations; currently equal to n_vars.
    !> @param[in] n_solid Number of independently transported solid classes.
    !> @param[in] n_add_gas Number of additional gas components, excluding ambient air.
    !> @param[in] n_stoch_vars Number of transported stochastic variables, zero or one.
    !> @param[in] n_pore_vars Number of transported pore-pressure variables, zero or one.
    !******************************************************************************

    SUBROUTINE initialize_model_layout(this, n_layers, n_vars, n_eqns, n_solid, n_add_gas,            &
        n_stoch_vars, n_pore_vars)

        CLASS(model_layout_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: n_layers
        INTEGER, INTENT(IN) :: n_vars
        INTEGER, INTENT(IN) :: n_eqns
        INTEGER, INTENT(IN) :: n_solid
        INTEGER, INTENT(IN) :: n_add_gas
        INTEGER, INTENT(IN) :: n_stoch_vars
        INTEGER, INTENT(IN) :: n_pore_vars
        INTEGER :: next_index

        CALL this%finalize

        IF ( .NOT. supported_layer_count(n_layers) ) THEN
            ERROR STOP 'Model layout requires N_LAYERS=1 or N_LAYERS=2'

        END IF

        IF ( n_solid < 0 .OR. n_add_gas < 0 ) THEN
            ERROR STOP 'Model layout component counts must be non-negative'

        END IF

        IF ( n_stoch_vars < 0 .OR. n_stoch_vars > 1 ) THEN
            ERROR STOP 'Model layout supports zero or one stochastic variable'

        END IF

        IF ( n_pore_vars < 0 .OR. n_pore_vars > 1 ) THEN
            ERROR STOP 'Model layout supports zero or one pore-pressure variable'

        END IF

        this%n_layers = n_layers
        this%n_vars = n_vars
        this%n_eqns = n_eqns
        this%n_solid = n_solid
        this%n_add_gas = n_add_gas
        this%n_stoch_vars = n_stoch_vars
        this%n_pore_vars = n_pore_vars
        this%idx_totMassEqn = 1
        this%idx_uEqn = 2
        this%idx_vEqn = 3
        this%idx_engyEqn = 4
        this%idx_T = 4

        this%idx_solid_first = 5
        this%idx_solid_last = 4 + n_solid
        this%idx_add_gas_first = this%idx_solid_last + 1
        this%idx_add_gas_last = this%idx_solid_last + n_add_gas

        next_index = this%idx_add_gas_last

        IF ( n_stoch_vars == 1 ) THEN
            next_index = next_index + 1
            this%idx_stoch = next_index

        END IF

        IF ( n_pore_vars == 1 ) THEN
            next_index = next_index + 1
            this%idx_pore = next_index

        END IF

        this%idx_u = 2
        this%idx_v = 3

        IF ( next_index /= n_vars ) THEN
            ERROR STOP 'Configured variable count is inconsistent with model layout'

        END IF

        IF ( n_eqns /= n_vars ) THEN
            ERROR STOP 'Current model requires one equation per conservative variable'

        END IF

    END SUBROUTINE initialize_model_layout


    !******************************************************************************
    !> @brief Reset all counts and indices to the uninitialized state.
    !>
    !> @param[in,out] this Canonical model variable-count and component-index descriptor.
    !******************************************************************************

    SUBROUTINE finalize_model_layout( this )

        CLASS(model_layout_type), INTENT(INOUT) :: this

        this%n_layers = 0
        this%n_vars = 0
        this%n_eqns = 0
        this%n_solid = 0
        this%n_add_gas = 0
        this%n_stoch_vars = 0
        this%n_pore_vars = 0
        this%idx_totMassEqn = 0
        this%idx_uEqn = 0
        this%idx_vEqn = 0
        this%idx_engyEqn = 0
        this%idx_T = 0
        this%idx_solid_first = 0
        this%idx_solid_last = 0
        this%idx_add_gas_first = 0
        this%idx_add_gas_last = 0
        this%idx_stoch = 0
        this%idx_pore = 0
        this%idx_u = 0
        this%idx_v = 0

    END SUBROUTINE finalize_model_layout


    !******************************************************************************
    !> @brief Check that all layout counts and indices are mutually consistent.
    !>
    !> @param[in] this Canonical model variable-count and component-index descriptor.
    !> @return True when the descriptor contains a complete supported layout.
    !******************************************************************************

    LOGICAL FUNCTION model_layout_is_initialized( this )

        CLASS(model_layout_type), INTENT(IN) :: this

        model_layout_is_initialized = supported_layer_count(this%n_layers)                                             &
            .AND. this%n_vars >= 4                                                                                &
            .AND. this%n_solid >= 0                                                                                    &
            .AND. this%n_add_gas >= 0                                                                           &
            .AND. this%n_stoch_vars >= 0 .AND. this%n_stoch_vars <= 1                                                  &
            .AND. this%n_pore_vars >= 0 .AND. this%n_pore_vars <= 1                                            &
            .AND. this%n_vars == 4 + this%n_solid                                                                 &
                    + this%n_add_gas + this%n_stoch_vars                                                        &
                    + this%n_pore_vars                                                                             &
            .AND. this%n_eqns == this%n_vars                                                                 &
            .AND. this%idx_totMassEqn == 1                                                                             &
            .AND. this%idx_uEqn == 2                                                                             &
            .AND. this%idx_vEqn == 3                                                                             &
            .AND. this%idx_T == 4                                                                                     &
            .AND. this%idx_engyEqn == 4                                                                                &
            .AND. this%idx_solid_first == 5                                                                            &
            .AND. this%idx_solid_last == 4 + this%n_solid                                                              &
            .AND. this%idx_add_gas_first == 5 + this%n_solid                                                    &
            .AND. this%idx_add_gas_last == 4 + this%n_solid                                                     &
                    + this%n_add_gas                                                                            &
            .AND. this%idx_stoch == MERGE(                                                                        &
                    5 + this%n_solid + this%n_add_gas, 0,                                                       &
                    this%n_stoch_vars == 1)                                                                            &
            .AND. this%idx_pore == MERGE(                                                                     &
                    5 + this%n_solid + this%n_add_gas                                                           &
                    + this%n_stoch_vars, 0, this%n_pore_vars == 1)                                                 &
            .AND. this%idx_u == 2                                                          &
            .AND. this%idx_v == 3

    END FUNCTION model_layout_is_initialized

END MODULE model_layout_2d
