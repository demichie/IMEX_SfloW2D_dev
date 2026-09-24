!********************************************************************************
!> \brief Immutable structural description of one configured model
!>
!> The layout snapshots variable counts and component indices before solver
!> workspaces are allocated. State arrays remain three-dimensional in this
!> phase; the explicit layer count prepares the later rank migration.
!********************************************************************************
MODULE model_layout_2d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: model_layout_type
  PUBLIC :: supported_layer_count

  TYPE :: model_layout_type

     INTEGER :: n_layers = 0
     INTEGER :: n_variables = 0
     INTEGER :: n_physical_variables = 0
     INTEGER :: n_equations = 0

     INTEGER :: n_solid = 0
     INTEGER :: n_additional_gas = 0
     INTEGER :: n_stochastic = 0
     INTEGER :: n_pore_pressure = 0
     LOGICAL :: has_liquid_fraction = .FALSE.

     INTEGER :: idx_mass = 0
     INTEGER :: idx_momentum_x = 0
     INTEGER :: idx_momentum_y = 0
     INTEGER :: idx_energy = 0
     INTEGER :: idx_solid_first = 0
     INTEGER :: idx_solid_last = 0
     INTEGER :: idx_additional_gas_first = 0
     INTEGER :: idx_additional_gas_last = 0
     INTEGER :: idx_stochastic = 0
     INTEGER :: idx_pore_pressure = 0
     INTEGER :: idx_liquid_fraction = 0
     INTEGER :: idx_velocity_x = 0
     INTEGER :: idx_velocity_y = 0

   CONTAINS

     PROCEDURE :: initialize => initialize_model_layout
     PROCEDURE :: finalize => finalize_model_layout
     PROCEDURE :: is_initialized => model_layout_is_initialized

  END TYPE model_layout_type

CONTAINS

  !******************************************************************************
  !> \brief Report whether a layer count is implemented by the current solver
  !******************************************************************************
  PURE LOGICAL FUNCTION supported_layer_count( n_layers )

    INTEGER, INTENT(IN) :: n_layers

    supported_layer_count = n_layers .EQ. 1

  END FUNCTION supported_layer_count


  !******************************************************************************
  !> \brief Build and validate the structural model layout
  !******************************************************************************
  SUBROUTINE initialize_model_layout( this, n_layers, n_variables,           &
       n_equations, n_solid, n_additional_gas, n_stochastic,                &
       n_pore_pressure, has_liquid_fraction )

    CLASS(model_layout_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: n_layers
    INTEGER, INTENT(IN) :: n_variables
    INTEGER, INTENT(IN) :: n_equations
    INTEGER, INTENT(IN) :: n_solid
    INTEGER, INTENT(IN) :: n_additional_gas
    INTEGER, INTENT(IN) :: n_stochastic
    INTEGER, INTENT(IN) :: n_pore_pressure
    LOGICAL, INTENT(IN) :: has_liquid_fraction

    INTEGER :: next_index

    CALL this%finalize

    IF ( .NOT. supported_layer_count(n_layers) ) THEN
       ERROR STOP 'Only N_LAYERS=1 is currently supported'
    END IF

    IF ( n_solid < 0 .OR. n_additional_gas < 0 ) THEN
       ERROR STOP 'Model layout component counts must be non-negative'
    END IF

    IF ( n_stochastic < 0 .OR. n_stochastic > 1 ) THEN
       ERROR STOP 'Model layout supports zero or one stochastic variable'
    END IF

    IF ( n_pore_pressure < 0 .OR. n_pore_pressure > 1 ) THEN
       ERROR STOP 'Model layout supports zero or one pore-pressure variable'
    END IF

    this%n_layers = n_layers
    this%n_variables = n_variables
    this%n_physical_variables = n_variables + 2
    this%n_equations = n_equations
    this%n_solid = n_solid
    this%n_additional_gas = n_additional_gas
    this%n_stochastic = n_stochastic
    this%n_pore_pressure = n_pore_pressure
    this%has_liquid_fraction = has_liquid_fraction

    this%idx_mass = 1
    this%idx_momentum_x = 2
    this%idx_momentum_y = 3
    this%idx_energy = 4

    this%idx_solid_first = 5
    this%idx_solid_last = 4 + n_solid
    this%idx_additional_gas_first = this%idx_solid_last + 1
    this%idx_additional_gas_last = this%idx_solid_last + n_additional_gas

    next_index = this%idx_additional_gas_last

    IF ( n_stochastic == 1 ) THEN
       next_index = next_index + 1
       this%idx_stochastic = next_index
    END IF

    IF ( n_pore_pressure == 1 ) THEN
       next_index = next_index + 1
       this%idx_pore_pressure = next_index
    END IF

    IF ( has_liquid_fraction ) THEN
       next_index = next_index + 1
       this%idx_liquid_fraction = next_index
    END IF

    this%idx_velocity_x = n_variables + 1
    this%idx_velocity_y = n_variables + 2

    IF ( next_index /= n_variables ) THEN
       ERROR STOP 'Configured variable count is inconsistent with model layout'
    END IF

    IF ( n_equations /= n_variables ) THEN
       ERROR STOP 'Current model requires one equation per conservative variable'
    END IF

  END SUBROUTINE initialize_model_layout


  !******************************************************************************
  !> \brief Reset the descriptor to its uninitialized state
  !******************************************************************************
  SUBROUTINE finalize_model_layout( this )

    CLASS(model_layout_type), INTENT(INOUT) :: this

    this%n_layers = 0
    this%n_variables = 0
    this%n_physical_variables = 0
    this%n_equations = 0
    this%n_solid = 0
    this%n_additional_gas = 0
    this%n_stochastic = 0
    this%n_pore_pressure = 0
    this%has_liquid_fraction = .FALSE.

    this%idx_mass = 0
    this%idx_momentum_x = 0
    this%idx_momentum_y = 0
    this%idx_energy = 0
    this%idx_solid_first = 0
    this%idx_solid_last = 0
    this%idx_additional_gas_first = 0
    this%idx_additional_gas_last = 0
    this%idx_stochastic = 0
    this%idx_pore_pressure = 0
    this%idx_liquid_fraction = 0
    this%idx_velocity_x = 0
    this%idx_velocity_y = 0

  END SUBROUTINE finalize_model_layout


  !******************************************************************************
  !> \brief Check that the descriptor contains a complete supported layout
  !******************************************************************************
  LOGICAL FUNCTION model_layout_is_initialized( this )

    CLASS(model_layout_type), INTENT(IN) :: this

    model_layout_is_initialized = supported_layer_count(this%n_layers)      &
         .AND. this%n_variables >= 4                                        &
         .AND. this%n_solid >= 0                                            &
         .AND. this%n_additional_gas >= 0                                   &
         .AND. this%n_stochastic >= 0 .AND. this%n_stochastic <= 1           &
         .AND. this%n_pore_pressure >= 0 .AND. this%n_pore_pressure <= 1     &
         .AND. this%n_variables == 4 + this%n_solid                          &
              + this%n_additional_gas + this%n_stochastic                   &
              + this%n_pore_pressure                                        &
              + MERGE(1, 0, this%has_liquid_fraction)                       &
         .AND. this%n_physical_variables == this%n_variables + 2            &
         .AND. this%n_equations == this%n_variables                         &
         .AND. this%idx_mass == 1                                           &
         .AND. this%idx_momentum_x == 2                                     &
         .AND. this%idx_momentum_y == 3                                     &
         .AND. this%idx_energy == 4                                         &
         .AND. this%idx_solid_first == 5                                    &
         .AND. this%idx_solid_last == 4 + this%n_solid                      &
         .AND. this%idx_additional_gas_first == 5 + this%n_solid            &
         .AND. this%idx_additional_gas_last == 4 + this%n_solid             &
              + this%n_additional_gas                                       &
         .AND. this%idx_stochastic == MERGE(                                &
              5 + this%n_solid + this%n_additional_gas, 0,                  &
              this%n_stochastic == 1)                                       &
         .AND. this%idx_pore_pressure == MERGE(                             &
              5 + this%n_solid + this%n_additional_gas                      &
              + this%n_stochastic, 0, this%n_pore_pressure == 1)            &
         .AND. this%idx_liquid_fraction == MERGE(this%n_variables, 0,       &
              this%has_liquid_fraction)                                     &
         .AND. this%idx_velocity_x == this%n_variables + 1                  &
         .AND. this%idx_velocity_y == this%n_variables + 2

  END FUNCTION model_layout_is_initialized

END MODULE model_layout_2d
