!********************************************************************************
!> \brief Stable source identifiers and interface capabilities for roadmap D17.
!> \details Contains metadata only: no source evaluation, clock or physical state.
!>          IDs follow the committed D17 inventory; dependent assemblies are not
!>          independently schedulable. Availability does not certify partial bundles.
!********************************************************************************
MODULE source_catalog_2d
  USE parameters_2d, ONLY : rheology_id_voellmy, &
       rheology_id_plastic, &
       rheology_id_lava, &
       rheology_id_lahar, &
       rheology_id_empirical, &
       rheology_id_quadratic, &
       rheology_id_froude, &
       rheology_id_lucas, &
       rheology_id_mui, &
       rheology_id_mui_voellmy
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_term_type, source_term_descriptor, rheology_source_terms
  PUBLIC :: SOURCE_TERM_COUNT, SOURCE_PHYSICAL_TERM_COUNT
  !> method explicit.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_EXPLICIT = 1
  !> method newton.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_NEWTON = 2
  !> method analytic.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_ANALYTIC = 3
  !> cadence stage.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_CADENCE_STAGE = 1
  !> cadence step.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_CADENCE_STEP = 2
  !> target velocity.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TARGET_VELOCITY = 1
  !> target temperature.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TARGET_TEMPERATURE = 2
  !> target pore pressure.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TARGET_PORE_PRESSURE = 3
  !> representation cell rate.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_REPRESENTATION_CELL_RATE = 1
  !> representation boundary state.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_REPRESENTATION_BOUNDARY_STATE = 2
  !> representation dependent.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_REPRESENTATION_DEPENDENT = 3
  !> role implicit.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_IMPLICIT = 1
  !> role semi implicit.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_SEMI_IMPLICIT = 2
  !> role thermal.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_THERMAL = 3
  !> role pore.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_PORE = 4
  !> role cell.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_CELL = 5
  !> role boundary.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_BOUNDARY = 6
  !> role mass.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_MASS = 7
  !> role derived.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_DERIVED = 8
  !> solve newton block.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_NEWTON_BLOCK = 1
  !> solve analytic velocity.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_VELOCITY = 2
  !> solve analytic temperature.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_TEMPERATURE = 3
  !> solve analytic pore pressure.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_PORE_PRESSURE = 4
  !> phase rk implicit.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_PHASE_RK_IMPLICIT = 1
  !> phase rk explicit.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_PHASE_RK_EXPLICIT = 2
  !> phase post rk.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_PHASE_POST_RK = 3
  !> rheology voellmy coulomb.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB = 1
  !> rheology voellmy turbulent drag.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG = 2
  !> rheology plastic yield.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD = 3
  !> rheology lava viscous.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS = 4
  !> rheology lava yield.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LAVA_YIELD = 5
  !> rheology lahar viscous.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS = 6
  !> rheology lahar turbulent drag.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG = 7
  !> rheology lahar yield.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LAHAR_YIELD = 8
  !> rheology empirical drag.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG = 9
  !> rheology quadratic drag.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG = 10
  !> rheology froude friction.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION = 11
  !> rheology lucas friction.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION = 12
  !> rheology mui friction.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_MUI_FRICTION = 13
  !> rheology mui voellmy friction.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION = 14
  !> rheology mui voellmy turbulent drag.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG = 15
  !> thermal radiative cooling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_THERMAL_RADIATIVE_COOLING = 16
  !> thermal convective cooling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_THERMAL_CONVECTIVE_COOLING = 17
  !> pore pressure dissipation.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_PORE_PRESSURE_DISSIPATION = 18
  !> bottom circular injection.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION = 19
  !> bottom fissural injection.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_BOTTOM_FISSURAL_INJECTION = 20
  !> radial injection.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_RADIAL_INJECTION = 21
  !> geometry curvature.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_GEOMETRY_CURVATURE = 22
  !> boundary injection.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_BOUNDARY_INJECTION = 23
  !> mass exchange erosion.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_EROSION = 24
  !> mass exchange deposition.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_DEPOSITION = 25
  !> mass exchange carrier liquid loss.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS = 26
  !> mass exchange carrier deposition loss.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS = 27
  !> mass exchange carrier darcy loss.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS = 28
  !> mass exchange entrainment.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT = 29
  !> mass exchange vaporization.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION = 30
  !> inlet thermal coupling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_INLET_THERMAL_COUPLING = 31
  !> inlet pore coupling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_INLET_PORE_COUPLING = 32
  !> mass pore coupling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_PORE_COUPLING = 33
  !> mass stochastic coupling.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TERM_MASS_STOCHASTIC_COUPLING = 34
  !> Number of catalogued physical and dependent contributions.
  INTEGER, PARAMETER :: SOURCE_TERM_COUNT = 34
  !> Physical terms eligible for scheduling (subject to capability/batch checks).
  INTEGER, PARAMETER :: SOURCE_PHYSICAL_TERM_COUNT = 30

  !> \brief Value descriptor of one physical or dependent source contribution.
  TYPE :: source_term_type
    INTEGER :: id = 0 !< Stable ID; zero denotes an invalid lookup.
    CHARACTER(LEN=64) :: name = "" !< Symbolic identifier.
    INTEGER :: role = 0 !< Rate/assembly role, not a separate temporal method.
    INTEGER :: representation = 0 !< Cell rate, boundary state or dependent assembly.
    INTEGER :: target = 0 !< Implicit target block, or zero for explicit/dependent work.
    INTEGER :: model_id = 0 !< Rheology family ID, zero for non-rheology entries.
    LOGICAL :: real_rate = .FALSE. !< Existing real rate evaluation.
    LOGICAL :: complex_rate = .FALSE. !< Existing complex rate evaluation.
    LOGICAL :: schedulable = .FALSE. !< False for dependent equation assembly.
  END TYPE source_term_type

CONTAINS

  !> \brief Look up immutable catalog metadata without accessing runtime state.
  !> \param[in] term_id Stable source ID.
  !> \return Descriptor; id=0 for an unknown identifier.
  PURE FUNCTION source_term_descriptor(term_id) RESULT(term)
    INTEGER, INTENT(IN) :: term_id
    TYPE(source_term_type) :: term
    SELECT CASE(term_id)
    CASE(SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB)
      term%name = 'SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_voellmy
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG)
      term%name = 'SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_voellmy
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD)
      term%name = 'SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_plastic
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS)
      term%name = 'SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lava
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LAVA_YIELD)
      term%name = 'SOURCE_TERM_RHEOLOGY_LAVA_YIELD'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lava
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS)
      term%name = 'SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lahar
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG)
      term%name = 'SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lahar
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LAHAR_YIELD)
      term%name = 'SOURCE_TERM_RHEOLOGY_LAHAR_YIELD'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lahar
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG)
      term%name = 'SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_empirical
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG)
      term%name = 'SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_quadratic
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION)
      term%name = 'SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_froude
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION)
      term%name = 'SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_lucas
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_MUI_FRICTION)
      term%name = 'SOURCE_TERM_RHEOLOGY_MUI_FRICTION'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_mui
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION)
      term%name = 'SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION'
      term%role = SOURCE_ROLE_SEMI_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_mui_voellmy
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG)
      term%name = 'SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG'
      term%role = SOURCE_ROLE_IMPLICIT
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_VELOCITY
      term%model_id = rheology_id_mui_voellmy
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_THERMAL_RADIATIVE_COOLING)
      term%name = 'SOURCE_TERM_THERMAL_RADIATIVE_COOLING'
      term%role = SOURCE_ROLE_THERMAL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_TEMPERATURE
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_THERMAL_CONVECTIVE_COOLING)
      term%name = 'SOURCE_TERM_THERMAL_CONVECTIVE_COOLING'
      term%role = SOURCE_ROLE_THERMAL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_TEMPERATURE
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_PORE_PRESSURE_DISSIPATION)
      term%name = 'SOURCE_TERM_PORE_PRESSURE_DISSIPATION'
      term%role = SOURCE_ROLE_PORE
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%target = SOURCE_TARGET_PORE_PRESSURE
      term%real_rate = .TRUE.
      term%complex_rate = .TRUE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION)
      term%name = 'SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION'
      term%role = SOURCE_ROLE_CELL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_BOTTOM_FISSURAL_INJECTION)
      term%name = 'SOURCE_TERM_BOTTOM_FISSURAL_INJECTION'
      term%role = SOURCE_ROLE_CELL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_RADIAL_INJECTION)
      term%name = 'SOURCE_TERM_RADIAL_INJECTION'
      term%role = SOURCE_ROLE_CELL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_GEOMETRY_CURVATURE)
      term%name = 'SOURCE_TERM_GEOMETRY_CURVATURE'
      term%role = SOURCE_ROLE_CELL
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_BOUNDARY_INJECTION)
      term%name = 'SOURCE_TERM_BOUNDARY_INJECTION'
      term%role = SOURCE_ROLE_BOUNDARY
      term%representation = SOURCE_REPRESENTATION_BOUNDARY_STATE
      term%real_rate = .FALSE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_EROSION)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_EROSION'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_DEPOSITION)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_DEPOSITION'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION)
      term%name = 'SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION'
      term%role = SOURCE_ROLE_MASS
      term%representation = SOURCE_REPRESENTATION_CELL_RATE
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .TRUE.
    CASE(SOURCE_TERM_INLET_THERMAL_COUPLING)
      term%name = 'SOURCE_TERM_INLET_THERMAL_COUPLING'
      term%role = SOURCE_ROLE_DERIVED
      term%representation = SOURCE_REPRESENTATION_DEPENDENT
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .FALSE.
    CASE(SOURCE_TERM_INLET_PORE_COUPLING)
      term%name = 'SOURCE_TERM_INLET_PORE_COUPLING'
      term%role = SOURCE_ROLE_DERIVED
      term%representation = SOURCE_REPRESENTATION_DEPENDENT
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .FALSE.
    CASE(SOURCE_TERM_MASS_PORE_COUPLING)
      term%name = 'SOURCE_TERM_MASS_PORE_COUPLING'
      term%role = SOURCE_ROLE_DERIVED
      term%representation = SOURCE_REPRESENTATION_DEPENDENT
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .FALSE.
    CASE(SOURCE_TERM_MASS_STOCHASTIC_COUPLING)
      term%name = 'SOURCE_TERM_MASS_STOCHASTIC_COUPLING'
      term%role = SOURCE_ROLE_DERIVED
      term%representation = SOURCE_REPRESENTATION_DEPENDENT
      term%real_rate = .TRUE.
      term%complex_rate = .FALSE.
      term%schedulable = .FALSE.
    CASE DEFAULT
      RETURN
    END SELECT
    term%id = term_id
  END FUNCTION source_term_descriptor

  !> \brief Return all rate IDs belonging to one selected rheology.
  !> \param[in] model_id Resolved rheology family ID.
  !> \return Ordered family bundle; an unknown/disabled family returns an empty list.
  FUNCTION rheology_source_terms(model_id) RESULT(ids)
    INTEGER, INTENT(IN) :: model_id
    INTEGER, ALLOCATABLE :: ids(:)
    TYPE(source_term_type) :: term
    INTEGER :: i
    ALLOCATE(ids(0))
    DO i=1,SOURCE_TERM_COUNT
      term=source_term_descriptor(i)
      IF (term%model_id==model_id .AND. model_id/=0) ids=[ids,i]
    END DO
  END FUNCTION rheology_source_terms

END MODULE source_catalog_2d
