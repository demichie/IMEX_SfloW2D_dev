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
  PUBLIC :: source_bundle_terms
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
  INTEGER, PARAMETER :: SOURCE_BUNDLE_LAHAR_DRAG = 1
  INTEGER, PARAMETER :: SOURCE_BUNDLE_THERMAL_COOLING = 2

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
    INTEGER :: bundle_id = 0 !< Nonzero for terms that must share one evaluation group.
  END TYPE source_term_type

  !> Single registration point for term capabilities and rheology membership.
  TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT) = [ &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB, &
            name='SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_voellmy, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_voellmy, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_plastic, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS, &
            name='SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lava, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAVA_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_LAVA_YIELD', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lava, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_LAHAR_DRAG), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_LAHAR_DRAG), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_YIELD', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_empirical, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_quadratic, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_froude, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lucas, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_MUI_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui_voellmy, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui_voellmy, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_THERMAL_RADIATIVE_COOLING, &
            name='SOURCE_TERM_THERMAL_RADIATIVE_COOLING', role=SOURCE_ROLE_THERMAL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_TEMPERATURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_THERMAL_COOLING), &
       source_term_type(id=SOURCE_TERM_THERMAL_CONVECTIVE_COOLING, &
            name='SOURCE_TERM_THERMAL_CONVECTIVE_COOLING', role=SOURCE_ROLE_THERMAL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_TEMPERATURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_THERMAL_COOLING), &
       source_term_type(id=SOURCE_TERM_PORE_PRESSURE_DISSIPATION, &
            name='SOURCE_TERM_PORE_PRESSURE_DISSIPATION', role=SOURCE_ROLE_PORE, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_PORE_PRESSURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION, &
            name='SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_BOTTOM_FISSURAL_INJECTION, &
            name='SOURCE_TERM_BOTTOM_FISSURAL_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RADIAL_INJECTION, &
            name='SOURCE_TERM_RADIAL_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_GEOMETRY_CURVATURE, &
            name='SOURCE_TERM_GEOMETRY_CURVATURE', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_BOUNDARY_INJECTION, &
            name='SOURCE_TERM_BOUNDARY_INJECTION', role=SOURCE_ROLE_BOUNDARY, &
            representation=SOURCE_REPRESENTATION_BOUNDARY_STATE, real_rate=.FALSE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_EROSION, &
            name='SOURCE_TERM_MASS_EXCHANGE_EROSION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_DEPOSITION, &
            name='SOURCE_TERM_MASS_EXCHANGE_DEPOSITION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT, &
            name='SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION, &
            name='SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE.), &
       source_term_type(id=SOURCE_TERM_INLET_THERMAL_COUPLING, &
            name='SOURCE_TERM_INLET_THERMAL_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE.), &
       source_term_type(id=SOURCE_TERM_INLET_PORE_COUPLING, &
            name='SOURCE_TERM_INLET_PORE_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_PORE_COUPLING, &
            name='SOURCE_TERM_MASS_PORE_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_STOCHASTIC_COUPLING, &
            name='SOURCE_TERM_MASS_STOCHASTIC_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE.) &
       ]

CONTAINS

  !> \brief Look up immutable catalog metadata without accessing runtime state.
  !> \param[in] term_id Stable source ID.
  !> \return Descriptor; id=0 for an unknown identifier.
  PURE FUNCTION source_term_descriptor(term_id) RESULT(term)
    INTEGER, INTENT(IN) :: term_id
    TYPE(source_term_type) :: term
    IF (term_id<1 .OR. term_id>SOURCE_TERM_COUNT) RETURN
    term=SOURCE_CATALOG(term_id)
  END FUNCTION source_term_descriptor

  !> \brief Return the registered terms that must share an evaluation group.
  !> \param[in] bundle_id Catalog batch ID; zero returns an empty list.
  !> \return Stable term IDs in catalog order.
  FUNCTION source_bundle_terms(bundle_id) RESULT(ids)
    INTEGER, INTENT(IN) :: bundle_id
    INTEGER, ALLOCATABLE :: ids(:)
    INTEGER :: i
    ALLOCATE(ids(0))
    IF (bundle_id<=0) RETURN
    DO i=1,SOURCE_TERM_COUNT
      IF (SOURCE_CATALOG(i)%bundle_id==bundle_id) ids=[ids,i]
    END DO
  END FUNCTION source_bundle_terms

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
