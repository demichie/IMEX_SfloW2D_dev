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
  PUBLIC :: source_term_type, source_term_descriptor
  PUBLIC :: source_bundle_names
  PUBLIC :: source_term_descriptor_by_name, source_catalog_names_unique
  PUBLIC :: rheology_source_name, canonical_rheology_catalog_unique
  PUBLIC :: canonical_source_name, rheology_has_role
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
    LOGICAL :: has_explicit = .FALSE. !< Explicit rate/update interface is available.
    LOGICAL :: has_newton = .FALSE. !< Both Newton-facing interfaces exist.
    LOGICAL :: has_analytic = .FALSE. !< Analytic local update is available.
    INTEGER :: bundle_id = 0 !< Nonzero for terms that must share one evaluation group.
  END TYPE source_term_type

  !> One public identity per complete rheology; the legacy terms remain internal.
  TYPE :: rheology_law_type
    CHARACTER(LEN=64) :: name=''
    INTEGER :: model_id=0
  END TYPE rheology_law_type
  TYPE(rheology_law_type), PARAMETER :: RHEOLOGY_LAWS(10) = [ &
       rheology_law_type('rheology_voellmy',rheology_id_voellmy), &
       rheology_law_type('rheology_plastic',rheology_id_plastic), &
       rheology_law_type('rheology_lava',rheology_id_lava), &
       rheology_law_type('rheology_lahar',rheology_id_lahar), &
       rheology_law_type('rheology_empirical',rheology_id_empirical), &
       rheology_law_type('rheology_quadratic',rheology_id_quadratic), &
       rheology_law_type('rheology_froude',rheology_id_froude), &
       rheology_law_type('rheology_lucas',rheology_id_lucas), &
       rheology_law_type('rheology_mui',rheology_id_mui), &
       rheology_law_type('rheology_mui_voellmy',rheology_id_mui_voellmy) &
       ]

  !> Single registration point for term capabilities and rheology membership.
  TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT) = [ &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB, &
            name='SOURCE_TERM_RHEOLOGY_VOELLMY_COULOMB', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_voellmy, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_VOELLMY_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_voellmy, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_PLASTIC_YIELD', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_plastic, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS, &
            name='SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lava, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAVA_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_LAVA_YIELD', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lava, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_VISCOUS', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_LAHAR_DRAG, &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_LAHAR_DRAG, &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LAHAR_YIELD, &
            name='SOURCE_TERM_RHEOLOGY_LAHAR_YIELD', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lahar, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_EMPIRICAL_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_empirical, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_QUADRATIC_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_quadratic, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_FROUDE_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_froude, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_LUCAS_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_lucas, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_MUI_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION, &
            name='SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_FRICTION', role=SOURCE_ROLE_SEMI_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui_voellmy, real_rate=.TRUE., complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG, &
            name='SOURCE_TERM_RHEOLOGY_MUI_VOELLMY_TURBULENT_DRAG', role=SOURCE_ROLE_IMPLICIT, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_VELOCITY, &
            model_id=rheology_id_mui_voellmy, real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_THERMAL_RADIATIVE_COOLING, &
            name='SOURCE_TERM_THERMAL_RADIATIVE_COOLING', role=SOURCE_ROLE_THERMAL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_TEMPERATURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_THERMAL_COOLING, &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_THERMAL_CONVECTIVE_COOLING, &
            name='SOURCE_TERM_THERMAL_CONVECTIVE_COOLING', role=SOURCE_ROLE_THERMAL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_TEMPERATURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            bundle_id=SOURCE_BUNDLE_THERMAL_COOLING, &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_PORE_PRESSURE_DISSIPATION, &
            name='SOURCE_TERM_PORE_PRESSURE_DISSIPATION', role=SOURCE_ROLE_PORE, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, target=SOURCE_TARGET_PORE_PRESSURE, &
            real_rate=.TRUE., complex_rate=.TRUE., schedulable=.TRUE., &
            has_explicit=.FALSE., has_newton=.TRUE., &
            has_analytic=.TRUE.), &
       source_term_type(id=SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION, &
            name='SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_BOTTOM_FISSURAL_INJECTION, &
            name='SOURCE_TERM_BOTTOM_FISSURAL_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_RADIAL_INJECTION, &
            name='SOURCE_TERM_RADIAL_INJECTION', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_GEOMETRY_CURVATURE, &
            name='SOURCE_TERM_GEOMETRY_CURVATURE', role=SOURCE_ROLE_CELL, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_BOUNDARY_INJECTION, &
            name='SOURCE_TERM_BOUNDARY_INJECTION', role=SOURCE_ROLE_BOUNDARY, &
            representation=SOURCE_REPRESENTATION_BOUNDARY_STATE, real_rate=.FALSE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_EROSION, &
            name='SOURCE_TERM_MASS_EXCHANGE_EROSION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_DEPOSITION, &
            name='SOURCE_TERM_MASS_EXCHANGE_DEPOSITION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS, &
            name='SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT, &
            name='SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION, &
            name='SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION', role=SOURCE_ROLE_MASS, &
            representation=SOURCE_REPRESENTATION_CELL_RATE, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.TRUE., &
            has_explicit=.TRUE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_INLET_THERMAL_COUPLING, &
            name='SOURCE_TERM_INLET_THERMAL_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE., &
            has_explicit=.FALSE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_INLET_PORE_COUPLING, &
            name='SOURCE_TERM_INLET_PORE_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE., &
            has_explicit=.FALSE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_PORE_COUPLING, &
            name='SOURCE_TERM_MASS_PORE_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE., &
            has_explicit=.FALSE., has_newton=.FALSE., &
            has_analytic=.FALSE.), &
       source_term_type(id=SOURCE_TERM_MASS_STOCHASTIC_COUPLING, &
            name='SOURCE_TERM_MASS_STOCHASTIC_COUPLING', role=SOURCE_ROLE_DERIVED, &
            representation=SOURCE_REPRESENTATION_DEPENDENT, real_rate=.TRUE., &
            complex_rate=.FALSE., schedulable=.FALSE., &
            has_explicit=.FALSE., has_newton=.FALSE., &
            has_analytic=.FALSE.) &
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


  !> Look up a source by its canonical name. This is the first name-based
  !> interface; numeric catalog IDs remain temporarily for existing callers.
  !> An unknown name returns an empty descriptor (id=0).
  PURE FUNCTION source_term_descriptor_by_name(name) RESULT(term)
    CHARACTER(*), INTENT(IN) :: name
    TYPE(source_term_type) :: term
    INTEGER :: i

    term=source_term_type()
    DO i=1,SIZE(SOURCE_CATALOG)
      IF (TRIM(SOURCE_CATALOG(i)%name)==TRIM(name)) THEN
        term=SOURCE_CATALOG(i)
        RETURN
      END IF
    END DO
    ! Complete rheological laws have a single scheduling identity.
    DO i=1,SIZE(RHEOLOGY_LAWS)
      IF (TRIM(RHEOLOGY_LAWS(i)%name)/=TRIM(name)) CYCLE
      term=source_term_type()
      term%name=RHEOLOGY_LAWS(i)%name
      term%target=SOURCE_TARGET_VELOCITY
      term%model_id=RHEOLOGY_LAWS(i)%model_id
      term%schedulable=.TRUE.
      term%has_newton=.TRUE.
      term%has_analytic=.TRUE.
      RETURN
    END DO
  END FUNCTION source_term_descriptor_by_name

  !> Guard the name-only identity contract against duplicate catalog entries.
  PURE LOGICAL FUNCTION source_catalog_names_unique() RESULT(unique)
    INTEGER :: i,j
    unique=.FALSE.
    DO i=1,SIZE(SOURCE_CATALOG)
      IF (LEN_TRIM(SOURCE_CATALOG(i)%name)==0) RETURN
      DO j=i+1,SIZE(SOURCE_CATALOG)
        IF (TRIM(SOURCE_CATALOG(i)%name)==TRIM(SOURCE_CATALOG(j)%name)) RETURN
      END DO
    END DO
    unique=.TRUE.
  END FUNCTION source_catalog_names_unique


  !> Temporary name-based adapter for atomic legacy evaluation batches.
  FUNCTION source_bundle_names(bundle_id) RESULT(names)
    INTEGER, INTENT(IN) :: bundle_id
    CHARACTER(LEN=64), ALLOCATABLE :: names(:)
    INTEGER :: i,n
    ALLOCATE(names(SOURCE_TERM_COUNT))
    n=0
    IF (bundle_id>0) THEN
      DO i=1,SOURCE_TERM_COUNT
        IF (SOURCE_CATALOG(i)%bundle_id/=bundle_id) CYCLE
        n=n+1
        names(n)=SOURCE_CATALOG(i)%name
      END DO
    END IF
    names=names(:n)
  END FUNCTION source_bundle_names

  !> Resolve any legacy mathematical component to the identity of its
  !> complete physical source. Keep the compatibility metadata internal.
  PURE FUNCTION canonical_source_name(component_name) RESULT(name)
    CHARACTER(*), INTENT(IN) :: component_name
    CHARACTER(LEN=64) :: name
    TYPE(source_term_type) :: component
    component=source_term_descriptor_by_name(component_name)
    name=component%name
    IF (component%model_id/=0) name=rheology_source_name(component%model_id)
  END FUNCTION canonical_source_name

  !> Whether a complete rheological law has contributions of this internal
  !> assembly role. This is metadata, never a separately schedulable source.
  PURE LOGICAL FUNCTION rheology_has_role(model_id,role) RESULT(found)
    INTEGER, INTENT(IN) :: model_id,role
    INTEGER :: i
    found=.FALSE.
    IF (model_id==0) RETURN
    DO i=1,SOURCE_TERM_COUNT
      IF (SOURCE_CATALOG(i)%model_id/=model_id) CYCLE
      IF (SOURCE_CATALOG(i)%role==role) THEN
        found=.TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION rheology_has_role

  !> Canonical identity of an entire rheological law, not its formula addends.
  PURE FUNCTION rheology_source_name(model_id) RESULT(name)
    INTEGER, INTENT(IN) :: model_id
    CHARACTER(LEN=64) :: name
    INTEGER :: i
    name=''
    DO i=1,SIZE(RHEOLOGY_LAWS)
      IF (RHEOLOGY_LAWS(i)%model_id==model_id) THEN
        name=RHEOLOGY_LAWS(i)%name
        RETURN
      END IF
    END DO
  END FUNCTION rheology_source_name

  PURE LOGICAL FUNCTION canonical_rheology_catalog_unique() RESULT(ok)
    INTEGER :: i,j
    ok=.FALSE.
    DO i=1,SIZE(RHEOLOGY_LAWS)
      IF (LEN_TRIM(RHEOLOGY_LAWS(i)%name)==0) RETURN
      DO j=i+1,SIZE(RHEOLOGY_LAWS)
        IF (RHEOLOGY_LAWS(i)%name==RHEOLOGY_LAWS(j)%name) RETURN
        IF (RHEOLOGY_LAWS(i)%model_id==RHEOLOGY_LAWS(j)%model_id) RETURN
      END DO
    END DO
    ok=.TRUE.
  END FUNCTION canonical_rheology_catalog_unique


END MODULE source_catalog_2d
