!> @brief Define the canonical identities and supported treatments of physical source terms.
!> @details Physical processes are selected by input configuration, but their admissible
!>          integration methods must be known before constructing any Lie source group.
!>          SOURCE_CATALOG provides that fixed lookup table: a canonical source name,
!>          explicit/Newton/analytic availability, the default method and permitted carrier phases.
!>          Entries include rheologies, cooling, pore-pressure dissipation, injections,
!>          geometric terms, environmental mass exchange and pair-wide interlayer processes.
!>
!>          During simulation initialization, source_configuration_2d resolves active physical
!>          laws against these descriptors and source_registry_2d validates their ordering,
!>          method compatibility and equation footprints. The resulting metadata is consumed
!>          by lie_source_executor_2d when advancing each cell. Numerical fluxes, actual
!>          source rates, Newton iterations and physical parameter values are NOT owned here:
!>          they belong to the spatial, physical-law and source-execution modules.
!>
!>          The catalogue IDs are internal dispatch and bitmask positions, not user-selected
!>          Lie GROUP numbers. The Newton equation sets describe structural coupling, not
!>          whether an instantaneous source rate happens to be nonzero.
MODULE source_catalog_2d
    USE parameters_2d, ONLY : idx_uEqn,idx_vEqn,idx_engyEqn,idx_poreEqn
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: source_term_type, SOURCE_CATALOG, SOURCE_TERM_COUNT
    PUBLIC :: source_term_descriptor_by_name, source_catalog_names_unique
    PUBLIC :: source_newton_affected_equations
    PUBLIC :: canonical_rheology_source_name, rheology_name_has_role
    PUBLIC :: rheology_source_id, source_term_id_by_name, source_id_mask
    ! Catalogue indices are stable internal dispatch IDs, not user-facing group IDs.
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_THERMAL_COOLING=11
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_BOTTOM_CIRCULAR_INJECTION=13
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_BOTTOM_FISSURAL_INJECTION=14
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_RADIAL_INJECTION=15
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_GEOMETRY_CURVATURE=16
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_LAVA_DEPOSITION=17
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_EROSION=18
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_DEPOSITION=19
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_CARRIER_LIQUID_LOSS=20
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_ENTRAINMENT=21
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_EXCHANGE_VAPORIZATION=22
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_INTERLAYER_MASS_EXCHANGE=23
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ID_INTERLAYER_DRAG=24
    INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_EXPLICIT=1, SOURCE_METHOD_NEWTON=2, SOURCE_METHOD_ANALYTIC=3
    INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_NEWTON_BLOCK=1, SOURCE_SOLVE_ANALYTIC_VELOCITY=2
    INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_TEMPERATURE=3, SOURCE_SOLVE_ANALYTIC_PORE_PRESSURE=4
    INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_IMPLICIT=1, SOURCE_ROLE_SEMI_IMPLICIT=2
    ! Equation-group membership is configured independently of source identity.
    !> @brief One entry in the fixed physical-source catalogue.
    !> @details The descriptor records which integration methods can be requested and which
    !>          carrier configurations may use a source. It contains no closure coefficients,
    !>          evolving state or callable source evaluator. Registry validation rejects a
    !>          configured method not enabled by these flags before the time integrator runs.
    TYPE :: source_term_type
        !> Unique catalog source name used to resolve configured Lie groups.
        CHARACTER(LEN=64) :: name = ''
        !> Whether an explicit rate evaluator is supported for this source.
        LOGICAL :: has_explicit = .FALSE.
        !> Whether this source supports the implicit Newton solver.
        LOGICAL :: has_newton = .FALSE.
        !> Whether this source supports an analytic time-integrated update.
        LOGICAL :: has_analytic = .FALSE.
        !> Default source integration method from the catalog method IDs.
        INTEGER :: default_method = SOURCE_METHOD_NEWTON
        ! Explicit per-entry permissions; absent permissions never enable a new source.
        !> Whether this source can operate with a gaseous carrier.
        LOGICAL :: allows_gas = .FALSE.
        !> Whether this source can operate with a liquid carrier.
        LOGICAL :: allows_liquid = .FALSE.
        !> Whether this source can operate in a lava-carrier configuration.
        LOGICAL :: allows_lava = .FALSE.
    END TYPE source_term_type
    ! Default method is explicit in each Lie source catalog entry; constants
    ! SOURCE_METHOD_* remain generic method IDs also used by the solvers.
    ! Rheology IDs correspond to the first entries of the canonical source catalog.
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_NONE=0
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_VOELLMY=1
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_PLASTIC=2
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_LAVA=3
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_LAHAR=4
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_EMPIRICAL=5
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_QUADRATIC=6
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_FROUDE=7
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_LUCAS=8
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_MUI=9
    INTEGER, PARAMETER, PUBLIC :: RHEOLOGY_ID_MUI_VOELLMY=10
    INTEGER, PARAMETER :: SOURCE_TERM_COUNT=24
    TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT)=[                                           &
            source_term_type('rheology_voellmy',.TRUE.,.TRUE.,.TRUE.,                                                  &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_plastic',.TRUE.,.TRUE.,.TRUE.,                                                  &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_lava',.TRUE.,.TRUE.,.TRUE.,                                                     &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.FALSE., allows_liquid=.TRUE., allows_lava=.TRUE.),                                    &
            source_term_type('rheology_lahar',.TRUE.,.TRUE.,.TRUE.,                                                    &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.FALSE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                    &
            source_term_type('rheology_empirical',.TRUE.,.TRUE.,.TRUE.,                                                &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_quadratic',.TRUE.,.TRUE.,.TRUE.,                                                &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_froude',.TRUE.,.FALSE.,.TRUE.,                                                  &
                    default_method=SOURCE_METHOD_EXPLICIT,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_lucas',.TRUE.,.FALSE.,.TRUE.,                                                   &
                    default_method=SOURCE_METHOD_EXPLICIT,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_mui',.TRUE.,.FALSE.,.TRUE.,                                                     &
                    default_method=SOURCE_METHOD_EXPLICIT,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('rheology_mui_voellmy',.TRUE.,.TRUE.,.TRUE.,                                              &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('thermal_cooling',.TRUE.,.TRUE.,.TRUE.,                                                   &
                    default_method=SOURCE_METHOD_NEWTON,                                                               &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.TRUE.),                                      &
            source_term_type('pore_pressure_dissipation',.FALSE.,.FALSE.,.TRUE.,                                       &
                    default_method=SOURCE_METHOD_ANALYTIC,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('bottom_circular_injection',.TRUE.,.FALSE.,.FALSE.,                                       &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('bottom_fissural_injection',.TRUE.,.FALSE.,.FALSE.,                                       &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.TRUE.),                                     &
            source_term_type('radial_injection',.TRUE.,.FALSE.,.FALSE.,                                                &
                    default_method=SOURCE_METHOD_EXPLICIT,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('geometry_curvature',.TRUE.,.FALSE.,.FALSE.,                                              &
                    default_method=SOURCE_METHOD_EXPLICIT,                                                             &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('mass_exchange_lava_deposition',.TRUE.,.FALSE.,.FALSE.,                                   &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.FALSE., allows_liquid=.FALSE., allows_lava=.TRUE.),                                    &
            source_term_type('mass_exchange_erosion',.TRUE.,.FALSE.,.FALSE.,                                           &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('mass_exchange_deposition',.TRUE.,.FALSE.,.FALSE.,                                        &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                     &
            source_term_type('mass_exchange_carrier_liquid_loss',.TRUE.,.FALSE.,.FALSE.,                               &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.FALSE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                    &
            source_term_type('mass_exchange_entrainment',.TRUE.,.FALSE.,.FALSE.,                                       &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.FALSE., allows_lava=.FALSE.),                                    &
            source_term_type('mass_exchange_vaporization',.TRUE.,.FALSE.,.FALSE.,                                      &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.FALSE., allows_lava=.FALSE.),                                    &
            source_term_type('interlayer_mass_exchange',.TRUE.,.FALSE.,.FALSE.,                                        &
                default_method=SOURCE_METHOD_EXPLICIT,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.),                                    &
            source_term_type('interlayer_drag',.FALSE.,.FALSE.,.TRUE., &
                default_method=SOURCE_METHOD_ANALYTIC,                                                                 &
                    allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.)                                     &
            ]
CONTAINS
    !> @brief Resolve a canonical source name to its catalog index during setup.
    !> @details Scans the fixed table in declaration order. The returned index is used as
    !>          a stable internal dispatch ID and as the bit position encoded by source_id_mask;
    !>          it must not be interpreted as the input GROUP or numerical method.
    !> @param[in] name Exact source identifier from the immutable registry.
    !> @return One-based catalog ID, or zero if no source matches.
    PURE INTEGER FUNCTION source_term_id_by_name(name) RESULT(id)
        CHARACTER(*), INTENT(IN) :: name
        INTEGER :: i
        id=0
        DO i=1,SOURCE_TERM_COUNT
            IF (TRIM(name)==TRIM(SOURCE_CATALOG(i)%name)) THEN
                id=i
                RETURN
            END IF
        END DO
    END FUNCTION source_term_id_by_name

    !> @brief Encode a source catalog ID as one bit of the active-source mask.
    !> @details The registry combines these bits to describe the configured members of a
    !>          group. The executor tests them to dispatch physical laws without repeatedly
    !>          comparing character names inside the OpenMP cell loop.
    !> @param[in] id One-based source ID to encode.
    !> @return The corresponding active bit, or zero when out of range.
    PURE INTEGER FUNCTION source_id_mask(id) RESULT(mask)
        INTEGER, INTENT(IN) :: id
        mask=0
        IF (id>=1 .AND. id<=SOURCE_TERM_COUNT) mask=IBSET(mask,id-1)
    END FUNCTION source_id_mask
    !> @brief Resolve the source ID of a named rheology at initialization.
    !> @param[in] name Canonical rheology name, or 'none' to disable rheology.
    !> @return Source ID; RHEOLOGY_ID_NONE for 'none', -1 when unknown.
    PURE INTEGER FUNCTION rheology_source_id(name) RESULT(id)
        CHARACTER(*), INTENT(IN) :: name
        INTEGER :: i
        id=-1
        IF (TRIM(name)=='none' .OR. LEN_TRIM(name)==0) THEN
            id=RHEOLOGY_ID_NONE
            RETURN
        END IF
        DO i=1,RHEOLOGY_ID_MUI_VOELLMY
            IF (TRIM(name)==TRIM(SOURCE_CATALOG(i)%name)) THEN
                id=i
                RETURN
            END IF
        END DO
    END FUNCTION rheology_source_id


    !> @brief Resolve a rheology model name to its canonical source identity.
    !> @details Converts a physical rheology model such as 'voellmy' to the corresponding
    !>          'rheology_voellmy' source identity. This separates the input model selection
    !>          from the names used by the source registry; it does not activate the law.
    !> @param[in] model_name Text identifier of the rheology constitutive law.
    !> @return Registered rheology source name, empty when not found or disabled.
    PURE FUNCTION canonical_rheology_source_name(model_name) RESULT(name)
        CHARACTER(*), INTENT(IN) :: model_name
        CHARACTER(LEN=64) :: name
        INTEGER :: i
        name=''
        IF (TRIM(model_name)=='none') RETURN

        DO i=1,SOURCE_TERM_COUNT
            IF (INDEX(TRIM(SOURCE_CATALOG(i)%name),'rheology_')/=1) CYCLE

            IF (TRIM(model_name)==TRIM(SOURCE_CATALOG(i)%name(10:))) THEN
                name=SOURCE_CATALOG(i)%name
                RETURN

            END IF

        END DO
    END FUNCTION canonical_rheology_source_name
    !> @brief Return conservative equation indices structurally affected by Newton.
    !> @details Newton groups solve only the conserved quantities coupled by their sources.
    !>          Supported rheology terms affect the two momentum equations; thermal cooling
    !>          affects thermal content. The returned footprint is independent of current cell
    !>          values, so a locally vanishing rate does not change the nonlinear system size.
    !>          source_configuration_2d unions these footprints for each proposed group, and
    !>          source_registry_2d verifies the exact set before nonlinear_solver_2d uses it.
    !>          Unknown or non-Newton sources return an empty set; an absent optional
    !>          equation is never inserted.
    !> @param[in] name Canonical source name.
    !> @return Allocatable list of conservative equation indices.
    FUNCTION source_newton_affected_equations(name) RESULT(ids)
        CHARACTER(*), INTENT(IN) :: name
        INTEGER, ALLOCATABLE :: ids(:)
        TYPE(source_term_type) :: descriptor
        descriptor=source_term_descriptor_by_name(name)
        ALLOCATE(ids(0))
        IF (.NOT.descriptor%has_newton) RETURN

        SELECT CASE(TRIM(name))
        CASE('rheology_voellmy','rheology_plastic','rheology_lava',                                                    &
            'rheology_lahar','rheology_empirical','rheology_quadratic',                                                &
            'rheology_mui_voellmy')
            ids=[idx_uEqn,idx_vEqn]
        CASE('thermal_cooling')
            ids=[idx_engyEqn]
        CASE('pore_pressure_dissipation')
            IF (idx_poreEqn>0) ids=[idx_poreEqn]

        END SELECT
    END FUNCTION source_newton_affected_equations


    !> @brief Look up one catalog descriptor by its canonical source name.
    !> @details Returns a value descriptor, not a pointer to mutable catalogue state.
    !>          Unknown names map to the default descriptor with all method and phase
    !>          capabilities disabled, allowing setup validation to reject the source.
    !> @param[in] name Exact registered source name.
    !> @return Descriptor, or default-disabled descriptor when unknown.
    PURE FUNCTION source_term_descriptor_by_name(name) RESULT(term)
        CHARACTER(*), INTENT(IN) :: name
        TYPE(source_term_type) :: term
        INTEGER :: i
        term=source_term_type()

        DO i=1,SOURCE_TERM_COUNT

            IF (TRIM(name)==TRIM(SOURCE_CATALOG(i)%name)) THEN
                term=SOURCE_CATALOG(i)
                RETURN

            END IF

        END DO
    END FUNCTION source_term_descriptor_by_name

    !> @brief Check uniqueness and nonempty names across source descriptors.
    !> @details Checks the fixed catalogue's identifying strings, independently of any
    !>          simulation input, to guard the name-to-ID resolution used by the registry.
    !> @return True if all registered source names are valid and unique.
    PURE LOGICAL FUNCTION source_catalog_names_unique() RESULT(ok)
        INTEGER :: i,j
        ok=.FALSE.

        DO i=1,SOURCE_TERM_COUNT
            IF (LEN_TRIM(SOURCE_CATALOG(i)%name)==0) RETURN

            DO j=i+1,SOURCE_TERM_COUNT
                IF (SOURCE_CATALOG(i)%name==SOURCE_CATALOG(j)%name) RETURN

            END DO

        END DO

        ok=.TRUE.
    END FUNCTION source_catalog_names_unique


    !> @brief Query whether a rheology implements a requested physical treatment.
    !> @details This query distinguishes the constitutive implicit and semi-implicit
    !>          contributions from the catalogue's group-level integration methods.
    !>          Callers use the result to select a supported rheological contribution;
    !>          the answer depends on the canonical source name rather than a model ID.
    !> @param[in] source_name Canonical rheology source name.
    !> @param[in] role Requested implicit or semi-implicit contribution identifier.
    !> @return True if that rheology supports the selected contribution.
    PURE LOGICAL FUNCTION rheology_name_has_role(source_name,role) RESULT(found)
        CHARACTER(*), INTENT(IN) :: source_name
        INTEGER, INTENT(IN) :: role
        found=.FALSE.

        SELECT CASE(role)
        CASE(SOURCE_ROLE_IMPLICIT)
            SELECT CASE(TRIM(source_name))
            CASE('rheology_voellmy','rheology_plastic','rheology_lava',                                                &
                'rheology_lahar','rheology_empirical','rheology_quadratic',                                            &
                'rheology_mui_voellmy')
                found=.TRUE.

            END SELECT
        CASE(SOURCE_ROLE_SEMI_IMPLICIT)
            SELECT CASE(TRIM(source_name))
            CASE('rheology_voellmy','rheology_lava','rheology_lahar',                                                  &
                'rheology_froude','rheology_lucas','rheology_mui',                                                     &
                'rheology_mui_voellmy')
                found=.TRUE.

            END SELECT

        END SELECT
    END FUNCTION rheology_name_has_role
END MODULE source_catalog_2d
