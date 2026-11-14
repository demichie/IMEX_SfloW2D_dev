!********************************************************************************
!> \brief Translate configured source activation into validated D17 metadata.
!> \details Runs once from physical source settings. Owns no
!>          physical state and performs no rate evaluation or temporal update.
!>          Per-cell wetness, packing and window gates remain in existing kernels.
!********************************************************************************
MODULE source_configuration_2d
    USE source_catalog_2d
    USE source_registry_2d, ONLY : source_registry_type, source_group_type
    USE source_input_2d
    USE parameters_2d, ONLY : wp, n_eqns, rheology_flag,                                                               &
            curvature_term_flag, bottom_radial_source_flag,                                        &
            bottom_fissural_source_flag, radial_source_flag,                                                           &
            liquid_flag, gas_flag, lava_flag,                                                &
            idx_pore
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE thermal_parameters_2d, ONLY : thermal_config_type
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: initialize_configured_sources
CONTAINS
    !> Identify active physical sources exactly once for both input syntaxes.
    !> Thermal activation is determined by the physical cooling flag.
    SUBROUTINE collect_active_source_names(rheology,mass_config,include_thermal,include_pore,active_names,error)
        TYPE(rheology_config_type), INTENT(IN) :: rheology
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_config
        LOGICAL, INTENT(IN) :: include_thermal,include_pore
        CHARACTER(LEN=64), ALLOCATABLE, INTENT(OUT) :: active_names(:)
        CHARACTER(*), INTENT(OUT) :: error
        ALLOCATE(active_names(0))
        error=''
        ! Exactly one continuous carrier may be active in the source registry.
        IF (gas_flag .AND. liquid_flag) THEN
        error='GAS_FLAG and LIQUID_FLAG are mutually exclusive'; RETURN

        END IF

        IF (lava_flag .AND. (gas_flag .OR. liquid_flag)) THEN
        error='LAVA_FLAG requires GAS_FLAG=F and LIQUID_FLAG=F'; RETURN

        END IF

        IF (.NOT.(lava_flag .OR. gas_flag .OR. liquid_flag)) THEN
        error='At least one physical phase flag must be true'; RETURN

        END IF

        IF (rheology_flag) THEN
        IF (LEN_TRIM(canonical_rheology_source_name(rheology%rheology_model))>0)                                    &
                active_names=[CHARACTER(LEN=64) :: active_names,                                                    &
                    canonical_rheology_source_name(rheology%rheology_model)]
        IF (LEN_TRIM(canonical_rheology_source_name(rheology%rheology_model))==0 .AND.                              &
                TRIM(rheology%rheology_model)/='none') THEN
            error='D17: unsupported configured rheology'; RETURN

        END IF

        END IF

        IF (include_thermal) active_names=[CHARACTER(LEN=64) :: active_names,                                          &
            'thermal_cooling']
        IF (include_pore) active_names=[CHARACTER(LEN=64) :: active_names,'pore_pressure_dissipation']
        IF (curvature_term_flag) active_names=[CHARACTER(LEN=64) :: active_names,'geometry_curvature']
        IF (bottom_radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_circular_injection']
        IF (bottom_fissural_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_fissural_injection']
        IF (radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'radial_injection']
        ! Lateral boundary injection is a spatial boundary condition, not a Lie cell source.
        IF (mass_config%erosion%coefficient>0.0_wp) THEN
                active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_erosion']

        END IF

        IF (mass_config%lava_deposition%flag) THEN

        IF (.NOT.lava_flag) THEN
            error='Lava deposition requires LAVA_FLAG=T'; RETURN

        END IF

        active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_lava_deposition']

        END IF

        IF (mass_config%deposition%settling_flag) THEN
        active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_deposition']

        END IF

        IF (liquid_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_liquid_loss']
        IF (mass_config%entrainment%flag) &
            active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_entrainment']
        IF (gas_flag .AND. mass_config%vaporization%liquid_vaporization_flag) &
            active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_vaporization']
        error=''
    END SUBROUTINE collect_active_source_names

    !> Create Lie groups directly from physical METHOD and GROUP configuration.
    SUBROUTINE initialize_configured_sources(registry,error,rheology,thermal,pore_pressure,mass_config, &
        pore_enabled)
        TYPE(source_registry_type), INTENT(INOUT) :: registry
        CHARACTER(*), INTENT(OUT) :: error
        TYPE(rheology_config_type), INTENT(IN) :: rheology
        TYPE(thermal_config_type), INTENT(IN) :: thermal
        TYPE(pore_pressure_config_type), INTENT(IN) :: pore_pressure
        TYPE(mass_exchange_config_type), INTENT(IN) :: mass_config
        LOGICAL, INTENT(IN) :: pore_enabled
        TYPE(source_group_type), ALLOCATABLE :: groups(:)
        TYPE(source_group_type) :: g,tmp
        TYPE(source_term_type) :: term
        CHARACTER(LEN=64), ALLOCATABLE :: names(:)
        CHARACTER(LEN=16) :: requested_method
        INTEGER :: i,j,k,id,method,pore,n,affected_idx
        INTEGER, ALLOCATABLE :: affected(:)
        ! A failed reconfiguration must never leave a previous registry active.
        CALL registry%finalize()
        CALL collect_active_source_names(rheology,mass_config,thermal%thermal_cooling_flag, &
            pore_enabled,names,error)
        IF (LEN_TRIM(error)>0) RETURN
        ALLOCATE(groups(1))
        groups(1)%id=0
        groups(1)%method=SOURCE_METHOD_EXPLICIT
        groups(1)%source_names=[CHARACTER(LEN=64)::]
        groups(1)%equation_ids=[(n,n=1,n_eqns)]

        DO i=1,SIZE(names)
        term=source_term_descriptor_by_name(names(i))
        method=term%default_method
        id=0
        requested_method=''
        BLOCK
            CHARACTER(LEN=64) :: active_names(1)
            active_names(1)=names(i)

        SELECT CASE(TRIM(active_names(1)))
        CASE('rheology_voellmy','rheology_plastic','rheology_lava',                                                    &
                'rheology_lahar','rheology_empirical','rheology_quadratic',                                            &
                'rheology_froude','rheology_lucas','rheology_mui',                                                     &
                'rheology_mui_voellmy','rheology_newmodel')
            requested_method=rheology%rheology_method
            id=rheology%rheology_group
        CASE('thermal_cooling')
            requested_method=thermal%thermal_method
            id=thermal%thermal_group
        CASE('pore_pressure_dissipation')
            requested_method=pore_pressure%pore_pressure_method
            id=pore_pressure%pore_pressure_group
        CASE('bottom_circular_injection','bottom_fissural_injection',                                                  &
                'radial_injection','geometry_curvature')

            SELECT CASE(TRIM(active_names(1)))
            CASE('bottom_circular_injection')
                requested_method=bottom_circular_method
                id=bottom_circular_group
            CASE('bottom_fissural_injection')
                requested_method=bottom_fissural_method
                id=bottom_fissural_group
            CASE('radial_injection')
                requested_method=radial_injection_method
                id=radial_injection_group
            CASE('geometry_curvature')
                requested_method=geometry_curvature_method
                id=geometry_curvature_group

            END SELECT
        CASE('mass_exchange_lava_deposition','mass_exchange_erosion',                                                  &
                'mass_exchange_deposition','mass_exchange_carrier_liquid_loss',                                        &
                'mass_exchange_entrainment',                                                                           &
                'mass_exchange_vaporization')
            requested_method=mass_config%method
            id=mass_config%group
        CASE DEFAULT
        END SELECT

        END BLOCK

        SELECT CASE(TRIM(requested_method))
        CASE('')
            CONTINUE
        CASE('EXPLICIT')
            method=SOURCE_METHOD_EXPLICIT
        CASE('NEWTON')
            method=SOURCE_METHOD_NEWTON
        CASE('ANALYTIC')
            method=SOURCE_METHOD_ANALYTIC
        CASE DEFAULT
            error='D19: invalid METHOD'; RETURN

        END SELECT

        IF ((lava_flag .AND. .NOT.term%allows_lava) .OR.                                                               &
                (gas_flag .AND. .NOT.term%allows_gas) .OR.                                                             &
                (liquid_flag .AND. .NOT.term%allows_liquid)) THEN
            error='D19: source incompatible with phase configuration'; RETURN

        END IF

        IF (id<0 .OR. (id==0 .AND. method/=SOURCE_METHOD_EXPLICIT)) THEN
            error='D19: GROUP 0 requires EXPLICIT'; RETURN

        END IF

        SELECT CASE(method)
        CASE(SOURCE_METHOD_EXPLICIT)
            IF (.NOT.term%has_explicit) THEN
                error='D19: unsupported EXPLICIT method'; RETURN

            END IF
        CASE(SOURCE_METHOD_NEWTON)
            IF (.NOT.term%has_newton) THEN
                error='D19: unsupported NEWTON method'; RETURN

            END IF
        CASE(SOURCE_METHOD_ANALYTIC)
            IF (.NOT.term%has_analytic) THEN
                error='D19: unsupported ANALYTIC method'; RETURN

            END IF

        END SELECT

        k=0

        DO j=1,SIZE(groups)
            IF (groups(j)%id==id) k=j

        END DO

        IF (k==0) THEN
            g%id=id
            g%method=method
            g%source_names=[CHARACTER(LEN=64)::]
            g%equation_ids=[INTEGER::]
            IF (method==SOURCE_METHOD_EXPLICIT) g%equation_ids=[(n,n=1,n_eqns)]
            groups=[groups,g]
            k=SIZE(groups)

        END IF

        IF (groups(k)%method/=method) THEN
            error='D19: GROUP cannot mix METHOD'; RETURN

        END IF

        IF (method==SOURCE_METHOD_ANALYTIC .AND. SIZE(groups(k)%source_names)>0) THEN
            error='D19: ANALYTIC group must be exclusive'; RETURN

        END IF

        groups(k)%source_names=[CHARACTER(LEN=64)::groups(k)%source_names,names(i)]

        IF (method==SOURCE_METHOD_NEWTON) THEN
            affected=source_newton_affected_equations(names(i))

            IF (SIZE(affected)==0) THEN
                error='D19: Newton source has no structurally affected equations'; RETURN

            END IF

            DO j=1,SIZE(affected)
                affected_idx=affected(j)
                IF (.NOT.ANY(groups(k)%equation_ids==affected_idx))                                                    &
                    groups(k)%equation_ids=[groups(k)%equation_ids,affected_idx]

            END DO

        END IF

        END DO
        ! Vaporization produces steam from material deposited in the same Lie
        ! transaction. A separate (or absent) deposition source has no rate to use.
        IF (ANY(names=='mass_exchange_vaporization') .AND.                                                             &
            .NOT.ANY(names=='mass_exchange_deposition')) THEN
        error='D19: vaporization requires deposition in the same mass-exchange group'; RETURN

        END IF

        DO i=1,SIZE(groups)

        DO j=i+1,SIZE(groups)

            IF (groups(j)%id<groups(i)%id) THEN
                tmp=groups(j); groups(j)=groups(i); groups(i)=tmp

            END IF

        END DO

        END DO

        pore=0
        IF (pore_enabled) pore=idx_pore
        CALL registry%initialize(n_eqns,pore,groups,error)
    END SUBROUTINE initialize_configured_sources

END MODULE source_configuration_2d
