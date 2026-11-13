!> @brief Validated, per-layer phase and rheology selection.
!> @details Phase and rheology selection are validated before grid/worker allocation.
!>          A descriptor is read via queries and is never used to switch
!>          module-global phase/rheology parameters inside OpenMP loops.
!>          Rheology and thermal configurations are referenced without persistent copies.
!>          Two-layer activation remains a later step.
MODULE layer_physics_2d
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE thermal_parameters_2d, ONLY : thermal_config_type
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_config_type
    USE source_catalog_2d, ONLY : source_term_type, source_term_descriptor_by_name, &
        canonical_rheology_source_name, rheology_source_id, RHEOLOGY_ID_NONE
    IMPLICIT NONE
    PRIVATE
    INTEGER, PARAMETER, PUBLIC :: LAYER_PHASE_GAS=1
    INTEGER, PARAMETER, PUBLIC :: LAYER_PHASE_LIQUID=2
    INTEGER, PARAMETER, PUBLIC :: LAYER_PHASE_LAVA=3
    PUBLIC :: layer_physics_type

    TYPE :: layer_physics_type
        PRIVATE
        INTEGER :: carrier_phase=0
        INTEGER :: rheology_id=RHEOLOGY_ID_NONE
        CHARACTER(LEN=64) :: rheology_source=''
        !> Non-owning canonical rheology configuration.
        TYPE(rheology_config_type), POINTER :: config => NULL()
        !> Non-owning reference to the canonical thermal configuration.
        TYPE(thermal_config_type), POINTER :: thermal => NULL()
        !> Non-owning reference to the canonical layer Darcy parameters.
        TYPE(pore_pressure_config_type), POINTER :: pore_pressure => NULL()
        LOGICAL :: initialized=.FALSE.
    CONTAINS
        PROCEDURE :: initialize => initialize_layer_physics
        PROCEDURE :: finalize => finalize_layer_physics
        PROCEDURE :: is_initialized => layer_physics_initialized
        PROCEDURE :: phase => layer_physics_phase
        PROCEDURE :: selected_rheology => layer_physics_rheology
        PROCEDURE :: selected_rheology_id => layer_physics_rheology_id
        PROCEDURE :: rheology_configuration => layer_physics_configuration
        PROCEDURE :: thermal_configuration => layer_physics_thermal_configuration
        PROCEDURE :: pore_pressure_configuration => layer_physics_pore_pressure_configuration
    END TYPE layer_physics_type

CONTAINS

    !> @brief Initialize the validated layer descriptor and bind canonical thermal settings.
    !> @param[in,out] this Layer descriptor to initialize.
    !> @param[in] phase Material carrier phase identifier.
    !> @param[in] rheology_enabled Rheology activation for this layer.
    !> @param[in] rheology_config Stable-lifetime canonical rheology settings.
    !> @param[in] thermal_config Stable-lifetime canonical thermal settings.
    !> @param[in] pore_config Stable-lifetime canonical Darcy settings.
    !> @param[out] error Empty on success; diagnostic on failure.
    SUBROUTINE initialize_layer_physics(this,phase,rheology_enabled,rheology_config,thermal_config,pore_config,error)
        CLASS(layer_physics_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: phase
        LOGICAL, INTENT(IN) :: rheology_enabled
        TYPE(rheology_config_type), TARGET, INTENT(IN) :: rheology_config
        TYPE(thermal_config_type), TARGET, INTENT(IN) :: thermal_config
        TYPE(pore_pressure_config_type), TARGET, INTENT(IN) :: pore_config
        CHARACTER(*), INTENT(OUT) :: error
        CHARACTER(LEN=64) :: source
        TYPE(source_term_type) :: term
        CALL this%finalize()
        error=''
        IF (phase<LAYER_PHASE_GAS .OR. phase>LAYER_PHASE_LAVA) THEN
            error='TL04: invalid layer carrier phase'
            RETURN
        END IF
        source=''
        IF (rheology_enabled) THEN
            IF (TRIM(rheology_config%rheology_model)/='none') THEN
                source=canonical_rheology_source_name(TRIM(rheology_config%rheology_model))
                IF (LEN_TRIM(source)==0) THEN
                    error='TL04: unknown layer rheology'
                    RETURN
                END IF
                term=source_term_descriptor_by_name(source)
                SELECT CASE(phase)
                CASE(LAYER_PHASE_GAS)
                    IF (.NOT.term%allows_gas) error='TL04: rheology incompatible with gas'
                CASE(LAYER_PHASE_LIQUID)
                    IF (.NOT.term%allows_liquid) error='TL04: rheology incompatible with liquid'
                CASE(LAYER_PHASE_LAVA)
                    IF (.NOT.term%allows_lava) error='TL04: rheology incompatible with lava'
                END SELECT
                IF (LEN_TRIM(error)>0) RETURN
            END IF
        END IF
        this%config=>rheology_config
        this%thermal=>thermal_config
        this%pore_pressure=>pore_config
        this%carrier_phase=phase
        this%rheology_source=source
        this%rheology_id=RHEOLOGY_ID_NONE
        IF (LEN_TRIM(source)>0) this%rheology_id=rheology_source_id(source)
        this%initialized=.TRUE.
    END SUBROUTINE initialize_layer_physics

    SUBROUTINE finalize_layer_physics(this)
        CLASS(layer_physics_type), INTENT(INOUT) :: this
        NULLIFY(this%config)
        NULLIFY(this%thermal)
        NULLIFY(this%pore_pressure)
        this%carrier_phase=0
        this%rheology_id=RHEOLOGY_ID_NONE
        this%rheology_source=''
        this%initialized=.FALSE.
    END SUBROUTINE finalize_layer_physics

    PURE LOGICAL FUNCTION layer_physics_initialized(this)
        CLASS(layer_physics_type), INTENT(IN) :: this
        layer_physics_initialized=this%initialized
    END FUNCTION layer_physics_initialized

    PURE INTEGER FUNCTION layer_physics_phase(this)
        CLASS(layer_physics_type), INTENT(IN) :: this
        layer_physics_phase=this%carrier_phase
    END FUNCTION layer_physics_phase

    PURE FUNCTION layer_physics_rheology(this) RESULT(name)
        CLASS(layer_physics_type), INTENT(IN) :: this
        CHARACTER(LEN=64) :: name
        name=this%rheology_source
    END FUNCTION layer_physics_rheology

    ! Return a value view; callers cannot modify the canonical configuration through this query.
    PURE FUNCTION layer_physics_configuration(this) RESULT(config)
        CLASS(layer_physics_type), INTENT(IN) :: this
        TYPE(rheology_config_type) :: config
        config=this%config
    END FUNCTION layer_physics_configuration

    !> @brief Read the referenced canonical thermal configuration without exposing its pointer.
    !> @param[in] this Initialized layer descriptor.
    !> @return Temporary value view of the canonical thermal configuration.
    PURE FUNCTION layer_physics_thermal_configuration(this) RESULT(config)
        CLASS(layer_physics_type), INTENT(IN) :: this
        TYPE(thermal_config_type) :: config
        config=this%thermal
    END FUNCTION layer_physics_thermal_configuration

    !> @brief Read the canonical per-layer Darcy configuration as a value view.
    !> @param[in] this Initialized layer descriptor.
    !> @return Temporary Darcy configuration, independent of the canonical owner.
    PURE FUNCTION layer_physics_pore_pressure_configuration(this) RESULT(config)
        CLASS(layer_physics_type), INTENT(IN) :: this
        TYPE(pore_pressure_config_type) :: config
        config=this%pore_pressure
    END FUNCTION layer_physics_pore_pressure_configuration

    PURE INTEGER FUNCTION layer_physics_rheology_id(this)
        CLASS(layer_physics_type), INTENT(IN) :: this
        layer_physics_rheology_id=this%rheology_id
    END FUNCTION layer_physics_rheology_id
END MODULE layer_physics_2d
