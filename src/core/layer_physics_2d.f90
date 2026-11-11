!> @brief Validated, per-layer phase and rheology selection.
!> @details Configuration is captured once before grid/worker allocation.
!>          A descriptor is read via queries and is never used to switch
!>          module-global phase/rheology parameters inside OpenMP loops.
!>          Physical coefficients and two-layer activation are later TL04 steps.
MODULE layer_physics_2d
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
        LOGICAL :: initialized=.FALSE.
    CONTAINS
        PROCEDURE :: initialize => initialize_layer_physics
        PROCEDURE :: finalize => finalize_layer_physics
        PROCEDURE :: is_initialized => layer_physics_initialized
        PROCEDURE :: phase => layer_physics_phase
        PROCEDURE :: selected_rheology => layer_physics_rheology
        PROCEDURE :: selected_rheology_id => layer_physics_rheology_id
    END TYPE layer_physics_type

CONTAINS

    SUBROUTINE initialize_layer_physics(this,phase,rheology_enabled,rheology_model,error)
        CLASS(layer_physics_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: phase
        LOGICAL, INTENT(IN) :: rheology_enabled
        CHARACTER(*), INTENT(IN) :: rheology_model
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
            IF (TRIM(rheology_model)/='none') THEN
                source=canonical_rheology_source_name(TRIM(rheology_model))
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
        this%carrier_phase=phase
        this%rheology_source=source
        this%rheology_id=RHEOLOGY_ID_NONE
        IF (LEN_TRIM(source)>0) this%rheology_id=rheology_source_id(source)
        this%initialized=.TRUE.
    END SUBROUTINE initialize_layer_physics

    SUBROUTINE finalize_layer_physics(this)
        CLASS(layer_physics_type), INTENT(INOUT) :: this
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

    PURE INTEGER FUNCTION layer_physics_rheology_id(this)
        CLASS(layer_physics_type), INTENT(IN) :: this
        layer_physics_rheology_id=this%rheology_id
    END FUNCTION layer_physics_rheology_id
END MODULE layer_physics_2d
