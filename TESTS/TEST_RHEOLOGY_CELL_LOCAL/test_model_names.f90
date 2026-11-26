!> @brief Verify canonical rheology identities, source catalog ordering and masks.
PROGRAM test_rheology_model_names
    USE source_catalog_2d
    IMPLICIT NONE

    CHARACTER(LEN=32), PARAMETER :: names(11)=[CHARACTER(LEN=32) ::                                        &
        'none','voellmy','plastic','lava','lahar','empirical',                                              &
        'quadratic','froude','lucas','mui','mui_voellmy']
    CHARACTER(LEN=64), PARAMETER :: source_names(13)=[CHARACTER(LEN=64) ::                               &
        'thermal_cooling','bottom_circular_injection','bottom_fissural_injection',                         &
        'radial_injection','geometry_curvature','mass_exchange_lava_deposition',                            &
        'mass_exchange_erosion','mass_exchange_deposition','mass_exchange_carrier_liquid_loss',             &
        'mass_exchange_entrainment','mass_exchange_vaporization',                                          &
        'interlayer_mass_exchange','interlayer_drag']
    INTEGER, PARAMETER :: source_ids(13)=[SOURCE_ID_THERMAL_COOLING,                                      &
        SOURCE_ID_BOTTOM_CIRCULAR_INJECTION,SOURCE_ID_BOTTOM_FISSURAL_INJECTION,                           &
        SOURCE_ID_RADIAL_INJECTION,SOURCE_ID_GEOMETRY_CURVATURE,SOURCE_ID_EXCHANGE_LAVA_DEPOSITION,        &
        SOURCE_ID_EXCHANGE_EROSION,SOURCE_ID_EXCHANGE_DEPOSITION,                                         &
        SOURCE_ID_EXCHANGE_CARRIER_LIQUID_LOSS,SOURCE_ID_EXCHANGE_ENTRAINMENT,                             &
        SOURCE_ID_EXCHANGE_VAPORIZATION,SOURCE_ID_INTERLAYER_MASS_EXCHANGE,SOURCE_ID_INTERLAYER_DRAG]
    TYPE(source_term_type) :: descriptor
    CHARACTER(LEN=64) :: source_name
    INTEGER :: model

    DO model=1,SIZE(names)
        source_name=canonical_rheology_source_name(names(model))

        IF (TRIM(names(model))=='none') THEN
            IF (LEN_TRIM(source_name)/=0) ERROR STOP 'none is not disabled'
        ELSE
            IF (TRIM(source_name)/='rheology_'//TRIM(names(model))) &
                ERROR STOP 'incorrect canonical rheology name'
            descriptor=source_term_descriptor_by_name(source_name)
            IF (descriptor%name/=source_name) ERROR STOP 'missing canonical source'
        END IF
    END DO

    ! Source ID values are internal, but each exported constant must index its own catalog entry.
    DO model=1,SIZE(source_ids)
        IF (TRIM(SOURCE_CATALOG(source_ids(model))%name)/=TRIM(source_names(model))) &
            ERROR STOP 'source ID mismatches catalog position'
    END DO

    DO model=1,SOURCE_TERM_COUNT
        IF (source_term_id_by_name(SOURCE_CATALOG(model)%name)/=model) &
            ERROR STOP 'source name to ID lookup mismatch'
        IF (source_id_mask(model)/=IBSET(0,model-1)) ERROR STOP 'source mask mismatch'
    END DO

    IF (source_id_mask(SOURCE_TERM_COUNT+1)/=0) ERROR STOP 'out-of-range source mask is nonzero'
    IF (LEN_TRIM(canonical_rheology_source_name('newmodel'))/=0) &
        ERROR STOP 'retired development fixture is still registered'
    IF (LEN_TRIM(canonical_rheology_source_name('unknown-law'))/=0) &
        ERROR STOP 'unknown rheology accepted'
END PROGRAM test_rheology_model_names
