!> @brief Read and validate source-related input settings and scheduling metadata.
!> @details This module owns the input parsing of environmental mass exchange and
!>          pair-wide interlayer settings. It binds the single-layer mass-exchange
!>          NAMELIST and the explicit L1 NAMELIST to the same canonical L1 object,
!>          while L2 uses its own configuration. Pair-wide interlayer coefficients
!>          belong to the interlayer modules, not to either individual layer.
!>
!>          It also provides the METHOD/GROUP settings for injection and curvature
!>          sources. After parsing and normalization, source_configuration_2d uses
!>          these settings to build the ordered source identities; source_registry_2d
!>          validates method compatibility, and lie_source_executor_2d applies the
!>          selected physical laws. Parsing does not evaluate any cell-local source
!>          rates or mutate the evolving conservative state.
MODULE source_input_2d
    USE parameters_2d, ONLY : wp, n_layers
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type, mass_exchange_l1, mass_exchange_l2
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type, interlayer_mass_exchange
    USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange_config_type, &
        interlayer_momentum_exchange
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE
    CHARACTER(LEN=16), PUBLIC :: bottom_circular_method=''
    INTEGER, PUBLIC :: bottom_circular_group=0
    CHARACTER(LEN=16), PUBLIC :: bottom_fissural_method=''
    INTEGER, PUBLIC :: bottom_fissural_group=0
    CHARACTER(LEN=16), PUBLIC :: radial_injection_method=''
    INTEGER, PUBLIC :: radial_injection_group=0
    CHARACTER(LEN=16), PUBLIC :: geometry_curvature_method=''
    INTEGER, PUBLIC :: geometry_curvature_group=0

    ! Two input spellings write the same object; no copied or synchronized coefficients.
    TYPE(mass_exchange_config_type), POINTER :: mass_exchange => mass_exchange_l1
    NAMELIST /mass_exchange_parameters/ mass_exchange
    NAMELIST /mass_exchange_l1_parameters/ mass_exchange_l1
    NAMELIST /mass_exchange_l2_parameters/ mass_exchange_l2
    ! The pair-wide NAMELIST references the authoritative interlayer objects.
    NAMELIST /interlayer_parameters/ interlayer_mass_exchange,interlayer_momentum_exchange
    PUBLIC :: read_mass_exchange_parameters, write_mass_exchange_parameters
    PUBLIC :: read_interlayer_parameters, write_interlayer_parameters
    PUBLIC :: uppercase_source_text

CONTAINS

    !> @brief Read layer-local mass-exchange options and validate the combined source.
    !> @details Scans the input for the accepted mass-exchange NAMELIST forms,
    !>          rejecting duplicates and incompatible layer configurations.
    !>          A single-layer run may use the short or explicit L1 name, both of
    !>          which configure the same L1 physical object. A two-layer run
    !>          requires separate L1/L2 groups. Each selected configuration is
    !>          normalized and checked before source group construction; the input
    !>          unit is rewound for subsequent parsing. This routine does not decide
    !>          which cell exchanges material or evaluate a mass-transfer rate.
    !> @param[in] input_unit Open formatted input unit.
    !> @param[out] error Empty on success; error message on invalid input.
    SUBROUTINE read_mass_exchange_parameters(input_unit,error)
        INTEGER, INTENT(IN) :: input_unit
        CHARACTER(*), INTENT(OUT) :: error
        CHARACTER(LEN=2048) :: line
        INTEGER :: ios,n_short,n_l1,n_l2

        error=''
        n_short=0
        n_l1=0
        n_l2=0
        REWIND(input_unit,IOSTAT=ios)
        IF (ios/=0) THEN
            error='Unable to scan mass-exchange NAMELISTs'
            RETURN
        END IF

        DO
            READ(input_unit,'(A)',IOSTAT=ios) line
            IF (ios<0) EXIT
            IF (ios>0) THEN
                error='Unable to scan mass-exchange NAMELISTs'
                RETURN
            END IF

            line=uppercase_source_text(ADJUSTL(line))
            IF (INDEX(line,'&SOURCE_GROUPS')==1 .OR. INDEX(line,'$SOURCE_GROUPS')==1) THEN
                error='SOURCE_GROUPS is removed; use physical METHOD and GROUP'
                RETURN
            END IF
            IF (INDEX(line,'&MASS_EXCHANGE_PARAMETERS')==1 .OR. &
                INDEX(line,'$MASS_EXCHANGE_PARAMETERS')==1) n_short=n_short+1
            IF (INDEX(line,'&MASS_EXCHANGE_L1_PARAMETERS')==1 .OR. &
                INDEX(line,'$MASS_EXCHANGE_L1_PARAMETERS')==1) n_l1=n_l1+1
            IF (INDEX(line,'&MASS_EXCHANGE_L2_PARAMETERS')==1 .OR. &
                INDEX(line,'$MASS_EXCHANGE_L2_PARAMETERS')==1) n_l2=n_l2+1
        END DO

        IF (n_layers==1) THEN
            IF (n_l2/=0) THEN
                error='MASS_EXCHANGE_L2_PARAMETERS requires N_LAYERS=2'
                RETURN
            END IF
            IF (n_short+n_l1>1) THEN
                error='Provide at most one MASS_EXCHANGE_PARAMETERS or MASS_EXCHANGE_L1_PARAMETERS'
                RETURN
            END IF
            IF (n_short==1) THEN
                REWIND(input_unit)
                READ(input_unit,NML=mass_exchange_parameters,IOSTAT=ios)
            ELSEIF (n_l1==1) THEN
                REWIND(input_unit)
                READ(input_unit,NML=mass_exchange_l1_parameters,IOSTAT=ios)
            ELSE
                ios=0
            END IF
            IF (ios/=0) THEN
                error='Invalid mass-exchange L1 NAMELIST'
                RETURN
            END IF
            CALL validate_mass_exchange_configuration(mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) RETURN
        ELSE
            IF (n_short/=0 .OR. n_l1/=1 .OR. n_l2/=1) THEN
                error='N_LAYERS=2 requires exactly one mass-exchange L1 and L2 NAMELIST'
                RETURN
            END IF
            REWIND(input_unit)
            READ(input_unit,NML=mass_exchange_l1_parameters,IOSTAT=ios)
            IF (ios/=0) THEN
                error='Invalid mass-exchange L1 NAMELIST'
                RETURN
            END IF
            REWIND(input_unit)
            READ(input_unit,NML=mass_exchange_l2_parameters,IOSTAT=ios)
            IF (ios/=0) THEN
                error='Invalid mass-exchange L2 NAMELIST'
                RETURN
            END IF
            CALL validate_mass_exchange_configuration(mass_exchange_l1,error)
            IF (LEN_TRIM(error)>0) RETURN
            CALL validate_mass_exchange_configuration(mass_exchange_l2,error)
            IF (LEN_TRIM(error)>0) RETURN
        END IF

        REWIND(input_unit,IOSTAT=ios)
        IF (ios/=0) error='Unable to rewind after mass-exchange NAMELISTs'
    END SUBROUTINE read_mass_exchange_parameters

    !> @brief Normalize and validate the configured combined mass-exchange source.
    !> @param[in,out] config Canonical per-layer configuration being validated.
    !> @param[out] error Empty on success; diagnostic otherwise.
    SUBROUTINE validate_mass_exchange_configuration(config,error)
        TYPE(mass_exchange_config_type), INTENT(INOUT) :: config
        CHARACTER(*), INTENT(OUT) :: error

        error=''
        config%method=uppercase_source_text(ADJUSTL(config%method))
        IF (config%method/='EXPLICIT' .AND. config%method/='NEWTON' .AND. &
            config%method/='ANALYTIC') THEN
            error='Unsupported mass-exchange METHOD'
            RETURN
        END IF
        IF (config%group<0) THEN
            error='Mass-exchange GROUP must not be negative'
            RETURN
        END IF

        config%lava_deposition%model=uppercase_source_text(ADJUSTL(config%lava_deposition%model))
        IF (config%lava_deposition%flag) THEN
            SELECT CASE(TRIM(config%lava_deposition%model))
            CASE('CONSTANT','TEMPERATURE','VISCOSITY')
            CASE DEFAULT
                error='Unsupported LAVA_DEPOSITION_MODEL'
                RETURN
            END SELECT
            IF (config%lava_deposition%rate<0.0_wp .OR. &
                config%lava_deposition%visc_ref<=0.0_wp .OR. &
                config%lava_deposition%exponent<=0.0_wp .OR. &
                config%lava_deposition%temp_width<=0.0_wp) THEN
                error='Invalid lava deposition parameters'
                RETURN
            END IF
        END IF
    END SUBROUTINE validate_mass_exchange_configuration


    !> @brief Read canonical pair-wide mass exchange settings from one optional NAMELIST.
    !> @details Single-layer input must not contain interlayer settings. A single
    !>          pair owns these settings; neither L1 nor L2 is a surrogate owner.
    !> @param[in] input_unit Open Fortran unit containing the pair-wide interlayer NAMELIST.
    !> @param[out] error Validation or error message.
    SUBROUTINE read_interlayer_parameters(input_unit,error)
        INTEGER, INTENT(IN) :: input_unit
        CHARACTER(*), INTENT(OUT) :: error
        CHARACTER(LEN=2048) :: line
        INTEGER :: ios,count

        error=''
        count=0
        REWIND(input_unit,IOSTAT=ios)
        IF (ios/=0) THEN
            error='Unable to scan INTERLAYER_PARAMETERS'
            RETURN
        END IF

        DO
            READ(input_unit,'(A)',IOSTAT=ios) line
            IF (ios<0) EXIT
            IF (ios>0) THEN
                error='Unable to scan INTERLAYER_PARAMETERS'
                RETURN
            END IF
            line=uppercase_source_text(ADJUSTL(line))
            IF (INDEX(line,'&INTERLAYER_PARAMETERS')==1 .OR. &
                INDEX(line,'$INTERLAYER_PARAMETERS')==1) count=count+1
        END DO

        IF (count>1) THEN
            error='Duplicate INTERLAYER_PARAMETERS NAMELIST'
            RETURN
        END IF
        IF (n_layers/=2 .AND. count/=0) THEN
            error='INTERLAYER_PARAMETERS requires N_LAYERS=2'
            RETURN
        END IF
        IF (count==1) THEN
            REWIND(input_unit)
            READ(input_unit,NML=interlayer_parameters,IOSTAT=ios)
            IF (ios/=0) THEN
                error='Invalid INTERLAYER_PARAMETERS NAMELIST'
                RETURN
            END IF
        END IF

        CALL validate_interlayer_exchange_configuration(interlayer_mass_exchange,error)
        IF (LEN_TRIM(error)>0) RETURN
        CALL validate_interlayer_momentum_exchange_configuration(interlayer_momentum_exchange,error)
        IF (LEN_TRIM(error)>0) RETURN
        REWIND(input_unit,IOSTAT=ios)
        IF (ios/=0) error='Unable to rewind after interlayer NAMELIST'
    END SUBROUTINE read_interlayer_parameters

    !> @brief Validate pair-wide co-lift parameters without enabling an unimplemented source.
    !> @param[in,out] config Physical configuration.
    !> @param[out] error Validation or error message.
    SUBROUTINE validate_interlayer_exchange_configuration(config,error)
        TYPE(interlayer_mass_exchange_config_type), INTENT(INOUT) :: config
        CHARACTER(*), INTENT(OUT) :: error
        error=''
        config%method=uppercase_source_text(ADJUSTL(config%method))
        config%elutriation_carrier_mode=uppercase_source_text(ADJUSTL(config%elutriation_carrier_mode))
        config%elutriation_law=uppercase_source_text(ADJUSTL(config%elutriation_law))
        IF (config%method/='EXPLICIT' .OR. config%group<0) THEN
            error='Interlayer mass exchange currently supports EXPLICIT groups only'
            RETURN
        END IF
        SELECT CASE(TRIM(config%elutriation_carrier_mode))
        CASE('LIMIT_SOLIDS','ALLOW_CONCENTRATED')
        CASE DEFAULT
            error='Invalid ELUTRIATION_CARRIER_MODE'
            RETURN
        END SELECT
        IF (.NOT.ieee_is_finite(config%settling_velocity_weight) .OR. &
            .NOT.ieee_is_finite(config%settling_plug_thickness) .OR. &
            config%settling_velocity_weight<0.0_wp .OR. config%settling_velocity_weight>1.0_wp .OR. &
            config%settling_plug_thickness<=0.0_wp) THEN
            error='Invalid interlayer settling transfer velocity parameters'
            RETURN
        END IF
        IF (config%elutriation_flag) THEN
            SELECT CASE(TRIM(config%elutriation_law))
            CASE('LINEAR_USTAR','SHEAR_POWER')
            CASE DEFAULT
                error='Unsupported ELUTRIATION_LAW'
                RETURN
            END SELECT
            IF (.NOT.ieee_is_finite(config%elutriation_efficiency) .OR. &
                .NOT.ieee_is_finite(config%shear_interlayer) .OR. &
                .NOT.ieee_is_finite(config%shear_surface) .OR. &
                .NOT.ieee_is_finite(config%upper_depth_gate) .OR. &
                .NOT.ieee_is_finite(config%rouse_critical) .OR. &
                .NOT.ieee_is_finite(config%rouse_exponent) .OR. &
                .NOT.ieee_is_finite(config%bed_availability_exponent) .OR. &
                .NOT.ieee_is_finite(config%sediment_diffusivity_ratio) .OR. &
                .NOT.ieee_is_finite(config%elutriation_donor_cap)) THEN
                error='Nonfinite elutriation coefficient'
                RETURN
            END IF
            IF (config%elutriation_efficiency<0.0_wp .OR. config%shear_interlayer<0.0_wp .OR. &
                config%shear_surface<0.0_wp .OR. config%upper_depth_gate<=0.0_wp .OR. &
                config%rouse_critical<=0.0_wp .OR. config%rouse_exponent<=0.0_wp .OR. &
                config%bed_availability_exponent<=0.0_wp .OR. &
                config%sediment_diffusivity_ratio<=0.0_wp .OR. &
                config%elutriation_donor_cap<=0.0_wp .OR. config%elutriation_donor_cap>1.0_wp) THEN
                error='Invalid elutriation physical coefficient'
                RETURN
            END IF
            IF (config%elutriation_law=='SHEAR_POWER') THEN
                IF (.NOT.ieee_is_finite(config%shields_critical) .OR. &
                    .NOT.ieee_is_finite(config%shields_exponent) .OR. &
                    .NOT.ieee_is_finite(config%lift_length_diameters)) THEN
                    error='Invalid SHEAR_POWER coefficient'
                    RETURN
                END IF
                IF (config%shields_critical<0.0_wp .OR. config%shields_exponent<=0.0_wp .OR. &
                    config%lift_length_diameters<=0.0_wp) THEN
                    error='SHEAR_POWER requires explicit Shields and lift-length coefficients'
                    RETURN
                END IF
            END IF
        END IF
        IF (config%elutriation_flag) THEN
            IF (.NOT.ieee_is_finite(config%elutriation_solid_volume_fraction)) THEN
                error='Invalid ELUTRIATION_SOLID_VOLUME_FRACTION'
                RETURN
            END IF
            IF (config%elutriation_solid_volume_fraction/=-1.0_wp .AND. &
                (config%elutriation_solid_volume_fraction<=0.0_wp .OR. &
                config%elutriation_solid_volume_fraction>1.0_wp)) THEN
                error='ELUTRIATION_SOLID_VOLUME_FRACTION must be -1 (local L1) or in (0,1]'
                RETURN
            END IF
        END IF
    END SUBROUTINE validate_interlayer_exchange_configuration

    !> @brief Validate the independently scheduled pair-wide Chezy momentum source.
    !> @param[in,out] config Canonical momentum-exchange settings; method and model normalized.
    !> @param[out] error Empty on valid configuration, otherwise a descriptive diagnostic.
    SUBROUTINE validate_interlayer_momentum_exchange_configuration(config,error)
        TYPE(interlayer_momentum_exchange_config_type), INTENT(INOUT) :: config
        CHARACTER(*), INTENT(OUT) :: error
        error=''
        config%model=uppercase_source_text(ADJUSTL(config%model))
        config%method=uppercase_source_text(ADJUSTL(config%method))
        IF (config%model/='CHEZY') THEN
            error='Unsupported interlayer momentum exchange model'
            RETURN
        END IF
        IF (config%method/='ANALYTIC' .OR. config%group<=0) THEN
            error='Interlayer Chezy drag requires ANALYTIC method and positive GROUP'
            RETURN
        END IF
        IF (.NOT.ieee_is_finite(config%chezy_drag_coeff) .OR. config%chezy_drag_coeff<0.0_wp) THEN
            error='Invalid CHEZY_DRAG_COEFF'
            RETURN
        END IF
        IF (config%chezy_drag_coeff>0.0_wp .AND. &
            config%group==interlayer_mass_exchange%group .AND. &
            (interlayer_mass_exchange%settling_flag .OR. interlayer_mass_exchange%elutriation_flag)) THEN
            error='Interlayer mass and momentum exchange require different Lie groups'
            RETURN
        END IF
    END SUBROUTINE validate_interlayer_momentum_exchange_configuration

    !> @brief Write the unique canonical interlayer configuration for two-layer restarts.
    !> @param[in] output_unit Open output unit receiving the canonical interlayer NAMELIST.
    SUBROUTINE write_interlayer_parameters(output_unit)
        INTEGER, INTENT(IN) :: output_unit
        IF (n_layers==2) WRITE(output_unit,NML=interlayer_parameters)
    END SUBROUTINE write_interlayer_parameters

    !> @brief Convert ASCII lowercase source-input text to uppercase.
    !> @param[in] text Input characters.
    !> @return ASCII-uppercase text.
    PURE FUNCTION uppercase_source_text(text) RESULT(upper)
        CHARACTER(*), INTENT(IN) :: text
        CHARACTER(LEN=LEN(text)) :: upper
        INTEGER :: i,code

        upper=text
        DO i=1,LEN(text)
            code=IACHAR(text(i:i))
            IF (code>=IACHAR('a') .AND. code<=IACHAR('z')) &
                upper(i:i)=ACHAR(code-IACHAR('a')+IACHAR('A'))
        END DO
    END FUNCTION uppercase_source_text

    !> @brief Emit the canonical L1 configuration using the short NAMELIST name.
    !> @param[in] output_unit Open formatted output unit.
    SUBROUTINE write_mass_exchange_parameters(output_unit)
        INTEGER, INTENT(IN) :: output_unit
        WRITE(output_unit,NML=mass_exchange_parameters)
    END SUBROUTINE write_mass_exchange_parameters
END MODULE source_input_2d
