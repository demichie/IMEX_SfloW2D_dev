!********************************************************************************
!> \brief Native source-group namelist storage shared by input and translation.
!> \details Fixed capacities are bounded by the source catalog: a schedulable
!>          physical term can occur in only one group, and each implicit solve
!>          can target only the three supported equation blocks.
!********************************************************************************
MODULE source_input_2d
  USE source_catalog_2d, ONLY : SOURCE_TERM_COUNT
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER, PUBLIC :: SOURCE_INPUT_MAX_GROUPS=SOURCE_TERM_COUNT
  INTEGER, PARAMETER, PUBLIC :: SOURCE_INPUT_MAX_TERMS=SOURCE_TERM_COUNT
  INTEGER, PARAMETER, PUBLIC :: SOURCE_INPUT_MAX_TARGETS=3

  INTEGER, PUBLIC :: source_group_count=0
  INTEGER, PUBLIC :: source_group_id(SOURCE_INPUT_MAX_GROUPS)=0
  CHARACTER(LEN=16), PUBLIC :: source_group_cadence(SOURCE_INPUT_MAX_GROUPS)=''
  CHARACTER(LEN=16), PUBLIC :: source_group_method(SOURCE_INPUT_MAX_GROUPS)=''
  INTEGER, PUBLIC :: source_group_term_count(SOURCE_INPUT_MAX_GROUPS)=0
  CHARACTER(LEN=64), PUBLIC :: source_group_terms( &
       SOURCE_INPUT_MAX_TERMS,SOURCE_INPUT_MAX_GROUPS)=''
  ! Final group-level target input; legacy per-solver fields remain for migration.
  INTEGER, PUBLIC :: source_group_target_count(SOURCE_INPUT_MAX_GROUPS)=0
  CHARACTER(LEN=16), PUBLIC :: source_group_targets( &
       SOURCE_INPUT_MAX_TARGETS,SOURCE_INPUT_MAX_GROUPS)=''
  ! D18: physical-source scheduling preferences, absent unless specified.
  CHARACTER(LEN=16), PUBLIC :: rheology_method='', rheology_cadence=''
  CHARACTER(LEN=16), PUBLIC :: thermal_method='', thermal_cadence=''
  CHARACTER(LEN=16), PUBLIC :: pore_pressure_method='', pore_pressure_cadence=''
  INTEGER, PUBLIC :: rheology_group=0, thermal_group=0, pore_pressure_group=0
  CHARACTER(LEN=16), PUBLIC :: mass_exchange_method='EXPLICIT'
  CHARACTER(LEN=16), PUBLIC :: mass_exchange_cadence='STEP'
  INTEGER, PUBLIC :: mass_exchange_group=0
  LOGICAL, PUBLIC :: mass_exchange_configuration_present=.FALSE.
  NAMELIST /mass_exchange_parameters/ mass_exchange_method, &
       mass_exchange_cadence, mass_exchange_group
  LOGICAL, PUBLIC :: native_source_configuration_present=.FALSE.

  NAMELIST /source_groups/ source_group_count, source_group_id, &
       source_group_cadence, source_group_method, source_group_term_count, source_group_terms, &
       source_group_target_count, source_group_targets

  PUBLIC :: clear_source_input
  PUBLIC :: read_native_source_groups
  PUBLIC :: read_mass_exchange_parameters
  PUBLIC :: uppercase_source_text
CONTAINS
  !> \brief Reset all parsed native source-group input to the absent state.
  SUBROUTINE clear_source_input
    source_group_count=0
    source_group_id=0
    source_group_cadence=''
    source_group_method=''
    source_group_term_count=0
    source_group_terms=''
    source_group_target_count=0
    source_group_targets=''
    native_source_configuration_present=.FALSE.
  END SUBROUTINE clear_source_input

  !> \brief Find and parse the optional native source-group namelist.
  !> \param[in] input_unit Rewoundable unit containing the complete solver input.
  !> \param[in] legacy_imex_explicit True when a legacy IMEX field was set explicitly.
  !> \param[out] error Empty on success; otherwise describes malformed or conflicting input.
  SUBROUTINE read_native_source_groups(input_unit,legacy_imex_explicit,error)
    INTEGER, INTENT(IN) :: input_unit
    LOGICAL, INTENT(IN) :: legacy_imex_explicit
    CHARACTER(*), INTENT(OUT) :: error
    INTEGER :: io_status,namelist_count
    CHARACTER(LEN=2048) :: input_line

    CALL clear_source_input()
    error=''
    namelist_count=0
    REWIND(input_unit,IOSTAT=io_status)
    IF (io_status/=0) THEN
      error='Unable to rewind input while scanning SOURCE_GROUPS'
      RETURN
    END IF
    DO
      READ(input_unit,'(A)',IOSTAT=io_status) input_line
      IF (io_status<0) EXIT
      IF (io_status>0) THEN
        error='Unable to scan input for SOURCE_GROUPS'
        RETURN
      END IF
      input_line=ADJUSTL(input_line)
      IF (LEN_TRIM(input_line)==0) CYCLE
      IF (input_line(1:1)=='!') CYCLE
      input_line=uppercase_source_text(input_line)
      IF (INDEX(input_line,'&SOURCE_GROUPS')==1 .OR. &
           INDEX(input_line,'$SOURCE_GROUPS')==1) namelist_count=namelist_count+1
    END DO
    IF (namelist_count>1) THEN
      error='SOURCE_GROUPS namelist may appear only once'
      RETURN
    END IF
    IF (namelist_count==1) THEN
      REWIND(input_unit,IOSTAT=io_status)
      IF (io_status/=0) THEN
        error='Unable to rewind input before reading SOURCE_GROUPS'
        RETURN
      END IF
      READ(input_unit,NML=source_groups,IOSTAT=io_status)
      IF (io_status/=0) THEN
        ! A failed namelist read may have assigned only a prefix of fields.
        ! Do not expose an invalid partial configuration to subsequent reads.
        CALL clear_source_input()
        error='Invalid SOURCE_GROUPS namelist'
        RETURN
      END IF
      native_source_configuration_present=.TRUE.
    END IF
    IF (native_source_configuration_present .AND. legacy_imex_explicit) THEN
      error='SOURCE_GROUPS cannot be combined with explicit legacy IMEX_* settings'
      RETURN
    END IF
    REWIND(input_unit,IOSTAT=io_status)
    IF (io_status/=0) error='Unable to rewind input after reading SOURCE_GROUPS'
  END SUBROUTINE read_native_source_groups

  ! Read the optional indivisible mass-exchange schedule.
  SUBROUTINE read_mass_exchange_parameters(input_unit,error)
    INTEGER, INTENT(IN) :: input_unit
    CHARACTER(*), INTENT(OUT) :: error
    CHARACTER(LEN=2048) :: line
    INTEGER :: ios,count
    error=''
    count=0
    mass_exchange_configuration_present=.FALSE.
    REWIND(input_unit,IOSTAT=ios)
    IF (ios/=0) THEN
       error='Unable to scan MASS_EXCHANGE_PARAMETERS'; RETURN
    END IF
    DO
       READ(input_unit,'(A)',IOSTAT=ios) line
       IF (ios<0) EXIT
       IF (ios>0) THEN
          error='Unable to scan MASS_EXCHANGE_PARAMETERS'; RETURN
       END IF
       line=uppercase_source_text(ADJUSTL(line))
       IF (INDEX(line,'&MASS_EXCHANGE_PARAMETERS')==1 .OR. &
           INDEX(line,'$MASS_EXCHANGE_PARAMETERS')==1) count=count+1
    END DO
    IF (count>1) THEN
       error='MASS_EXCHANGE_PARAMETERS may appear only once'; RETURN
    END IF
    IF (count==1) THEN
       REWIND(input_unit)
       READ(input_unit,NML=mass_exchange_parameters,IOSTAT=ios)
       IF (ios/=0) THEN
          error='Invalid MASS_EXCHANGE_PARAMETERS namelist'; RETURN
       END IF
       mass_exchange_configuration_present=.TRUE.
       mass_exchange_method=uppercase_source_text(ADJUSTL(mass_exchange_method))
       mass_exchange_cadence=uppercase_source_text(ADJUSTL(mass_exchange_cadence))
       IF (mass_exchange_method/='EXPLICIT' .AND. &
            mass_exchange_method/='NEWTON' .AND. &
            mass_exchange_method/='ANALYTIC') THEN
          error='Unsupported mass-exchange METHOD'; RETURN
       END IF
       IF (mass_exchange_cadence/='STEP' .AND. mass_exchange_cadence/='STAGE') THEN
          error='Unsupported mass-exchange CADENCE'; RETURN
       END IF
       IF (mass_exchange_group<0) THEN
          error='Mass-exchange GROUP must not be negative'; RETURN
       END IF
    END IF
    REWIND(input_unit,IOSTAT=ios)
    IF (ios/=0) error='Unable to rewind after MASS_EXCHANGE_PARAMETERS'
  END SUBROUTINE read_mass_exchange_parameters

  !> \brief Convert ASCII lowercase letters in source input to uppercase.
  !> \param[in] text Input text.
  !> \return Same-length text with lowercase ASCII letters converted to uppercase.
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
END MODULE source_input_2d
