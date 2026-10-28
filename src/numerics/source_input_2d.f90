!********************************************************************************
!> \brief Physical-source scheduling input.
!> \details Physical METHOD, CADENCE and GROUP are specified in the owning
!>          namelists; SOURCE_GROUPS is no longer accepted.
!********************************************************************************
MODULE source_input_2d
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       loss_rate, entrainment_flag
  USE parameters_2d, ONLY : wp, gas_loss_flag, liquid_vaporization_flag
  IMPLICIT NONE
  PRIVATE
  ! D18: physical-source scheduling preferences, absent unless specified.
  CHARACTER(LEN=16), PUBLIC :: rheology_method='', rheology_cadence=''
  CHARACTER(LEN=16), PUBLIC :: thermal_method='', thermal_cadence=''
  CHARACTER(LEN=16), PUBLIC :: pore_pressure_method='', pore_pressure_cadence=''
  INTEGER, PUBLIC :: rheology_group=0, thermal_group=0, pore_pressure_group=0
  CHARACTER(LEN=16), PUBLIC :: bottom_circular_method='', bottom_circular_cadence=''
  INTEGER, PUBLIC :: bottom_circular_group=0
  CHARACTER(LEN=16), PUBLIC :: bottom_fissural_method='', bottom_fissural_cadence=''
  INTEGER, PUBLIC :: bottom_fissural_group=0
  CHARACTER(LEN=16), PUBLIC :: radial_injection_method='', radial_injection_cadence=''
  INTEGER, PUBLIC :: radial_injection_group=0
  CHARACTER(LEN=16), PUBLIC :: geometry_curvature_method='', geometry_curvature_cadence=''
  INTEGER, PUBLIC :: geometry_curvature_group=0
  CHARACTER(LEN=16), PUBLIC :: boundary_injection_method='', boundary_injection_cadence=''
  INTEGER, PUBLIC :: boundary_injection_group=0
  CHARACTER(LEN=16), PUBLIC :: mass_exchange_method='EXPLICIT'
  CHARACTER(LEN=16), PUBLIC :: mass_exchange_cadence='STEP'
  INTEGER, PUBLIC :: mass_exchange_group=0
  LOGICAL, PUBLIC :: lava_deposition_flag=.FALSE.
  CHARACTER(LEN=16), PUBLIC :: lava_deposition_model='CONSTANT'
  REAL(wp), PUBLIC :: lava_deposition_rate=0.0_wp
  REAL(wp), PUBLIC :: lava_deposition_visc_ref=1.0E-2_wp
  REAL(wp), PUBLIC :: lava_deposition_exponent=1.0_wp
  REAL(wp), PUBLIC :: lava_deposition_temp_ref=900.0_wp
  REAL(wp), PUBLIC :: lava_deposition_temp_width=100.0_wp
  LOGICAL, PUBLIC :: mass_exchange_configuration_present=.FALSE.
  NAMELIST /mass_exchange_parameters/ mass_exchange_method, &
       mass_exchange_cadence, mass_exchange_group, erosion_coeff, &
       settling_flag, loss_rate, entrainment_flag, gas_loss_flag, &
       liquid_vaporization_flag, lava_deposition_flag, lava_deposition_model, &
       lava_deposition_rate, lava_deposition_visc_ref, lava_deposition_exponent, &
       lava_deposition_temp_ref, lava_deposition_temp_width
  PUBLIC :: read_mass_exchange_parameters, write_mass_exchange_parameters
  PUBLIC :: eval_lava_deposition_rate, uppercase_source_text
CONTAINS
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
       IF (INDEX(line,'&SOURCE_GROUPS')==1 .OR. &
           INDEX(line,'$SOURCE_GROUPS')==1) THEN
          error='SOURCE_GROUPS is removed; use physical METHOD, CADENCE, GROUP'
          RETURN
       END IF
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
       lava_deposition_model=uppercase_source_text(ADJUSTL(lava_deposition_model))
       IF (lava_deposition_flag) THEN
          SELECT CASE(TRIM(lava_deposition_model))
          CASE('CONSTANT','TEMPERATURE','VISCOSITY')
          CASE DEFAULT
             error='Unsupported LAVA_DEPOSITION_MODEL'; RETURN
          END SELECT
          IF (lava_deposition_rate<0.0_wp .OR. lava_deposition_visc_ref<=0.0_wp .OR. &
               lava_deposition_exponent<=0.0_wp .OR. lava_deposition_temp_width<=0.0_wp) THEN
             error='Invalid lava deposition parameters'; RETURN
          END IF
       END IF
       IF (mass_exchange_group<0) THEN
          error='Mass-exchange GROUP must not be negative'; RETURN
       END IF
    END IF
    REWIND(input_unit,IOSTAT=ios)
    IF (ios/=0) error='Unable to rewind after MASS_EXCHANGE_PARAMETERS'
  END SUBROUTINE read_mass_exchange_parameters

  SUBROUTINE write_mass_exchange_parameters(output_unit)
    INTEGER, INTENT(IN) :: output_unit
    WRITE(output_unit,NML=mass_exchange_parameters)
  END SUBROUTINE write_mass_exchange_parameters

  !> \brief Requested immobile lava deposition thickness rate.
  !> \details Caller enforces available thickness over dt.
  !> \param[in] temperature Current mobile temperature [K].
  !> \param[in] viscosity Lava kinematic viscosity [m2/s].
  !> \return Prescribed deposition rate [m/s].
  PURE REAL(wp) FUNCTION eval_lava_deposition_rate(temperature,viscosity) RESULT(rate)
    REAL(wp), INTENT(IN) :: temperature,viscosity
    REAL(wp) :: factor
    factor=1.0_wp
    SELECT CASE(TRIM(lava_deposition_model))
    CASE('TEMPERATURE')
       factor=MIN(1.0_wp,MAX(0.0_wp,(lava_deposition_temp_ref-temperature)/ &
            lava_deposition_temp_width))
    CASE('VISCOSITY')
       factor=MIN(1.0_wp,(MAX(0.0_wp,viscosity)/lava_deposition_visc_ref)** &
            lava_deposition_exponent)
    END SELECT
    rate=MAX(0.0_wp,lava_deposition_rate)*factor
  END FUNCTION eval_lava_deposition_rate

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
