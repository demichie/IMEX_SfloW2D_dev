!> @brief Verify both namelist spellings update one canonical L1 configuration.
!> @details Exercises the actual namelist groups declared by inpout_2d.
PROGRAM test_rheology_namelist_alias

    USE parameters_2d, ONLY : wp
    USE inpout_2d, ONLY : rheology_parameters, rheology_l1_parameters, rheology
    USE rheology_parameters_2d, ONLY : rheology_config_type, rheology_l1, rheology_l2

    IMPLICIT NONE

    INTEGER :: input_unit, ios
    REAL(wp), PARAMETER :: tolerance=1.0E-12_wp

    rheology_l1=rheology_config_type()
    rheology_l2=rheology_config_type()

    IF (.NOT.ASSOCIATED(rheology,rheology_l1)) ERROR STOP 'Legacy name is not an alias of L1'

    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='readwrite')
    WRITE(input_unit,'(A)') '&RHEOLOGY_PARAMETERS'
    WRITE(input_unit,'(A)') " RHEOLOGY%RHEOLOGY_MODEL='quadratic',"
    WRITE(input_unit,'(A)') ' RHEOLOGY%FRICTION_FACTOR=0.02D0,'
    WRITE(input_unit,'(A)') " RHEOLOGY%RHEOLOGY_METHOD='NEWTON',"
    WRITE(input_unit,'(A)') ' RHEOLOGY%RHEOLOGY_GROUP=1,'
    WRITE(input_unit,'(A)') '/'
    REWIND(input_unit)
    READ(input_unit,NML=rheology_parameters,IOSTAT=ios)
    CLOSE(input_unit)

    IF (ios/=0) ERROR STOP 'Could not read RHEOLOGY_PARAMETERS'
    IF (TRIM(rheology_l1%rheology_model)/='quadratic') ERROR STOP 'Legacy model mismatch'
    IF (ABS(rheology_l1%friction_factor-0.02_wp)>tolerance) ERROR STOP 'Legacy coefficient mismatch'
    IF (TRIM(rheology_l1%rheology_method)/='NEWTON') ERROR STOP 'Legacy method mismatch'
    IF (rheology_l1%rheology_group/=1) ERROR STOP 'Legacy group mismatch'
    IF (TRIM(rheology_l2%rheology_model)/='none') ERROR STOP 'L2 changed during legacy read'

    rheology_l1=rheology_config_type()
    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='readwrite')
    WRITE(input_unit,'(A)') '&RHEOLOGY_L1_PARAMETERS'
    WRITE(input_unit,'(A)') " RHEOLOGY_L1%RHEOLOGY_MODEL='voellmy',"
    WRITE(input_unit,'(A)') ' RHEOLOGY_L1%MU=0.3D0,'
    WRITE(input_unit,'(A)') ' RHEOLOGY_L1%XI=250.0D0,'
    WRITE(input_unit,'(A)') '/'
    REWIND(input_unit)
    READ(input_unit,NML=rheology_l1_parameters,IOSTAT=ios)
    CLOSE(input_unit)

    IF (ios/=0) ERROR STOP 'Could not read RHEOLOGY_L1_PARAMETERS'
    IF (TRIM(rheology_l1%rheology_model)/='voellmy') ERROR STOP 'L1 model mismatch'
    IF (ABS(rheology_l1%mu-0.3_wp)>tolerance) ERROR STOP 'L1 mu mismatch'
    IF (ABS(rheology_l1%xi-250.0_wp)>tolerance) ERROR STOP 'L1 xi mismatch'
    IF (TRIM(rheology%rheology_model)/='voellmy') ERROR STOP 'Legacy alias lost association'
    IF (TRIM(rheology_l2%rheology_model)/='none') ERROR STOP 'L2 changed during L1 read'

    WRITE(*,*) 'PASS: single-layer rheology NAMELIST aliases share the canonical L1 object'

END PROGRAM test_rheology_namelist_alias
