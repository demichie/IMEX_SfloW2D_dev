!> @brief Test canonical L1/L2 mass-exchange NAMELIST acceptance and rejection.
PROGRAM test_mass_exchange_namelist_alias
    USE parameters_2d, ONLY : wp, n_layers
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type, mass_exchange_l1, mass_exchange_l2
    USE source_input_2d, ONLY : read_mass_exchange_parameters
    IMPLICIT NONE

    INTEGER :: unit
    CHARACTER(LEN=256) :: error

    n_layers=1
    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE%EROSION%COEFFICIENT=0.125,'
    WRITE(unit,'(A)') ' MASS_EXCHANGE%DEPOSITION%SETTLING_FLAG=T,'
    WRITE(unit,'(A)') ' MASS_EXCHANGE%METHOD="EXPLICIT", MASS_EXCHANGE%GROUP=4'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Short L1 mass-exchange NAMELIST rejected'
    IF (ABS(mass_exchange_l1%erosion%coefficient-0.125_wp)>EPSILON(1.0_wp)) &
        ERROR STOP 'Short NAMELIST did not bind canonical L1'
    IF (.NOT.mass_exchange_l1%deposition%settling_flag) ERROR STOP 'Nested flag lost'
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_L1_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L1%EROSION%COEFFICIENT=0.25,'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L1%GROUP=3'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Explicit L1 mass-exchange NAMELIST rejected'
    IF (ABS(mass_exchange_l1%erosion%coefficient-0.25_wp)>EPSILON(1.0_wp)) &
        ERROR STOP 'Explicit NAMELIST did not share canonical L1'
    IF (mass_exchange_l1%group/=3) ERROR STOP 'L1 group lost'
    IF (mass_exchange_l2%group/=4) ERROR STOP 'L2 configuration was modified'
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE%GROUP=4'
    WRITE(unit,'(A)') '/'
    WRITE(unit,'(A)') '&MASS_EXCHANGE_L1_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L1%GROUP=4'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (INDEX(error,'at most one')==0) ERROR STOP 'Duplicate L1 group accepted'
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_L2_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L2%GROUP=4'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (INDEX(error,'requires N_LAYERS=2')==0) ERROR STOP 'Premature L2 accepted'
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_PARAMETERS'
    WRITE(unit,'(A)') ' EROSION_COEFF=1.0'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (LEN_TRIM(error)==0) ERROR STOP 'Retired flat scalar accepted'
    CLOSE(unit)

    n_layers=2
    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE%GROUP=4'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (LEN_TRIM(error)==0) ERROR STOP 'Short group accepted for two layers'
    CLOSE(unit)

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&MASS_EXCHANGE_L1_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L1%EROSION%COEFFICIENT=0.75'
    WRITE(unit,'(A)') '/'
    WRITE(unit,'(A)') '&MASS_EXCHANGE_L2_PARAMETERS'
    WRITE(unit,'(A)') ' MASS_EXCHANGE_L2%EROSION%COEFFICIENT=0.5'
    WRITE(unit,'(A)') '/'
    CALL read_mass_exchange_parameters(unit,error)
    IF (LEN_TRIM(error)>0) ERROR STOP 'Two separate layer configurations rejected'
    IF (ABS(mass_exchange_l1%erosion%coefficient-0.75_wp)>EPSILON(1.0_wp)) &
        ERROR STOP 'Wrong L1 erosion coefficient'
    IF (ABS(mass_exchange_l2%erosion%coefficient-0.5_wp)>EPSILON(1.0_wp)) &
        ERROR STOP 'Wrong L2 erosion coefficient'
    CLOSE(unit)
    n_layers=1

    mass_exchange_l1=mass_exchange_config_type()
    mass_exchange_l2=mass_exchange_config_type()
    PRINT '(A)', 'PASS: mass-exchange single-layer NAMELIST aliases and rejection'
END PROGRAM test_mass_exchange_namelist_alias
