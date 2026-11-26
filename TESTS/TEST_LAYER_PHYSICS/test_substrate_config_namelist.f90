!> @brief Verify canonical substrate NAMELIST input and shared erodible fractions.
PROGRAM test_substrate_config_namelist
    USE parameters_2d, ONLY : wp
    USE substrate_parameters_2d, ONLY : substrate_config_type, substrate, erodible_fract, coeff_porosity
    IMPLICIT NONE
    NAMELIST /substrate_parameters/ substrate, erodible_fract
    INTEGER :: unit, ios

    substrate=substrate_config_type()
    ALLOCATE(erodible_fract(2))
    erodible_fract=-1.0_wp

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&SUBSTRATE_PARAMETERS'
    WRITE(unit,'(A)') ' SUBSTRATE%T_GROUND=317.0, SUBSTRATE%ERODIBLE_POROSITY=0.25,'
    WRITE(unit,'(A)') ' SUBSTRATE%TOPOGRAPHY_FILE="DEM/topography.asc",'
    WRITE(unit,'(A)') ' SUBSTRATE%ERODIBLE_FILE="DEM/erodible.asc",'
    WRITE(unit,'(A)') ' SUBSTRATE%TOPO_CHANGE_FLAG=T, SUBSTRATE%ERODIBLE_DEPOSIT_FLAG=F,'
    WRITE(unit,'(A)') ' ERODIBLE_FRACT=0.6,0.4'
    WRITE(unit,'(A)') '/'
    REWIND(unit)
    READ(unit,NML=substrate_parameters,IOSTAT=ios)
    IF (ios/=0) ERROR STOP 'Canonical substrate NAMELIST read failed'
    IF (substrate%T_ground/=317.0_wp) ERROR STOP 'Incorrect ground temperature'
    IF (substrate%erodible_porosity/=0.25_wp) ERROR STOP 'Incorrect substrate porosity'
    IF (.NOT.substrate%topo_change_flag) ERROR STOP 'Incorrect topography policy'
    IF (substrate%erodible_deposit_flag) ERROR STOP 'Incorrect erodible-deposit policy'
    IF (TRIM(substrate%topography_file)/='DEM/topography.asc') ERROR STOP 'Incorrect DEM path'
    IF (TRIM(substrate%erodible_file)/='DEM/erodible.asc') ERROR STOP 'Incorrect erodible path'
    IF (ANY(ABS(erodible_fract-[0.6_wp,0.4_wp])>EPSILON(1.0_wp))) &
        ERROR STOP 'Incorrect erodible solid fractions'
    coeff_porosity=substrate%erodible_porosity/(1.0_wp-substrate%erodible_porosity)
    IF (coeff_porosity<=0.0_wp) ERROR STOP 'Derived porosity coefficient missing'
    CLOSE(unit)

    ! Retired input names must not be silently interpreted as aliases.
    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&SUBSTRATE_PARAMETERS T_ERODIBLE=400.0 /'
    REWIND(unit)
    READ(unit,NML=substrate_parameters,IOSTAT=ios)
    IF (ios==0) ERROR STOP 'Obsolete T_ERODIBLE input unexpectedly accepted'
    CLOSE(unit)

    PRINT *, 'PASS: canonical shared substrate and retired-name rejection'
END PROGRAM test_substrate_config_namelist
