PROGRAM test_pore_pressure_namelist_alias
    USE parameters_2d, ONLY : wp
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_config_type, pore_pressure_l1, pore_pressure_l2, &
        pore_pressure_cache_l1, precompute_pascal_coefficient
    USE inpout_2d, ONLY : pore_pressure, pore_pressure_parameters, pore_pressure_l1_parameters, &
        pore_pressure_l2_parameters
    IMPLICIT NONE
    INTEGER :: ios, input_unit
    REAL(wp), PARAMETER :: tolerance=1.0e-14_wp

    pore_pressure_l1=pore_pressure_config_type()
    pore_pressure_l2=pore_pressure_config_type()
    IF (.NOT.ASSOCIATED(pore_pressure,pore_pressure_l1)) ERROR STOP 'Legacy Darcy alias does not target L1'

    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='READWRITE')
    WRITE(input_unit,'(A)') '&PORE_PRESSURE_PARAMETERS'
    WRITE(input_unit,'(A)') ' PORE_PRESSURE%HYDRAULIC_PERMEABILITY=1.0D-12,'
    WRITE(input_unit,'(A)') " PORE_PRESSURE%F_INHIBIT_MODE='STATIC',"
    WRITE(input_unit,'(A)') ' PORE_PRESSURE%N_INH=2, PORE_PRESSURE%ALPHA_TRANS=0.2,'
    WRITE(input_unit,'(A)') " PORE_PRESSURE%PORE_PRESSURE_METHOD='ANALYTIC',"
    WRITE(input_unit,'(A)') ' PORE_PRESSURE%PORE_PRESSURE_GROUP=3,'
    WRITE(input_unit,'(A)') '/'
    REWIND(input_unit)
    READ(input_unit,NML=pore_pressure_parameters,IOSTAT=ios)
    CLOSE(input_unit)
    IF (ios/=0) ERROR STOP 'Cannot parse legacy PORE_PRESSURE_PARAMETERS'
    IF (ABS(pore_pressure_l1%hydraulic_permeability-1.0e-12_wp)>tolerance*1.0e-12_wp) &
        ERROR STOP 'Legacy Darcy permeability mismatch'
    IF (pore_pressure_l1%N_inh/=2) ERROR STOP 'Legacy Darcy N_inh mismatch'
    IF (pore_pressure_l1%f_inhibit_mode/='STATIC') ERROR STOP 'Legacy Darcy inhibition mismatch'
    IF (pore_pressure_l1%pore_pressure_method/='ANALYTIC') ERROR STOP 'Legacy Darcy method mismatch'
    IF (pore_pressure_l1%pore_pressure_group/=3) ERROR STOP 'Legacy Darcy group mismatch'

    CALL precompute_pascal_coefficient(pore_pressure_l1,pore_pressure_cache_l1)
    IF (.NOT.ALLOCATED(pore_pressure_cache_l1%pascal_coeff)) ERROR STOP 'Pascal cache not initialized'
    IF (LBOUND(pore_pressure_cache_l1%pascal_coeff,1)/=0 .OR. &
        UBOUND(pore_pressure_cache_l1%pascal_coeff,1)/=2) ERROR STOP 'Pascal cache bounds mismatch'

    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='READWRITE')
    WRITE(input_unit,'(A)') '&PORE_PRESSURE_L1_PARAMETERS'
    WRITE(input_unit,'(A)') ' PORE_PRESSURE_L1%HYDRAULIC_PERMEABILITY=3.0D-12,'
    WRITE(input_unit,'(A)') '/'
    REWIND(input_unit)
    READ(input_unit,NML=pore_pressure_l1_parameters,IOSTAT=ios)
    CLOSE(input_unit)
    IF (ios/=0) ERROR STOP 'Cannot parse PORE_PRESSURE_L1_PARAMETERS'
    IF (pore_pressure%hydraulic_permeability/=pore_pressure_l1%hydraulic_permeability) &
        ERROR STOP 'Darcy L1 broke legacy alias'
    OPEN(NEWUNIT=input_unit,STATUS='scratch',ACTION='READWRITE')
    WRITE(input_unit,'(A)') '&PORE_PRESSURE_L2_PARAMETERS'
    WRITE(input_unit,'(A)') ' PORE_PRESSURE_L2%HYDRAULIC_PERMEABILITY=5.0D-12,'
    WRITE(input_unit,'(A)') '/'
    REWIND(input_unit)
    READ(input_unit,NML=pore_pressure_l2_parameters,IOSTAT=ios)
    CLOSE(input_unit)
    IF (ios/=0) ERROR STOP 'Cannot parse PORE_PRESSURE_L2_PARAMETERS'
    IF (pore_pressure_l2%hydraulic_permeability==pore_pressure_l1%hydraulic_permeability) &
        ERROR STOP 'Darcy L1 and L2 are not independent'
    PRINT *, 'PASS: canonical Darcy NAMELIST L1/L2 and legacy alias'
END PROGRAM test_pore_pressure_namelist_alias
