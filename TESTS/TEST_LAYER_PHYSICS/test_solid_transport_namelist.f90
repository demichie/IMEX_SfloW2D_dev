!> @brief Verify canonical shared solid species NAMELIST and derived density cache.
PROGRAM test_solid_transport_namelist
    USE parameters_2d, ONLY : wp
    USE solid_transport_parameters_2d, ONLY : rho_s, inv_rho_s, diam_s, &
        sphericity_s, sp_heat_s, maximum_solid_packing
    IMPLICIT NONE

    NAMELIST /solid_transport_parameters/ rho_s, diam_s, sphericity_s, &
        sp_heat_s, maximum_solid_packing

    INTEGER :: unit, ios
    ALLOCATE(rho_s(2),inv_rho_s(2),diam_s(2),sphericity_s(2),sp_heat_s(2))
    rho_s=-1.0_wp
    diam_s=-1.0_wp
    sp_heat_s=-1.0_wp
    sphericity_s=1.0_wp
    inv_rho_s=0.0_wp
    maximum_solid_packing=-1.0_wp

    OPEN(NEWUNIT=unit,STATUS='SCRATCH',ACTION='READWRITE')
    WRITE(unit,'(A)') '&SOLID_TRANSPORT_PARAMETERS'
    WRITE(unit,'(A)') ' RHO_S=2500.0,1800.0,'
    WRITE(unit,'(A)') ' DIAM_S=1.0E-4,2.0E-4,'
    WRITE(unit,'(A)') ' SPHERICITY_S=0.9,0.8,'
    WRITE(unit,'(A)') ' SP_HEAT_S=900.0,1100.0,'
    WRITE(unit,'(A)') ' MAXIMUM_SOLID_PACKING=0.61'
    WRITE(unit,'(A)') '/'
    REWIND(unit)
    READ(unit,NML=solid_transport_parameters,IOSTAT=ios)
    IF (ios/=0) ERROR STOP 'Canonical solid transport NAMELIST read failed'
    IF (ANY(ABS(rho_s-[2500.0_wp,1800.0_wp])>EPSILON(1.0_wp))) ERROR STOP 'Wrong solid density'
    IF (ANY(ABS(diam_s-[1.0E-4_wp,2.0E-4_wp])>EPSILON(1.0_wp))) ERROR STOP 'Wrong solid diameter'
    IF (ANY(ABS(sphericity_s-[0.9_wp,0.8_wp])>EPSILON(1.0_wp))) ERROR STOP 'Wrong solid sphericity'
    IF (ANY(ABS(sp_heat_s-[900.0_wp,1100.0_wp])>EPSILON(1.0_wp))) ERROR STOP 'Wrong solid heat capacity'
    IF (ABS(maximum_solid_packing-0.61_wp)>EPSILON(1.0_wp)) ERROR STOP 'Wrong common packing limit'
    IF (ANY(ABS(inv_rho_s)>TINY(1.0_wp))) ERROR STOP 'Derived cache was parsed from input'
    inv_rho_s=1.0_wp/rho_s
    IF (ANY(ABS(inv_rho_s*rho_s-1.0_wp)>16.0_wp*EPSILON(1.0_wp))) &
        ERROR STOP 'Wrong derived reciprocal solid density'
    CLOSE(unit)

    PRINT *, 'PASS: shared solid transport NAMELIST and derived density cache'
END PROGRAM test_solid_transport_namelist
