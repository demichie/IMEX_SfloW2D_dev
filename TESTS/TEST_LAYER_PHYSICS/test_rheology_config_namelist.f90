PROGRAM test_rheology_config_namelist
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    IMPLICIT NONE

    TYPE(rheology_config_type) :: rheology, rheology_l1, rheology_l2
    NAMELIST /rheology_parameters/ rheology
    NAMELIST /rheology_l1_parameters/ rheology_l1
    NAMELIST /rheology_l2_parameters/ rheology_l2
    CHARACTER(LEN=180) :: legacy(4),lower(5),upper(5)
    INTEGER :: ios

    legacy=[CHARACTER(LEN=180) :: '&RHEOLOGY_PARAMETERS', &
        "RHEOLOGY%RHEOLOGY_MODEL='voellmy',", &
        'RHEOLOGY%MU=0.3, RHEOLOGY%XI=300.0,', '/']
    READ(legacy,NML=rheology_parameters,IOSTAT=ios)
    CALL check(ios==0,'single-layer derived-type namelist')
    rheology_l1=rheology
    CALL check(rheology_l1%rheology_model=='voellmy','copy legacy into L1')
    CALL check(close_real(rheology_l1%mu,0.3_wp) .AND. &
        close_real(rheology_l1%xi,300.0_wp),'L1 coefficients')

    lower=[CHARACTER(LEN=180) :: '&RHEOLOGY_L1_PARAMETERS', &
        "RHEOLOGY_L1%RHEOLOGY_MODEL='plastic',", &
        'RHEOLOGY_L1%TAU=120.0,', 'RHEOLOGY_L1%RHEOLOGY_GROUP=1,', '/']
    upper=[CHARACTER(LEN=180) :: '&RHEOLOGY_L2_PARAMETERS', &
        "RHEOLOGY_L2%RHEOLOGY_MODEL='froude',", &
        'RHEOLOGY_L2%MU_0=0.1, RHEOLOGY_L2%MU_INF=0.5,', &
        "RHEOLOGY_L2%FR_0=1.0, RHEOLOGY_L2%RHEOLOGY_METHOD='ANALYTIC',", '/']
    READ(lower,NML=rheology_l1_parameters,IOSTAT=ios)
    CALL check(ios==0,'layer 1 namelist')
    READ(upper,NML=rheology_l2_parameters,IOSTAT=ios)
    CALL check(ios==0,'layer 2 namelist')
    CALL check(close_real(rheology_l1%tau,120.0_wp),'L1 plastic coefficient')
    CALL check(close_real(rheology_l2%Fr_0,1.0_wp),'L2 Froude coefficient')
    CALL check(close_real(rheology_l2%mu,-1.0_wp) .AND. &
        close_real(rheology_l1%mu,0.3_wp), &
        'layer configurations are independent')
    PRINT *, 'PASS: NAMELIST derived-type rheology configuration L1/L2'
CONTAINS
    PURE LOGICAL FUNCTION close_real(actual,expected)
        REAL(wp), INTENT(IN) :: actual,expected
        close_real=ABS(actual-expected)<=32.0_wp*EPSILON(1.0_wp)* &
            MAX(1.0_wp,ABS(actual),ABS(expected))
    END FUNCTION close_real

    SUBROUTINE check(ok,label)
        LOGICAL, INTENT(IN) :: ok
        CHARACTER(*), INTENT(IN) :: label
        IF (.NOT.ok) THEN
            PRINT *, 'FAIL: ',label
            ERROR STOP 1
        END IF
    END SUBROUTINE check
END PROGRAM test_rheology_config_namelist
