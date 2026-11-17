!> @brief Cell-local paired transfer budgets and co-lift input regressions.
PROGRAM test_interlayer_exchange
    USE parameters_2d, ONLY : wp,n_layers
    USE interlayer_exchange_parameters_2d, ONLY : interlayer_exchange_config_type, interlayer_exchange
    USE interlayer_exchange_2d, ONLY : limit_interlayer_mass_transfer
    USE source_input_2d, ONLY : read_interlayer_exchange_parameters
    IMPLICIT NONE
    TYPE(interlayer_exchange_config_type) :: config
    REAL(wp) :: s1(2),s2(2),settle(2),elut(2),accepted_s(2),accepted_e(2)
    REAL(wp) :: down,up,dt,carrier1,carrier2,pack,strict_lift,concentrated_lift
    CHARACTER(LEN=256) :: error
    INTEGER :: unit

    dt=1.0_wp
    pack=0.6_wp
    config%settling_flag=.TRUE.
    config%elutriation_flag=.TRUE.
    config%elutriation_solid_volume_fraction=0.1_wp

    ! Settling creates L1 from a dry state: carrier comes from L2, not nowhere.
    s1=0.0_wp; carrier1=0.0_wp
    s2=[0.12_wp,0.08_wp]; carrier2=0.5_wp
    settle=[0.03_wp,0.02_wp]; elut=0.0_wp
    CALL execute(.TRUE.)
    CALL require(ALL(accepted_s>0.0_wp) .AND. ALL(accepted_e==0.0_wp),'L1 settling cold-start')
    CALL close(down,(1.0_wp-pack)/pack*SUM(accepted_s),'minimal downward carrier')
    CALL close(up,0.0_wp,'no settling co-lift')

    ! Elutriation creates L2 from a dry state at the parcel composition.
    s1=[0.12_wp,0.08_wp]; carrier1=1.0_wp
    s2=0.0_wp; carrier2=0.0_wp
    settle=0.0_wp; elut=[0.01_wp,0.02_wp]
    CALL execute(.TRUE.)
    CALL close(up,9.0_wp*SUM(accepted_e),'prescribed elutriated parcel composition')
    CALL close(down,0.0_wp,'no need for downward carrier on elutriation')

    ! At limited carrier inventory, LIMIT_SOLIDS preserves parcel composition
    ! and also preserves the donor packing constraint.
    s1=[0.2_wp,0.0_wp]; carrier1=0.8_wp
    s2=0.0_wp; carrier2=0.0_wp
    settle=0.0_wp; elut=[0.2_wp,0.0_wp]
    CALL execute(.TRUE.)
    strict_lift=SUM(accepted_e)
    CALL require(strict_lift>0.0_wp .AND. strict_lift<0.2_wp,'strict mode limits solids')
    CALL close(up,9.0_wp*strict_lift,'strict mode keeps solid fraction')

    ! Relaxing parcel composition can increase solids lifted, but never permits
    ! packing violation or donor over-withdrawal.
    config%elutriation_carrier_mode='ALLOW_CONCENTRATED'
    CALL execute(.TRUE.)
    concentrated_lift=SUM(accepted_e)
    CALL require(concentrated_lift>strict_lift,'concentrated mode permits more solid lift')
    CALL require(up<=carrier1/dt,'gross L1 carrier inventory respected')
    CALL require(up<9.0_wp*concentrated_lift,'concentrated mode relaxes parcel composition')

    ! The coordinator limits simultaneous opposed transfers on the same snapshot.
    config%elutriation_carrier_mode='LIMIT_SOLIDS'
    s1=[0.2_wp,0.05_wp]; carrier1=0.3_wp
    s2=[0.12_wp,0.05_wp]; carrier2=0.2_wp
    settle=[0.1_wp,0.02_wp]; elut=[0.08_wp,0.01_wp]
    CALL execute(.TRUE.)
    CALL require(SUM(accepted_s)>0.0_wp .AND. SUM(accepted_e)>0.0_wp, &
        'both physical transfers remain possible')
    CALL require(down>=0.0_wp .AND. up>=0.0_wp,'two-way carrier fluxes nonnegative')

    ! Zero fluxes and zero layer inventories are valid.
    s1=0.0_wp; s2=0.0_wp; carrier1=0.0_wp; carrier2=0.0_wp
    settle=0.0_wp; elut=0.0_wp
    CALL execute(.TRUE.)
    CALL require(ALL(accepted_s==0.0_wp) .AND. ALL(accepted_e==0.0_wp) &
        .AND. down==0.0_wp .AND. up==0.0_wp,'fully dry column is stable')

    ! Invalid parameters fail closed, and no previously accepted rate survives.
    config%elutriation_solid_volume_fraction=0.0_wp
    CALL execute(.FALSE.)
    CALL require(ALL(accepted_s==0.0_wp) .AND. ALL(accepted_e==0.0_wp), &
        'invalid configuration resets transfer rates')
    config%elutriation_solid_volume_fraction=0.1_wp
    config%elutriation_carrier_mode='WRONG'
    CALL execute(.FALSE.)

    ! The physical input belongs to one pair; reject it for single-layer runs.
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_EXCHANGE_PARAMETERS'
    WRITE(unit,'(A)') " interlayer_exchange%elutriation_flag = T,"
    WRITE(unit,'(A)') " interlayer_exchange%elutriation_solid_volume_fraction = 0.15,"
    WRITE(unit,'(A)') " interlayer_exchange%elutriation_carrier_mode = 'ALLOW_CONCENTRATED' "
    WRITE(unit,'(A)') '/'
    n_layers=1
    CALL read_interlayer_exchange_parameters(unit,error)
    CALL require(INDEX(error,'requires N_LAYERS=2')>0,'reject interlayer input in single-layer')
    n_layers=2
    CALL read_interlayer_exchange_parameters(unit,error)
    CALL require(LEN_TRIM(error)==0,'read pair-wide interlayer input: '//TRIM(error))
    CALL close(interlayer_exchange%elutriation_solid_volume_fraction,0.15_wp,'read lifted solid fraction')
    CALL require(interlayer_exchange%elutriation_carrier_mode=='ALLOW_CONCENTRATED', &
        'read and normalize carrier mode')
    CLOSE(unit)

    PRINT *, 'PASS paired interlayer inventory, packing, co-lift modes and input validation'

CONTAINS

    SUBROUTINE execute(expected_valid)
        LOGICAL, INTENT(IN) :: expected_valid
        LOGICAL :: valid
        REAL(wp) :: original_l1(2),original_l2(2),new_l1(2),new_l2(2),new_c1,new_c2
        CALL limit_interlayer_mass_transfer(s1,s2,carrier1,carrier2, &
            settle,elut,dt,pack,config,accepted_s,accepted_e,down,up,valid)
        CALL require(valid.EQV.expected_valid,'limiter accepts/rejects expected case')
        IF (.NOT.valid) RETURN
        new_l1=s1+dt*(accepted_s-accepted_e)
        new_l2=s2+dt*(accepted_e-accepted_s)
        new_c1=carrier1+dt*(down-up)
        new_c2=carrier2+dt*(up-down)
        original_l1=s1; original_l2=s2
        CALL require(ALL(new_l1>=-1.0E-12_wp) .AND. ALL(new_l2>=-1.0E-12_wp), &
            'nonnegative species after both transfers')
        CALL require(new_c1>=-1.0E-12_wp .AND. new_c2>=-1.0E-12_wp, &
            'nonnegative carrier after both transfers')
        CALL require(SUM(new_l1)<=pack*(SUM(new_l1)+new_c1)+1.0E-12_wp, &
            'L1 final packing')
        CALL require(SUM(new_l2)<=pack*(SUM(new_l2)+new_c2)+1.0E-12_wp, &
            'L2 final packing')
        CALL require(ALL(ABS(new_l1+new_l2-original_l1-original_l2)<1.0E-12_wp), &
            'pairwise conservation for each solid species')
        CALL close(new_c1+new_c2,carrier1+carrier2,'pairwise carrier conservation')
        CALL require(ALL(dt*accepted_s<=s2+1.0E-12_wp) .AND. &
            ALL(dt*accepted_e<=s1+1.0E-12_wp), &
            'gross species withdrawals bounded by initial donors')
        CALL require(dt*up<=carrier1+1.0E-12_wp .AND. &
            dt*down<=carrier2+1.0E-12_wp,'gross carrier withdrawal bound')
    END SUBROUTINE execute

    SUBROUTINE close(actual,expected,message)
        REAL(wp), INTENT(IN) :: actual,expected
        CHARACTER(*), INTENT(IN) :: message
        CALL require(ABS(actual-expected)<=1.0E-11_wp*MAX(1.0_wp,ABS(expected)),message)
    END SUBROUTINE close

    SUBROUTINE require(condition,message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(*), INTENT(IN) :: message
        IF (.NOT.condition) THEN
            PRINT *, 'FAIL interlayer: ',TRIM(message)
            ERROR STOP 1
        END IF
    END SUBROUTINE require
END PROGRAM test_interlayer_exchange
