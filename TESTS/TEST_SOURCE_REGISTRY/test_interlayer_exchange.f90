!> @brief Cell-local paired transfer budgets and co-lift input regressions.
PROGRAM test_interlayer_exchange
    USE parameters_2d, ONLY : wp,n_layers,n_solid,n_add_gas,n_vars,n_eqns
    USE parameters_2d, ONLY : n_stoch_vars,n_pore_vars,gas_flag,sutherland_flag
    USE parameters_2d, ONLY : configure_canonical_state_indices,idx_solidEqn_first
    USE interlayer_exchange_parameters_2d, ONLY : interlayer_exchange_config_type, interlayer_exchange
    USE interlayer_exchange_2d, ONLY : limit_interlayer_mass_transfer,eval_interlayer_exchange_explicit
    USE interlayer_elutriation_2d, ONLY : eval_interlayer_elutriation_explicit
    USE state_types, ONLY : physical_state_real,initialize_physical_state_real
    USE solid_transport_parameters_2d, ONLY : rho_s,sp_heat_s,diam_s,maximum_solid_packing
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE constitutive_parameters_2d, ONLY : grav,kin_visc_c,rho_a_amb,vonK
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
    CALL test_density_asymmetry()
    CALL test_physical_sources()

    PRINT *, 'PASS paired interlayer inventory, packing, co-lift modes and input validation'

CONTAINS


    !> @brief Verify carrier MASS conservation with distinct carrier densities.
    SUBROUTINE test_density_asymmetry()
        LOGICAL :: valid
        REAL(wp) :: sd(2),ed(2),cd,cu
        s1=0.0_wp
        s2=[0.12_wp,0.08_wp]
        carrier1=0.0_wp
        carrier2=0.5_wp
        settle=[0.03_wp,0.02_wp]
        elut=0.0_wp
        config%elutriation_carrier_mode='LIMIT_SOLIDS'
        CALL limit_interlayer_mass_transfer(s1,s2,carrier1,carrier2,settle,elut, &
            dt,pack,config,sd,ed,cd,cu,valid,rho_carrier_l1=2.0_wp,rho_carrier_l2=1.0_wp)
        CALL require(valid,'unequal carrier densities are admissible')
        CALL close(cd,2.0_wp*(1.0_wp-pack)/pack*SUM(sd),'unequal density cold-start packing')
        CALL close(2.0_wp*(carrier1+0.5_wp*dt*cd-dt*cu)+ &
            (carrier2+dt*(2.0_wp*cu-cd)),2.0_wp*carrier1+carrier2, &
            'unequal density carrier mass conservation')
    END SUBROUTINE test_density_asymmetry

    !> @brief Exercise physics and antisymmetric paired mass, momentum and energy.
    SUBROUTINE test_physical_sources()
        TYPE(physical_state_real) :: l1,l2
        TYPE(interlayer_exchange_config_type) :: phys_config
        TYPE(mass_exchange_config_type) :: mass_l2
        TYPE(rheology_config_type) :: rheology_l2
        REAL(wp) :: source1(6),source2(6),elut_rates(2)
        LOGICAL :: valid
        INTEGER :: i

        n_solid=2; n_add_gas=0; n_stoch_vars=0; n_pore_vars=0
        n_vars=6; n_eqns=6
        CALL configure_canonical_state_indices()
        gas_flag=.TRUE.
        sutherland_flag=.FALSE.
        rho_a_amb=1.0_wp; grav=9.81_wp; vonK=0.41_wp
        kin_visc_c=1.0E-5_wp; maximum_solid_packing=0.6_wp
        ALLOCATE(rho_s(2),sp_heat_s(2),diam_s(2))
        rho_s=[2200.0_wp,2500.0_wp]
        sp_heat_s=[800.0_wp,900.0_wp]
        diam_s=[2.0E-5_wp,5.0E-5_wp]
        CALL initialize_physical_state_real(l1)
        CALL initialize_physical_state_real(l2)
        l1%h=1.0_wp; l1%alphas=[0.15_wp,0.05_wp]
        l1%rho_c=1.0_wp; l1%rho_m=455.8_wp
        l1%xs=rho_s*l1%alphas/l1%rho_m
        l1%cp_c=1000.0_wp; l1%T=500.0_wp; l1%u=10.0_wp; l1%v=1.0_wp
        l2%h=0.0_wp; l2%alphas=0.0_wp
        l2%rho_c=1.0_wp; l2%rho_m=1.0_wp
        l2%cp_c=1000.0_wp; l2%T=300.0_wp
        phys_config%elutriation_flag=.TRUE.
        phys_config%elutriation_solid_volume_fraction=0.1_wp

        CALL eval_interlayer_elutriation_explicit(phys_config,mass_l2,l1,l2,dt,elut_rates)
        CALL require(ALL(elut_rates>0.0_wp),'cold-start elutriation generates both grain rates')
        phys_config%elutriation_law='SHEAR_POWER'
        phys_config%shields_critical=0.05_wp
        phys_config%shields_exponent=1.0_wp
        phys_config%lift_length_diameters=10.0_wp
        CALL eval_interlayer_elutriation_explicit(phys_config,mass_l2,l1,l2,dt,elut_rates)
        CALL require(ALL(elut_rates>0.0_wp),'shear-power and Shields rate in dry L2 cold-start')
        phys_config%elutriation_law='LINEAR_USTAR'
        CALL eval_interlayer_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'cold-start L2 paired physical transfer accepted')
        CALL require(source1(idx_solidEqn_first)<0.0_wp .AND. &
            source2(idx_solidEqn_first)>0.0_wp,'elutriation solid donor and recipient')
        CALL require(source1(1)<0.0_wp .AND. source2(1)>0.0_wp,'elutriation gas included')

        DO i=1,6
            CALL close(source1(i)+source2(i),0.0_wp,'elutriation antisymmetric conservative rate')
        END DO
        CALL close(source2(2),-source1(1)*l1%u,'elutriation donor velocity')
        CALL close(source2(3),-source1(1)*l1%v,'elutriation donor transverse velocity')

        ! Settling is independent of external basal deposition activation.
        phys_config%elutriation_flag=.FALSE.
        phys_config%settling_flag=.TRUE.
        l1%h=0.0_wp; l1%alphas=0.0_wp; l1%xs=0.0_wp
        l2%h=1.0_wp; l2%alphas=[0.10_wp,0.05_wp]
        l2%rho_m=1.0_wp+SUM((rho_s-1.0_wp)*l2%alphas)
        l2%xs=rho_s*l2%alphas/l2%rho_m
        l2%u=3.0_wp; l2%v=-1.0_wp
        CALL eval_interlayer_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'cold-start L1 settling accepted')
        CALL require(source1(idx_solidEqn_first)>0.0_wp,'settling creates L1 solids')
        CALL require(source1(1)>0.0_wp .AND. source2(1)<0.0_wp,'settling carrier is withdrawn')
        DO i=1,6
            CALL close(source1(i)+source2(i),0.0_wp,'settling antisymmetric conservative rate')
        END DO

        ! Both kernels can be evaluated at the same snapshot without state mutation.
        phys_config%elutriation_flag=.TRUE.
        l1%h=1.0_wp; l1%alphas=[0.15_wp,0.05_wp]
        l1%rho_m=1.0_wp+SUM((rho_s-1.0_wp)*l1%alphas)
        l1%xs=rho_s*l1%alphas/l1%rho_m
        CALL eval_interlayer_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'simultaneous source-rate assembly accepted')
        DO i=1,6
            CALL close(source1(i)+source2(i),0.0_wp,'simultaneous antisymmetric conservative rate')
        END DO

        ! The additional gas species follows the donor carrier composition.
        BLOCK
            REAL(wp) :: augmented1(7),augmented2(7)
            n_add_gas=1; n_vars=7; n_eqns=7
            CALL configure_canonical_state_indices()
            CALL initialize_physical_state_real(l1)
            CALL initialize_physical_state_real(l2)
            l1%h=1.0_wp; l1%alphas=[0.15_wp,0.05_wp]
            l1%rho_c=1.0_wp; l1%rho_m=455.8_wp
            l1%xs=rho_s*l1%alphas/l1%rho_m
            l1%xg=[0.0005_wp]
            l1%cp_c=1000.0_wp; l1%T=500.0_wp; l1%u=10.0_wp; l1%v=1.0_wp
            l2%h=0.0_wp; l2%rho_c=1.0_wp; l2%rho_m=1.0_wp
            l2%alphas=0.0_wp; l2%cp_c=1000.0_wp; l2%T=300.0_wp
            phys_config%settling_flag=.FALSE.
            CALL eval_interlayer_exchange_explicit(phys_config,rheology_l2,mass_l2, &
                l1,l2,l1,l2,dt,augmented1,augmented2,valid)
            CALL require(valid,'co-lift with additional gas species accepted')
            CALL require(augmented1(7)<0.0_wp .AND. augmented2(7)>0.0_wp, &
                'additional gas constituent leaves donor and enters recipient')
            DO i=1,7
                CALL close(augmented1(i)+augmented2(i),0.0_wp, &
                    'additional gas pairwise conservation')
            END DO
        END BLOCK
    END SUBROUTINE test_physical_sources

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
