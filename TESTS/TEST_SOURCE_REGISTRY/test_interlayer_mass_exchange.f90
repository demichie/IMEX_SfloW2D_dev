!> @brief Cell-local paired transfer budgets and co-lift input regressions.
PROGRAM test_interlayer_mass_exchange
    USE parameters_2d, ONLY : wp,n_layers,n_solid,n_add_gas,n_vars,n_eqns
    USE parameters_2d, ONLY : n_stoch_vars,n_pore_vars,gas_flag,sutherland_flag
    USE parameters_2d, ONLY : configure_canonical_state_indices,idx_solidEqn_first
    USE parameters_2d, ONLY : idx_totMassEqn,idx_uEqn,idx_vEqn,idx_engyEqn,idx_addGasEqn_first
    USE interlayer_mass_exchange_parameters_2d, ONLY : interlayer_mass_exchange_config_type, interlayer_mass_exchange
    USE interlayer_mass_exchange_2d, ONLY : limit_interlayer_mass_transfer,eval_interlayer_mass_exchange_explicit
    USE interlayer_elutriation_2d, ONLY : eval_interlayer_elutriation_explicit
    USE interlayer_settling_2d, ONLY : eval_interlayer_settling_explicit
    USE state_types, ONLY : physical_state_real,initialize_physical_state_real
    USE solid_transport_parameters_2d, ONLY : rho_s,sp_heat_s,diam_s,maximum_solid_packing
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE mass_exchange_parameters_2d, ONLY : mass_exchange_config_type
    USE constitutive_parameters_2d, ONLY : grav,kin_visc_c,rho_a_amb,vonK
    USE source_input_2d, ONLY : read_interlayer_parameters,write_interlayer_parameters
    USE interlayer_momentum_exchange_2d, ONLY : interlayer_momentum_exchange
    IMPLICIT NONE
    TYPE(interlayer_mass_exchange_config_type) :: config
    REAL(wp) :: s1(2),s2(2),settle(2),elut(2),accepted_s(2),accepted_e(2)
    REAL(wp) :: down,up,dt,carrier1,carrier2,pack,strict_lift,concentrated_lift
    CHARACTER(LEN=256) :: error
    INTEGER :: unit,unit_read_status
    CHARACTER(LEN=256) :: first_line

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

    ! An omitted lifted fraction takes the current donor concentration.
    ! Test successive L1 inventories without changing the physical solid rate.
    config%elutriation_solid_volume_fraction=-1.0_wp
    s1=[0.12_wp,0.08_wp]; carrier1=0.8_wp
    CALL execute(.TRUE.)
    CALL require(SUM(accepted_e)>0.0_wp,'dynamic L1 lift remains active')
    CALL close(up,4.0_wp*SUM(accepted_e),'dynamic 20 percent L1 co-lift')
    carrier1=1.8_wp
    CALL execute(.TRUE.)
    CALL close(up,9.0_wp*SUM(accepted_e),'dynamic 10 percent L1 co-lift')
    ! A carrier-only or fully dry L1 cannot lift solids or create gas.
    s1=0.0_wp; carrier1=0.5_wp
    CALL execute(.TRUE.)
    CALL require(ALL(accepted_e==0.0_wp) .AND. up==0.0_wp, &
        'dynamic co-lift with no donor solids is a no-op')
    s1=0.0_wp; carrier1=0.0_wp
    elut=0.0_wp
    CALL execute(.TRUE.)
    CALL require(up==0.0_wp,'dynamic fully dry donor is a no-op')

    ! At limited carrier inventory, LIMIT_SOLIDS preserves parcel composition
    ! and also preserves the donor packing constraint.
    config%elutriation_solid_volume_fraction=0.1_wp
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
    config%elutriation_solid_volume_fraction=-0.5_wp
    CALL execute(.FALSE.)
    config%elutriation_solid_volume_fraction=1.1_wp
    CALL execute(.FALSE.)
    config%elutriation_solid_volume_fraction=0.0_wp
    CALL execute(.FALSE.)
    CALL require(ALL(accepted_s==0.0_wp) .AND. ALL(accepted_e==0.0_wp), &
        'invalid configuration resets transfer rates')
    config%elutriation_solid_volume_fraction=0.1_wp
    config%elutriation_carrier_mode='WRONG'
    CALL execute(.FALSE.)

    ! The physical input belongs to one pair; reject it for single-layer runs.
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') " interlayer_mass_exchange%elutriation_flag = T,"
    WRITE(unit,'(A)') " interlayer_mass_exchange%elutriation_solid_volume_fraction = 0.15,"
    WRITE(unit,'(A)') " interlayer_mass_exchange%elutriation_carrier_mode = 'ALLOW_CONCENTRATED' "
    WRITE(unit,'(A)') '/'
    n_layers=1
    CALL read_interlayer_parameters(unit,error)
    CALL require(INDEX(error,'requires N_LAYERS=2')>0,'reject interlayer input in single-layer')
    n_layers=2
    CALL read_interlayer_parameters(unit,error)
    CALL require(LEN_TRIM(error)==0,'read pair-wide interlayer input: '//TRIM(error))
    CALL close(interlayer_mass_exchange%elutriation_solid_volume_fraction,0.15_wp,'read lifted solid fraction')
    CALL require(interlayer_mass_exchange%elutriation_carrier_mode=='ALLOW_CONCENTRATED', &
        'read and normalize carrier mode')
    CLOSE(unit)

    ! NAMELIST omission retains the canonical -1 sentinel for local L1 co-lift.
    interlayer_mass_exchange=interlayer_mass_exchange_config_type()
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') ' interlayer_mass_exchange%elutriation_flag = T'
    WRITE(unit,'(A)') '/'
    CALL read_interlayer_parameters(unit,error)
    CALL require(LEN_TRIM(error)==0,'omitted co-lift fraction is accepted: '//TRIM(error))
    CALL close(interlayer_mass_exchange%elutriation_solid_volume_fraction,-1.0_wp, &
        'omitted co-lift fraction selects local L1 state')
    CLOSE(unit)

    ! An arbitrary negative value must not be interpreted as the sentinel.
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') ' interlayer_mass_exchange%elutriation_solid_volume_fraction = -0.5'
    WRITE(unit,'(A)') '/'
    CALL read_interlayer_parameters(unit,error)
    CALL require(INDEX(error,'ELUTRIATION_SOLID_VOLUME_FRACTION')>0, &
        'reject negative non-sentinel co-lift input')
    CLOSE(unit)
    interlayer_mass_exchange=interlayer_mass_exchange_config_type()

    ! The canonical pair-wide NAMELIST schedules mass and momentum separately.
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') " interlayer_mass_exchange%settling_flag = T,"
    WRITE(unit,'(A)') " interlayer_mass_exchange%method = 'EXPLICIT',"
    WRITE(unit,'(A)') ' interlayer_mass_exchange%group = 4,'
    WRITE(unit,'(A)') " interlayer_momentum_exchange%model = 'chezy',"
    WRITE(unit,'(A)') ' interlayer_momentum_exchange%chezy_drag_coeff = 0.03,'
    WRITE(unit,'(A)') " interlayer_momentum_exchange%method = 'analytic',"
    WRITE(unit,'(A)') ' interlayer_momentum_exchange%group = 5'
    WRITE(unit,'(A)') '/'
    CALL read_interlayer_parameters(unit,error)
    CALL require(LEN_TRIM(error)==0,'canonical interlayer input: '//TRIM(error))
    CALL require(interlayer_mass_exchange%settling_flag .AND. interlayer_mass_exchange%group==4, &
        'canonical input modifies the only mass exchange object')
    CALL require(interlayer_momentum_exchange%group==5 .AND. &
        interlayer_momentum_exchange%method=='ANALYTIC' .AND. &
        interlayer_momentum_exchange%model=='CHEZY', 'independent analytic drag configuration')
    CALL close(interlayer_momentum_exchange%chezy_drag_coeff,0.03_wp,'Chezy coefficient')

    ! The backup always emits the unified canonical group, not both spellings.
    CLOSE(unit)
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    CALL write_interlayer_parameters(unit)
    REWIND(unit)
    READ(unit,'(A)',IOSTAT=unit_read_status) first_line
    CALL require(unit_read_status==0 .AND. INDEX(first_line,'&INTERLAYER_PARAMETERS')>0, &
        'canonical interlayer backup group')
    CLOSE(unit)

    ! Multiple occurrences of the one canonical pair-wide group are invalid.
    OPEN(NEWUNIT=unit,STATUS='scratch',FORM='formatted')
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') '/'
    WRITE(unit,'(A)') '&INTERLAYER_PARAMETERS'
    WRITE(unit,'(A)') '/'
    CALL read_interlayer_parameters(unit,error)
    CALL require(INDEX(error,'Duplicate INTERLAYER_PARAMETERS')>0, &
        'reject duplicate canonical interlayer groups')
    CLOSE(unit)
    CALL test_density_asymmetry()
    CALL test_gross_carrier_exchange()
    CALL test_independent_process_limits()
    CALL test_concentrated_carrier_limit()
    CALL test_physical_sources()
    CALL test_bidirectional_physical_budget()

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

    !> @brief Verify both gross carrier streams with unequal donor and recipient densities.
    !> @details Donor caps apply to withdrawals before the opposite stream arrives.
    !>          Recipient packing uses transferred MASS divided by recipient density,
    !>          not a fictitious shared carrier volume or a net gross flux.
    SUBROUTINE test_gross_carrier_exchange()
        TYPE(interlayer_mass_exchange_config_type) :: pair
        REAL(wp) :: solids1(2),solids2(2),rates_down(2),rates_up(2)
        REAL(wp) :: accepted_down(2),accepted_up(2),flow_down,flow_up
        REAL(wp) :: carriers1,carriers2,dt_local,packing_local,alpha
        REAL(wp) :: donor_rho1,donor_rho2,receiver_rho1,receiver_rho2
        REAL(wp) :: new_solids1(2),new_solids2(2),new_carrier1,new_carrier2
        REAL(wp) :: transferred_carrier_mass
        REAL(wp) :: strict_solid,concentrated_solid
        LOGICAL :: valid
        INTEGER :: mode,fraction_mode

        pair%settling_flag=.TRUE.
        pair%elutriation_flag=.TRUE.
        dt_local=1.0_wp
        packing_local=0.6_wp
        solids1=[0.06_wp,0.03_wp]
        solids2=[0.04_wp,0.02_wp]
        carriers1=0.12_wp
        carriers2=0.16_wp
        rates_down=[0.010_wp,0.005_wp]
        rates_up=[0.018_wp,0.008_wp]
        donor_rho1=1.7_wp
        donor_rho2=1.2_wp
        receiver_rho1=1.5_wp
        receiver_rho2=1.4_wp

        ! Exercise both a prescribed parcel and the local L1 donor mixture
        ! with two opposing solid streams and four distinct carrier densities.
        DO fraction_mode=1,2
            IF (fraction_mode==1) THEN
                alpha=0.1_wp
                pair%elutriation_solid_volume_fraction=alpha
            ELSE
                alpha=SUM(solids1)/(SUM(solids1)+carriers1)
                pair%elutriation_solid_volume_fraction=-1.0_wp
            END IF
            DO mode=1,2
                IF (mode==1) THEN
                    pair%elutriation_carrier_mode='LIMIT_SOLIDS'
                ELSE
                    pair%elutriation_carrier_mode='ALLOW_CONCENTRATED'
                END IF

                CALL limit_interlayer_mass_transfer(solids1,solids2,carriers1,carriers2, &
                    rates_down,rates_up,dt_local,packing_local,pair,accepted_down,accepted_up, &
                    flow_down,flow_up,valid, &
                    rho_carrier_l1=donor_rho1,rho_carrier_l2=donor_rho2, &
                    candidate_rho_l1=receiver_rho1,candidate_rho_l2=receiver_rho2)
                CALL require(valid,'two gross carrier streams with distinct densities accepted')
                CALL require(SUM(accepted_down)>0.0_wp .AND. SUM(accepted_up)>0.0_wp, &
                    'both settling and elutriation remain active with each parcel mode')
                CALL require(ALL(accepted_down>=0.0_wp) .AND. ALL(accepted_up>=0.0_wp), &
                    'nonnegative per-species gross solid rates')
                CALL require(flow_down>=0.0_wp .AND. flow_up>=0.0_wp, &
                    'nonnegative independent gross carrier streams')
                CALL require(ALL(dt_local*accepted_down<=solids2+1.0E-12_wp), &
                    'settling cannot spend simultaneous elutriation receipts')
                CALL require(ALL(dt_local*accepted_up<=solids1+1.0E-12_wp), &
                    'elutriation cannot spend simultaneous settling receipts')
                CALL require(dt_local*flow_up<=carriers1+1.0E-12_wp .AND. &
                    dt_local*flow_down<=carriers2+1.0E-12_wp, &
                    'donor carrier gross withdrawals do not spend simultaneous receipts')

                new_solids1=solids1+dt_local*(accepted_down-accepted_up)
                new_solids2=solids2-dt_local*(accepted_down-accepted_up)
                new_carrier1=carriers1+dt_local* &
                    (donor_rho2*flow_down-donor_rho1*flow_up)/receiver_rho1
                new_carrier2=carriers2+dt_local* &
                    (donor_rho1*flow_up-donor_rho2*flow_down)/receiver_rho2
                transferred_carrier_mass=dt_local*(donor_rho2*flow_down-donor_rho1*flow_up)

                CALL require(ALL(new_solids1>=-1.0E-12_wp) .AND. ALL(new_solids2>=-1.0E-12_wp), &
                    'species remain admissible with two gross streams')
                CALL require(new_carrier1>=-1.0E-12_wp .AND. new_carrier2>=-1.0E-12_wp, &
                    'recipient carrier inventories remain nonnegative')
                CALL require(SUM(new_solids1)<=packing_local*(SUM(new_solids1)+new_carrier1)+1.0E-12_wp, &
                    'L1 packing with recipient carrier density')
                CALL require(SUM(new_solids2)<=packing_local*(SUM(new_solids2)+new_carrier2)+1.0E-12_wp, &
                    'L2 packing with recipient carrier density')
                CALL close(SUM(new_solids1+new_solids2),SUM(solids1+solids2), &
                    'pair conservation of species volume under dissimilar carrier densities')
                CALL close(receiver_rho1*(new_carrier1-carriers1),transferred_carrier_mass, &
                    'L1 receives actual gross carrier mass')
                CALL close(receiver_rho2*(new_carrier2-carriers2),-transferred_carrier_mass, &
                    'L2 receives exactly the opposite carrier mass')
                CALL close(receiver_rho1*(new_carrier1-carriers1)+ &
                    receiver_rho2*(new_carrier2-carriers2),0.0_wp, &
                    'carrier mass exchange is antisymmetric')

                IF (mode==1) THEN
                    strict_solid=SUM(accepted_up)
                    CALL close(flow_up,(1.0_wp-alpha)/alpha*strict_solid, &
                        'LIMIT_SOLIDS preserves upward parcel composition')
                ELSE
                    concentrated_solid=SUM(accepted_up)
                    CALL require(concentrated_solid>=strict_solid-1.0E-12_wp, &
                        'ALLOW_CONCENTRATED does not forbid otherwise admissible solid lifting')
                    CALL require(flow_up<=(1.0_wp-alpha)/alpha*concentrated_solid+1.0E-12_wp, &
                        'concentrated carrier lift does not exceed prescribed composition')
                END IF
            END DO
        END DO

    END SUBROUTINE test_gross_carrier_exchange

    !> @brief Require independent settling and elutriation reductions at a packing boundary.
    !> @details The L2 packing inequality limits the upward rate, while
    !>          full downward settling remains admissible. A common reduction
    !>          factor would unnecessarily suppress the descending flux.
    SUBROUTINE test_independent_process_limits()
        TYPE(interlayer_mass_exchange_config_type) :: pair
        REAL(wp) :: solids1(2),solids2(2),settling_rate(2),elutriation_rate(2)
        REAL(wp) :: accepted_down(2),accepted_up(2),flow_down,flow_up
        REAL(wp) :: new_l1(2),new_l2(2),carrier_l1,carrier_l2
        REAL(wp) :: down_factor,up_factor
        LOGICAL :: valid

        pair%settling_flag=.TRUE.
        pair%elutriation_flag=.TRUE.
        pair%elutriation_solid_volume_fraction=1.0_wp
        pair%elutriation_carrier_mode='LIMIT_SOLIDS'
        solids1=[0.12_wp,0.08_wp]
        solids2=[0.05_wp,0.03_wp]
        carrier_l1=0.4_wp
        carrier_l2=0.12_wp
        settling_rate=[0.006_wp,0.004_wp]
        elutriation_rate=[0.072_wp,0.048_wp]

        CALL limit_interlayer_mass_transfer(solids1,solids2,carrier_l1,carrier_l2, &
            settling_rate,elutriation_rate,1.0_wp,0.6_wp,pair, &
            accepted_down,accepted_up,flow_down,flow_up,valid)
        CALL require(valid,'packing-limited unequal reductions accepted')
        down_factor=SUM(accepted_down)/SUM(settling_rate)
        up_factor=SUM(accepted_up)/SUM(elutriation_rate)
        CALL require(down_factor>up_factor+1.0E-4_wp, &
            'settling and elutriation must not share a needless common limiter factor')
        CALL close(down_factor,1.0_wp,'admissible settling remains at its proposed physical rate')
        CALL require(up_factor<1.0_wp .AND. up_factor>0.0_wp, &
            'L2 packing reduces only the competing upward solid stream')
        CALL close(flow_up,0.0_wp,'pure-solid lift requires no upward carrier')
        CALL close(flow_down,0.0_wp,'L1 has adequate carrier for accepted settling')

        new_l1=solids1+accepted_down-accepted_up
        new_l2=solids2+accepted_up-accepted_down
        CALL require(ALL(new_l1>=0.0_wp) .AND. ALL(new_l2>=0.0_wp), &
            'independently limited solid inventories are nonnegative')
        CALL require(SUM(new_l1)<=0.6_wp*(SUM(new_l1)+carrier_l1)+1.0E-12_wp, &
            'L1 final packing with separately limited streams')
        CALL require(SUM(new_l2)<=0.6_wp*(SUM(new_l2)+carrier_l2)+1.0E-12_wp, &
            'L2 final packing with separately limited streams')
        CALL close(SUM(new_l1+new_l2),SUM(solids1+solids2), &
            'independently limited gross solid exchange is conservative')
    END SUBROUTINE test_independent_process_limits

    !> @brief Keep both solid laws active while reducing only upward carrier.
    !> @details In ALLOW_CONCENTRATED a packed donor can accept full gross
    !>          solid transfer even if a full co-lift would exhaust the
    !>          independent downward-carrier donor.
    SUBROUTINE test_concentrated_carrier_limit()
        TYPE(interlayer_mass_exchange_config_type) :: pair
        REAL(wp) :: lower(2),upper(2),settle_rate(2),elut_rate(2)
        REAL(wp) :: down_rate(2),up_rate(2),carrier_down,carrier_up
        REAL(wp) :: lower_after,upper_after,lower_carrier,upper_carrier
        LOGICAL :: valid

        pair%settling_flag=.TRUE.
        pair%elutriation_flag=.TRUE.
        pair%elutriation_solid_volume_fraction=0.2_wp
        pair%elutriation_carrier_mode='ALLOW_CONCENTRATED'
        lower=[0.12_wp,0.08_wp]
        upper=[0.05_wp,0.03_wp]
        lower_carrier=0.14_wp
        upper_carrier=0.12_wp
        settle_rate=[0.05_wp,0.03_wp]
        elut_rate=[0.03_wp,0.02_wp]

        CALL limit_interlayer_mass_transfer(lower,upper,lower_carrier,upper_carrier, &
            settle_rate,elut_rate,1.0_wp,0.6_wp,pair, &
            down_rate,up_rate,carrier_down,carrier_up,valid)
        CALL require(valid,'independent concentrated carrier bound accepted')
        CALL close(SUM(down_rate),SUM(settle_rate),'all descending solids accepted')
        CALL close(SUM(up_rate),SUM(elut_rate),'all ascending solids accepted')
        CALL require(carrier_up<lower_carrier-1.0E-4_wp .AND. carrier_up>0.0_wp, &
            'upward carrier independently reduced without rejecting solids')
        CALL require(carrier_up<4.0_wp*SUM(up_rate), &
            'ALLOW_CONCENTRATED permits a more concentrated lifted parcel')
        CALL require(carrier_down<=upper_carrier+1.0E-12_wp, &
            'downward carrier cannot spend simultaneous upward receipts')
        lower_after=SUM(lower)+SUM(down_rate)-SUM(up_rate)
        upper_after=SUM(upper)-SUM(down_rate)+SUM(up_rate)
        CALL require(lower_after<=0.6_wp*(lower_after+lower_carrier+carrier_down-carrier_up)+1.0E-12_wp, &
            'L1 packing with separately limited upward carrier')
        CALL require(upper_after<=0.6_wp*(upper_after+upper_carrier+carrier_up-carrier_down)+1.0E-12_wp, &
            'L2 packing with separately limited upward carrier')
    END SUBROUTINE test_concentrated_carrier_limit

    !> @brief Exercise physics and antisymmetric paired mass, momentum and energy.
    SUBROUTINE test_physical_sources()
        TYPE(physical_state_real) :: l1,l2,l1_recipient
        TYPE(interlayer_mass_exchange_config_type) :: phys_config
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
        CALL initialize_physical_state_real(l1_recipient)
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
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'cold-start L2 paired physical transfer accepted')
        CALL require(source1(idx_solidEqn_first)<0.0_wp .AND. &
            source2(idx_solidEqn_first)>0.0_wp,'elutriation solid donor and recipient')
        CALL require(source1(1)<0.0_wp .AND. source2(1)>0.0_wp,'elutriation gas included')

        ! The physical assembly derives co-lift from the donor L1 composition
        ! when the input fraction is omitted, retaining pairwise conservation.
        phys_config%elutriation_solid_volume_fraction=-1.0_wp
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'dynamic cold-start L2 paired physical transfer accepted')
        CALL close(source2(1)-SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)), &
            4.0_wp*SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)/rho_s), &
            'physical donor L1 gas-to-solid volume ratio')
        CALL close(source2(4),l1%T*( &
            SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)*sp_heat_s)+ &
            (source2(1)-SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)))*l1%cp_c), &
            'dynamic co-lift transports donor sensible energy')

        ! Group-zero physical rates may use pre-spatial L1, but the co-lift
        ! composition must use the current post-spatial donor inventory.
        l1_recipient=l1
        l1_recipient%alphas=[0.08_wp,0.02_wp]
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1_recipient,l2,dt,source1,source2,valid)
        CALL require(valid,'post-spatial donor composition accepted')
        CALL close(source2(1)-SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)), &
            9.0_wp*SUM(source2(idx_solidEqn_first:idx_solidEqn_first+1)/rho_s), &
            'dynamic co-lift uses recipient not pre-spatial donor')
        phys_config%elutriation_solid_volume_fraction=0.1_wp

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
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
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
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'simultaneous source-rate assembly accepted')
        DO i=1,6
            CALL close(source1(i)+source2(i),0.0_wp,'simultaneous antisymmetric conservative rate')
        END DO

        ! Both donors can have different carrier densities, heat capacities,
        ! temperatures and velocities. Rates remain paired per conservative
        ! equation; an energy or momentum credit cannot use recipient values.
        l1%rho_c=1.5_wp
        l2%rho_c=0.8_wp
        l1%rho_m=SUM(rho_s*l1%alphas)+l1%rho_c*(1.0_wp-SUM(l1%alphas))
        l2%rho_m=SUM(rho_s*l2%alphas)+l2%rho_c*(1.0_wp-SUM(l2%alphas))
        l1%xs=rho_s*l1%alphas/l1%rho_m
        l2%xs=rho_s*l2%alphas/l2%rho_m
        l1%T=420.0_wp
        l2%T=280.0_wp
        l1%cp_c=1000.0_wp
        l2%cp_c=1200.0_wp
        l1%u=12.0_wp
        l2%u=-3.0_wp
        l1%v=2.0_wp
        l2%v=-1.0_wp
        CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
            l1,l2,l1,l2,dt,source1,source2,valid)
        CALL require(valid,'two-way exchange with unequal donor thermodynamics accepted')
        CALL require(ABS(source1(1))>0.0_wp,'unequal carrier density induces nonzero paired mass exchange')
        CALL require(ABS(source1(4))>0.0_wp,'unequal donor temperatures transfer sensible heat')
        DO i=1,6
            CALL close(source1(i)+source2(i),0.0_wp, &
                'unequal donor density and temperature preserve each conservative pair')
        END DO

        ! The additional gas species follows the donor carrier composition.
        BLOCK
            TYPE(physical_state_real) :: gas_l1,gas_l2
            REAL(wp) :: augmented1(7),augmented2(7)
            n_add_gas=1; n_vars=7; n_eqns=7
            CALL configure_canonical_state_indices()
            ! A physical-state workspace has an immutable component layout.
            ! Construct separate states after changing the test-only layout.
            CALL initialize_physical_state_real(gas_l1)
            CALL initialize_physical_state_real(gas_l2)
            gas_l1%h=1.0_wp; gas_l1%alphas=[0.15_wp,0.05_wp]
            gas_l1%rho_c=1.0_wp; gas_l1%rho_m=455.8_wp
            gas_l1%xs=rho_s*gas_l1%alphas/gas_l1%rho_m
            gas_l1%xg=[0.0005_wp]
            gas_l1%cp_c=1000.0_wp; gas_l1%T=500.0_wp; gas_l1%u=10.0_wp; gas_l1%v=1.0_wp
            gas_l2%h=0.0_wp; gas_l2%rho_c=1.0_wp; gas_l2%rho_m=1.0_wp
            gas_l2%alphas=0.0_wp; gas_l2%cp_c=1000.0_wp; gas_l2%T=300.0_wp
            phys_config%settling_flag=.FALSE.
            CALL eval_interlayer_mass_exchange_explicit(phys_config,rheology_l2,mass_l2, &
                gas_l1,gas_l2,gas_l1,gas_l2,dt,augmented1,augmented2,valid)
            CALL require(valid,'co-lift with additional gas species accepted')
            CALL require(augmented1(7)<0.0_wp .AND. augmented2(7)>0.0_wp, &
                'additional gas constituent leaves donor and enters recipient')
            DO i=1,7
                CALL close(augmented1(i)+augmented2(i),0.0_wp, &
                    'additional gas pairwise conservation')
            END DO
        END BLOCK
    END SUBROUTINE test_physical_sources

    !> @brief Verify full physical source budgets for simultaneous opposed streams.
    !> @details Independent law and limiter evaluations provide expected phase
    !>          fluxes. Two different donor gas densities and heat capacities
    !>          must produce conservative mass, momentum, energy and species
    !>          rates without requiring carrier-volume conservation.
    SUBROUTINE test_bidirectional_physical_budget()
        TYPE(physical_state_real) :: l1,l2
        TYPE(interlayer_mass_exchange_config_type) :: pair
        TYPE(mass_exchange_config_type) :: mass_l2
        TYPE(rheology_config_type) :: rheology_l2
        REAL(wp) :: source1(7),source2(7),expected(7)
        REAL(wp) :: proposed_down(2),proposed_up(2),accepted_down(2),accepted_up(2)
        REAL(wp) :: solid1(2),solid2(2),carrier1,carrier2
        REAL(wp) :: down,up,down_mass,up_mass,dt_local,alpha_lift
        REAL(wp) :: chi,interface_u,interface_v,carrier1_after,carrier2_after
        REAL(wp) :: solid1_after(2),solid2_after(2),volume_change
        LOGICAL :: valid
        INTEGER :: fraction_mode,carrier_mode,i

        ! The preceding physical tests allocate both shared solid species and
        ! the one optional gas species. Rebuild fresh cell-local states.
        n_layers=2
        n_solid=2
        n_add_gas=1
        n_vars=7
        n_eqns=7
        CALL configure_canonical_state_indices()
        CALL initialize_physical_state_real(l1)
        CALL initialize_physical_state_real(l2)
        l1%h=0.01_wp
        l1%alphas=[0.35_wp,0.24_wp]
        l1%rho_c=1.5_wp
        l1%rho_m=SUM(rho_s*l1%alphas)+l1%rho_c*(1.0_wp-SUM(l1%alphas))
        l1%xs=rho_s*l1%alphas/l1%rho_m
        l1%xg=[0.0001_wp]
        l1%cp_c=900.0_wp
        l1%T=410.0_wp
        l1%u=8.0_wp
        l1%v=0.5_wp
        l2%h=1.0_wp
        l2%alphas=[0.11_wp,0.07_wp]
        l2%rho_c=0.8_wp
        l2%rho_m=SUM(rho_s*l2%alphas)+l2%rho_c*(1.0_wp-SUM(l2%alphas))
        l2%xs=rho_s*l2%alphas/l2%rho_m
        l2%xg=[0.0002_wp]
        l2%cp_c=1100.0_wp
        l2%T=290.0_wp
        l2%u=-2.0_wp
        l2%v=-0.25_wp
        pair%settling_flag=.TRUE.
        pair%elutriation_flag=.TRUE.
        pair%rouse_critical=100.0_wp
        dt_local=0.5_wp
        solid1=l1%h*l1%alphas
        solid2=l2%h*l2%alphas
        carrier1=l1%h*(1.0_wp-SUM(l1%alphas))
        carrier2=l2%h*(1.0_wp-SUM(l2%alphas))
        CALL eval_interlayer_settling_explicit(pair,rheology_l2,mass_l2,l2,dt_local,proposed_down)
        CALL eval_interlayer_elutriation_explicit(pair,mass_l2,l1,l2,dt_local,proposed_up)
        CALL require(ALL(proposed_down>0.0_wp) .AND. ALL(proposed_up>0.0_wp), &
            'both physical solid laws propose positive rates')

        DO fraction_mode=1,2
            IF (fraction_mode==1) THEN
                alpha_lift=0.15_wp
                pair%elutriation_solid_volume_fraction=alpha_lift
            ELSE
                alpha_lift=SUM(solid1)/(SUM(solid1)+carrier1)
                pair%elutriation_solid_volume_fraction=-1.0_wp
            END IF
            DO carrier_mode=1,2
                IF (carrier_mode==1) THEN
                    pair%elutriation_carrier_mode='LIMIT_SOLIDS'
                ELSE
                    pair%elutriation_carrier_mode='ALLOW_CONCENTRATED'
                END IF
                CALL limit_interlayer_mass_transfer(solid1,solid2,carrier1,carrier2, &
                    proposed_down,proposed_up,dt_local,maximum_solid_packing,pair, &
                    accepted_down,accepted_up,down,up,valid, &
                    rho_carrier_l1=l1%rho_c,rho_carrier_l2=l2%rho_c, &
                    candidate_rho_l1=l1%rho_c,candidate_rho_l2=l2%rho_c)
                CALL require(valid,'bidirectional donor budgets accepted')
                CALL require(ALL(accepted_down>0.0_wp) .AND. ALL(accepted_up>0.0_wp), &
                    'both accepted solid rates remain positive')
                CALL require(down>0.0_wp .AND. up>0.0_wp, &
                    'both gross carrier streams remain positive')
                IF (carrier_mode==1) CALL close(up, &
                    (1.0_wp-alpha_lift)/alpha_lift*SUM(accepted_up), &
                    'strict carrier-to-solid lift ratio')

                CALL eval_interlayer_mass_exchange_explicit(pair,rheology_l2,mass_l2, &
                    l1,l2,l1,l2,dt_local,source1,source2,valid)
                CALL require(valid,'physical bidirectional source accepted')
                down_mass=l2%rho_c*down
                up_mass=l1%rho_c*up
                chi=pair%settling_velocity_weight+ &
                    (1.0_wp-pair%settling_velocity_weight)* &
                    pair%settling_plug_thickness/(l2%h+pair%settling_plug_thickness)
                interface_u=l1%u+chi*(l2%u-l1%u)
                interface_v=l1%v+chi*(l2%v-l1%v)
                expected=0.0_wp
                expected(idx_solidEqn_first:idx_solidEqn_first+1)= &
                    rho_s*(accepted_down-accepted_up)
                expected(idx_totMassEqn)=SUM(expected(idx_solidEqn_first:idx_solidEqn_first+1))+ &
                    down_mass-up_mass
                expected(idx_uEqn)=SUM(rho_s*accepted_down)*interface_u- &
                    SUM(rho_s*accepted_up)*l1%u+down_mass*l2%u-up_mass*l1%u
                expected(idx_vEqn)=SUM(rho_s*accepted_down)*interface_v- &
                    SUM(rho_s*accepted_up)*l1%v+down_mass*l2%v-up_mass*l1%v
                expected(idx_engyEqn)=SUM(rho_s*sp_heat_s*accepted_down)*l2%T- &
                    SUM(rho_s*sp_heat_s*accepted_up)*l1%T+ &
                    down_mass*l2%cp_c*l2%T-up_mass*l1%cp_c*l1%T
                expected(idx_addGasEqn_first)=down_mass*l2%xg(1)/(1.0_wp-SUM(l2%xs))- &
                    up_mass*l1%xg(1)/(1.0_wp-SUM(l1%xs))
                DO i=1,n_eqns
                    CALL close(source1(i),expected(i), &
                        'resolved physical source by conservative component')
                    CALL close(source1(i)+source2(i),0.0_wp, &
                        'opposed physical streams preserve pair budget')
                END DO

                solid1_after=solid1+dt_local*(accepted_down-accepted_up)
                solid2_after=solid2+dt_local*(accepted_up-accepted_down)
                carrier1_after=carrier1+dt_local*(down_mass-up_mass)/l1%rho_c
                carrier2_after=carrier2-dt_local*(down_mass-up_mass)/l2%rho_c
                CALL require(ALL(solid1_after>=0.0_wp) .AND. ALL(solid2_after>=0.0_wp), &
                    'solid species remain nonnegative')
                CALL require(carrier1_after>=0.0_wp .AND. carrier2_after>=0.0_wp, &
                    'gas inventories remain nonnegative')
                CALL require(SUM(solid1_after)<=maximum_solid_packing* &
                    (SUM(solid1_after)+carrier1_after)+1.0E-12_wp, &
                    'L1 packing after bidirectional source')
                CALL require(SUM(solid2_after)<=maximum_solid_packing* &
                    (SUM(solid2_after)+carrier2_after)+1.0E-12_wp, &
                    'L2 packing after bidirectional source')
                CALL close(l1%rho_c*(carrier1_after-carrier1)+ &
                    l2%rho_c*(carrier2_after-carrier2),0.0_wp, &
                    'two-temperature carrier mass conservation')
                volume_change=carrier1_after+carrier2_after-carrier1-carrier2
                CALL close(volume_change,dt_local*(down_mass-up_mass)* &
                    (1.0_wp/l1%rho_c-1.0_wp/l2%rho_c), &
                    'unequal carrier densities map mass into different volumes')
            END DO
        END DO
    END SUBROUTINE test_bidirectional_physical_budget

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
END PROGRAM test_interlayer_mass_exchange
