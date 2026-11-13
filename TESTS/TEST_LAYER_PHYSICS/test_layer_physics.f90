PROGRAM test_layer_physics
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE thermal_parameters_2d, ONLY : thermal_config_type
    USE pore_pressure_parameters_2d, ONLY : pore_pressure_config_type
    USE layer_physics_2d, ONLY : layer_physics_type, LAYER_PHASE_GAS, &
        LAYER_PHASE_LIQUID, LAYER_PHASE_LAVA
    IMPLICIT NONE
    TYPE(layer_physics_type) :: l1,l2
    TYPE(rheology_config_type), TARGET :: config1,config2
    TYPE(rheology_config_type) :: returned
    TYPE(thermal_config_type), TARGET :: thermal1,thermal2
    TYPE(thermal_config_type) :: thermal_view
    TYPE(pore_pressure_config_type), TARGET :: pore1,pore2
    TYPE(pore_pressure_config_type) :: pore_view
    CHARACTER(LEN=256) :: err
    thermal1%T_ground=300.0_wp
    thermal2%T_ground=300.0_wp
    thermal1%emissivity=0.8_wp
    thermal2%emissivity=0.3_wp
    config1%rheology_model='voellmy'
    config1%mu=0.3_wp
    config1%xi=300.0_wp
    CALL l1%initialize(LAYER_PHASE_LIQUID,.TRUE.,config1,thermal1,pore1,err)
    CALL require(LEN_TRIM(err)==0,'liquid voellmy accepted')
    CALL require(l1%is_initialized(),'liquid initialized')
    CALL require(l1%phase()==LAYER_PHASE_LIQUID,'liquid phase retained')
    CALL require(l1%selected_rheology()=='rheology_voellmy','canonical rheology')
    CALL require(l1%selected_rheology_id()>0,'cached dispatch id')
    CALL l2%initialize(LAYER_PHASE_GAS,.FALSE.,config2,thermal2,pore2,err)
    CALL require(LEN_TRIM(err)==0,'independent gas initialized')
    CALL require(l2%phase()==LAYER_PHASE_GAS,'second phase retained')
    CALL require(l2%selected_rheology_id()==0,'none dispatch id')
    CALL require(l1%selected_rheology_id()>0,'first layer unchanged')
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-0.3_wp)<EPSILON(1.0_wp),'layer retains L1 coefficient')
    returned%mu=5.0_wp
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-0.3_wp)<EPSILON(1.0_wp),'value view is read-only')
    config1%mu=7.0_wp
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-7.0_wp)<EPSILON(1.0_wp),'binding tracks canonical rheology input')
    thermal_view=l1%thermal_configuration()
    CALL require(thermal_view%emissivity==0.8_wp,'L1 reads canonical thermal input')
    thermal_view%emissivity=0.0_wp
    thermal_view=l1%thermal_configuration()
    CALL require(thermal_view%emissivity==0.8_wp,'view does not modify canonical thermal input')
    thermal1%emissivity=0.6_wp
    thermal_view=l1%thermal_configuration()
    CALL require(thermal_view%emissivity==0.6_wp,'binding tracks canonical thermal input')
    thermal_view=l2%thermal_configuration()
    CALL require(thermal_view%emissivity==0.3_wp,'L2 independent thermal configuration')
    pore1%hydraulic_permeability=1.0e-10_wp
    pore2%hydraulic_permeability=1.0e-12_wp
    pore_view=l1%pore_pressure_configuration()
    CALL require(pore_view%hydraulic_permeability==pore1%hydraulic_permeability, &
        'L1 Darcy configuration follows its canonical owner')
    pore_view%hydraulic_permeability=0.0_wp
    pore_view=l1%pore_pressure_configuration()
    CALL require(pore_view%hydraulic_permeability==pore1%hydraulic_permeability, &
        'Darcy value view does not modify the canonical owner')
    pore_view=l2%pore_pressure_configuration()
    CALL require(pore_view%hydraulic_permeability==pore2%hydraulic_permeability, &
        'Darcy L2 remains independently configured')
    config2%rheology_model='lava'
    config2%nu_ref=0.02_wp
    CALL l2%initialize(LAYER_PHASE_LAVA,.TRUE.,config2,thermal2,pore2,err)
    CALL require(LEN_TRIM(err)==0,'lava rheology allowed')
    CALL l2%initialize(LAYER_PHASE_GAS,.TRUE.,config2,thermal2,pore2,err)
    CALL require(LEN_TRIM(err)>0,'lava rheology rejected for gas')
    CALL require(.NOT.l2%is_initialized(),'invalid reconfiguration clears state')
    CALL l2%initialize(0,.FALSE.,config2,thermal2,pore2,err)
    CALL require(LEN_TRIM(err)>0,'invalid phase rejected')
    config2%rheology_model='unknown'
    CALL l2%initialize(LAYER_PHASE_GAS,.TRUE.,config2,thermal2,pore2,err)
    CALL require(LEN_TRIM(err)>0,'unknown rheology rejected')
    CALL l1%finalize()
    CALL require(.NOT.l1%is_initialized(),'finalize resets configuration')
    PRINT *, 'PASS: independent per-layer rheology/phase descriptors'
CONTAINS
    SUBROUTINE require(ok,label)
        LOGICAL, INTENT(IN) :: ok
        CHARACTER(*), INTENT(IN) :: label
        IF (.NOT.ok) THEN
            PRINT *, 'FAIL: ',label
            ERROR STOP 1
        END IF
    END SUBROUTINE require
END PROGRAM test_layer_physics
