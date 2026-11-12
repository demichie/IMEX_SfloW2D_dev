PROGRAM test_layer_physics
    USE parameters_2d, ONLY : wp
    USE rheology_parameters_2d, ONLY : rheology_config_type
    USE layer_physics_2d, ONLY : layer_physics_type, LAYER_PHASE_GAS, &
        LAYER_PHASE_LIQUID, LAYER_PHASE_LAVA
    IMPLICIT NONE
    TYPE(layer_physics_type) :: l1,l2
    TYPE(rheology_config_type) :: config1,config2,returned
    CHARACTER(LEN=256) :: err
    config1%rheology_model='voellmy'
    config1%mu=0.3_wp
    config1%xi=300.0_wp
    CALL l1%initialize(LAYER_PHASE_LIQUID,.TRUE.,config1,err)
    CALL require(LEN_TRIM(err)==0,'liquid voellmy accepted')
    CALL require(l1%is_initialized(),'liquid initialized')
    CALL require(l1%phase()==LAYER_PHASE_LIQUID,'liquid phase retained')
    CALL require(l1%selected_rheology()=='rheology_voellmy','canonical rheology')
    CALL require(l1%selected_rheology_id()>0,'cached dispatch id')
    CALL l2%initialize(LAYER_PHASE_GAS,.FALSE.,config2,err)
    CALL require(LEN_TRIM(err)==0,'independent gas initialized')
    CALL require(l2%phase()==LAYER_PHASE_GAS,'second phase retained')
    CALL require(l2%selected_rheology_id()==0,'none dispatch id')
    CALL require(l1%selected_rheology_id()>0,'first layer unchanged')
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-0.3_wp)<EPSILON(1.0_wp),'layer retains L1 coefficient')
    returned%mu=5.0_wp
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-0.3_wp)<EPSILON(1.0_wp),'snapshot is read-only')
    config1%mu=7.0_wp
    returned=l1%rheology_configuration()
    CALL require(ABS(returned%mu-0.3_wp)<EPSILON(1.0_wp),'layer independent of input variable')
    config2%rheology_model='lava'
    config2%nu_ref=0.02_wp
    CALL l2%initialize(LAYER_PHASE_LAVA,.TRUE.,config2,err)
    CALL require(LEN_TRIM(err)==0,'lava rheology allowed')
    CALL l2%initialize(LAYER_PHASE_GAS,.TRUE.,config2,err)
    CALL require(LEN_TRIM(err)>0,'lava rheology rejected for gas')
    CALL require(.NOT.l2%is_initialized(),'invalid reconfiguration clears state')
    CALL l2%initialize(0,.FALSE.,config2,err)
    CALL require(LEN_TRIM(err)>0,'invalid phase rejected')
    config2%rheology_model='unknown'
    CALL l2%initialize(LAYER_PHASE_GAS,.TRUE.,config2,err)
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
