PROGRAM test_source_registry
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  USE source_catalog_2d
  USE source_registry_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_NEWTON, IM_GROUP_METHOD_ANALYTIC
  USE source_configuration_2d, ONLY : initialize_configured_sources
  USE source_input_2d
  USE parameters_2d, ONLY : wp, n_eqns, idx_pore, &
       rheology_flag, rheology_model, pore_pressure_flag, thermal_cooling_flag, &
       curvature_term_flag, bottom_radial_source_flag, bottom_fissural_source_flag, &
       radial_source_flag, lateral_source_flag, liquid_flag, gas_flag, &
       gas_loss_flag, liquid_vaporization_flag
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       coeff_porosity, entrainment_flag
  IMPLICIT NONE
  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: registry,invalid
  TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:),wrong(:)
  TYPE(source_term_type) :: term
  CHARACTER(LEN=64), ALLOCATABLE :: active(:)
  CHARACTER(LEN=256) :: error
  INTEGER :: i,mode,a,b,c,bits,ids(3),methods(3),membership(7),nm
  LOGICAL :: mask(7),enabled(3)
  CALL require(SOURCE_TERM_COUNT==24,'physical catalog contains 24 indivisible sources including extension fixture')
  CALL require(source_catalog_names_unique(),'unique physical source names')
  CALL require(canonical_rheology_catalog_unique(),'unique complete rheology names')
  term=source_term_descriptor_by_name('rheology_newmodel')
  CALL require(term%has_newton .AND. .NOT.term%has_explicit .AND. &
       .NOT.term%has_analytic,'Newton-only extension capabilities')
  DO i=1,SOURCE_TERM_COUNT
    term=source_term_descriptor_by_name(SOURCE_CATALOG(i)%name)
    CALL require(term%name==SOURCE_CATALOG(i)%name,'physical catalog lookup')
    CALL require(LEN_TRIM(term%target)>0,'physical target declared')
  END DO
  term=source_term_descriptor_by_name('SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS')
  CALL require(term%name=='','legacy mathematical addend removed')
  term=source_term_descriptor_by_name('mass_exchange_carrier_deposition_loss')
  CALL require(term%name=='','dependent deposition carrier loss is not independently schedulable')
  term=source_term_descriptor_by_name('SOURCE_TERM_THERMAL_RADIATIVE_COOLING')
  CALL require(term%name=='','thermal formula component is not schedulable')
  DO mode=0,3
    enabled=[.TRUE.,BTEST(mode,0),BTEST(mode,1)]
    mask=.FALSE.;mask(2:3)=.TRUE.;mask(4)=enabled(2);mask(7)=enabled(3)
    active=[CHARACTER(LEN=64)::'rheology_lava','geometry_curvature', &
         'bottom_circular_injection','boundary_injection','mass_exchange_deposition']
    IF (enabled(2)) active=[active,'thermal_cooling']
    IF (enabled(3)) active=[active,'pore_pressure_dissipation']
    DO a=1,3
    DO b=1,3
    DO c=1,3
    DO bits=0,7
      ids=[a,b,c]
      DO i=1,3
        methods(i)=MERGE(IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON,BTEST(bits,i-1))
      END DO
      membership=0;nm=0
      BLOCK
        INTEGER :: selected(3)
        DO i=1,3
          IF (.NOT.ANY(enabled .AND. ids==i)) CYCLE
          nm=nm+1;selected(nm)=methods(i)
          IF (enabled(1) .AND. ids(1)==i) membership(2:3)=nm
          IF (enabled(2) .AND. ids(2)==i) membership(4)=nm
          IF (enabled(3) .AND. ids(3)==i) membership(7)=nm
        END DO
        CALL partition%initialize(mask,membership,selected(1:nm))
      END BLOCK
      CALL configure_lava_case()
      thermal_cooling_flag=enabled(2)
      pore_pressure_flag=enabled(3)
      rheology_group=ids(1); thermal_group=ids(2); pore_pressure_group=ids(3)
      rheology_method=MERGE('ANALYTIC','NEWTON  ',methods(ids(1))==IM_GROUP_METHOD_ANALYTIC)
      thermal_method=MERGE('ANALYTIC','NEWTON  ',methods(ids(2))==IM_GROUP_METHOD_ANALYTIC)
      pore_pressure_method=MERGE('ANALYTIC','NEWTON  ',methods(ids(3))==IM_GROUP_METHOD_ANALYTIC)
      CALL initialize_configured_sources(registry,partition,error)
      CALL require(error=='' .AND. registry%is_initialized(), &
           'physical triplet combination normalizes: '//TRIM(error))
      stages=registry%stage_groups();steps=registry%step_groups()
      CALL require(SIZE(steps)==1,'one post-RK source group')
      CALL require(group_contains_name(steps(1),'mass_exchange_deposition'),'mass source ownership')
      CALL require(.NOT.group_contains_name(steps(1),'rheology_lava'),'no stage source in step')
      CALL registry%validate(error)
      CALL require(error=='','repeat physical registry validation')
      wrong=stages
      DO i=1,SIZE(wrong)
        IF (wrong(i)%method==SOURCE_METHOD_EXPLICIT) CYCLE
        IF (.NOT.ANY(wrong(i)%equation_ids==2)) CYCLE
        wrong(i)%source_names=[CHARACTER(LEN=64)::'SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS']
        CALL invalid%initialize(7,MERGE(7,0,enabled(3)),wrong,steps,error)
        CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
             'mathematical addend rejected in velocity group')
        EXIT
      END DO
    END DO
    END DO
    END DO
    END DO
  END DO
  CALL verify_default_schedule()
  CALL verify_thermal_activation_matrix()
  CALL verify_analytic_step_schedule()
  CALL verify_newton_explicit_step_schedule()
  WRITE(*,'(A)') 'D18 physical catalog and 864 physical triplet mappings verified'
CONTAINS
  SUBROUTINE configure_lava_case()
    n_eqns=7;idx_pore=7
    rheology_flag=.TRUE.;rheology_model='lava'
    thermal_cooling_flag=.TRUE.
    pore_pressure_flag=.TRUE.
    curvature_term_flag=.FALSE.
    bottom_radial_source_flag=.FALSE.
    bottom_fissural_source_flag=.FALSE.
    radial_source_flag=.FALSE.
    lateral_source_flag=.TRUE.
    liquid_flag=.FALSE.;gas_flag=.FALSE.
    gas_loss_flag=.FALSE.;liquid_vaporization_flag=.FALSE.
    erosion_coeff=0.0_wp;settling_flag=.TRUE.
    coeff_porosity=0.0_wp;entrainment_flag=.FALSE.
    rheology_method='NEWTON'; rheology_cadence='STAGE'; rheology_group=3
    thermal_method='ANALYTIC'; thermal_cadence='STAGE'; thermal_group=1
    pore_pressure_method='ANALYTIC'; pore_pressure_cadence='STAGE'; pore_pressure_group=2
  END SUBROUTINE configure_lava_case

  SUBROUTINE verify_default_schedule()
    TYPE(source_group_type), ALLOCATABLE :: g(:),s(:)
    LOGICAL :: enabled_eq(7)
    INTEGER :: assigned(7)
    CALL configure_lava_case()
    enabled_eq=.FALSE.;enabled_eq(2:4)=.TRUE.;enabled_eq(7)=.TRUE.
    assigned=0;assigned(2:3)=3;assigned(4)=1;assigned(7)=2
    CALL partition%initialize(enabled_eq,assigned, &
         [IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON])
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','configured default physical schedule: '//TRIM(error))
    g=registry%stage_groups();s=registry%step_groups()
    CALL require(SIZE(g)==4 .AND. SIZE(s)==1,'default schedule group count')
    CALL require(SIZE(g(1)%source_names)==1 .AND. &
         group_contains_name(g(1),'thermal_cooling'),'one complete cooling source')
    CALL require(SIZE(g(2)%source_names)==1 .AND. &
         group_contains_name(g(2),'pore_pressure_dissipation'),'pore dissipation scheduling')
    CALL require(SIZE(g(3)%source_names)==1 .AND. &
         group_contains_name(g(3),'rheology_lava'),'complete lava law scheduling')
    CALL require(group_contains_name(g(4),'boundary_injection'), &
         'explicit boundary scheduling')
    CALL require(group_contains_name(s(1),'mass_exchange_deposition'), &
         'post-RK deposition scheduling')
  END SUBROUTINE verify_default_schedule

  SUBROUTINE verify_analytic_step_schedule()
    TYPE(source_group_type), ALLOCATABLE :: stage(:),post(:)
    CALL configure_lava_case()
    thermal_method='ANALYTIC'; thermal_cadence='STEP'; thermal_group=2
    pore_pressure_method='ANALYTIC'; pore_pressure_cadence='STEP'; pore_pressure_group=1
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','analytic STEP source configuration: '//TRIM(error))
    stage=registry%stage_groups(); post=registry%step_groups()
    CALL require(SIZE(post)==3,'two analytic STEP groups and one mass exchange batch')
    CALL require(post(1)%method==SOURCE_METHOD_ANALYTIC .AND. &
         group_contains_name(post(1),'pore_pressure_dissipation'), &
         'pore STEP first by GROUP')
    CALL require(post(2)%method==SOURCE_METHOD_ANALYTIC .AND. &
         group_contains_name(post(2),'thermal_cooling'), &
         'cooling STEP second by GROUP')
    CALL require(post(3)%method==SOURCE_METHOD_EXPLICIT .AND. &
         group_contains_name(post(3),'mass_exchange_deposition'), &
         'mass exchange remains atomic and last')
    CALL require(partition%n_IM_groups==1, &
         'STEP analytic sources do not enter stage implicit partition')
  END SUBROUTINE verify_analytic_step_schedule

  SUBROUTINE verify_newton_explicit_step_schedule()
    TYPE(source_group_type), ALLOCATABLE :: g(:),post(:)
    CALL configure_lava_case()
    rheology_method='NEWTON'; rheology_cadence='STEP'; rheology_group=1
    thermal_method='EXPLICIT'; thermal_cadence='STEP'; thermal_group=2
    pore_pressure_method='ANALYTIC'; pore_pressure_cadence='STAGE'; pore_pressure_group=1
    CALL initialize_configured_sources(registry,partition,error)
    CALL require(error=='','Newton and explicit STEP scheduling: '//TRIM(error))
    g=registry%stage_groups();post=registry%step_groups()
    CALL require(partition%n_IM_groups==1,'only pore pressure is implicit during RK')
    CALL require(SIZE(post)==3,'Newton, explicit and mass-exchange STEP groups')
    CALL require(post(1)%method==SOURCE_METHOD_NEWTON .AND. &
         group_contains_name(post(1),'rheology_lava'), &
         'Newton STEP complete rheology group')
    CALL require(post(2)%method==SOURCE_METHOD_EXPLICIT .AND. &
         group_contains_name(post(2),'thermal_cooling'), &
         'explicit STEP thermal source')
    CALL require(group_contains_name(post(3),'mass_exchange_deposition'), &
         'atomic mass-exchange STEP group preserved')
  END SUBROUTINE verify_newton_explicit_step_schedule

  SUBROUTINE verify_thermal_activation_matrix()
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    LOGICAL :: enabled_eq(7),has_cooling
    INTEGER :: assigned(7),model,choice,method(2),i
    DO model=1,2
      DO choice=1,2
        CALL configure_lava_case()
        rheology_model=MERGE('lava   ','voellmy',model==2)
        thermal_cooling_flag=(choice==2)
        pore_pressure_flag=.FALSE.
        lateral_source_flag=.FALSE.
        settling_flag=.FALSE.
        enabled_eq=.FALSE.; enabled_eq(2:3)=.TRUE.
        enabled_eq(4)=thermal_cooling_flag
        assigned=0
        IF (thermal_cooling_flag) THEN
          assigned(4)=1
          assigned(2:3)=2
          method=[IM_GROUP_METHOD_ANALYTIC,IM_GROUP_METHOD_NEWTON]
          CALL partition%initialize(enabled_eq,assigned,method)
        ELSE
          assigned(2:3)=1
          CALL partition%initialize(enabled_eq,assigned,[IM_GROUP_METHOD_NEWTON])
        END IF
        CALL initialize_configured_sources(registry,partition,error)
        CALL require(error=='','cooling on/off independent of rheology: '//TRIM(error))
        groups=registry%stage_groups()
        has_cooling=.FALSE.
        DO i=1,SIZE(groups)
          has_cooling=has_cooling .OR. group_contains_name(groups(i),'thermal_cooling')
        END DO
        CALL require(has_cooling .EQV. thermal_cooling_flag, &
             'physical cooling membership follows dedicated flag')
        CALL require(COUNT(enabled_eq)==MERGE(3,2,thermal_cooling_flag), &
             'legacy equation layout changes only when cooling active')
      END DO
    END DO
  END SUBROUTINE verify_thermal_activation_matrix

  SUBROUTINE require(valid,message)
    LOGICAL, INTENT(IN) :: valid
    CHARACTER(*), INTENT(IN) :: message
    IF (valid) RETURN
    WRITE(error_unit,'(A)') 'FAIL: '//message
    ERROR STOP 1
  END SUBROUTINE require
END PROGRAM test_source_registry
