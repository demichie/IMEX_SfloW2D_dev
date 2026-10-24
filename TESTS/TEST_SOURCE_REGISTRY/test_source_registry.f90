PROGRAM test_source_registry
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  USE source_catalog_2d
  USE source_registry_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_NEWTON, IM_GROUP_METHOD_ANALYTIC
  USE source_configuration_2d, ONLY : initialize_native_source_configuration, initialize_configured_sources
  USE source_input_2d
  USE parameters_2d, ONLY : wp, n_eqns, idx_pore, &
       rheology_flag, rheology_model, pore_pressure_flag, thermal_cooling_flag, &
       curvature_term_flag, bottom_radial_source_flag, bottom_fissural_source_flag, &
       radial_source_flag, lateral_source_flag, liquid_flag, gas_flag, &
       gas_loss_flag, liquid_vaporization_flag, imex_velocity_group, &
       imex_energy_group, imex_pore_pressure_group, imex_group_method
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
  CALL require(SOURCE_TERM_COUNT==24,'physical catalog contains exactly 24 sources')
  CALL require(source_catalog_names_unique(),'unique physical source names')
  CALL require(canonical_rheology_catalog_unique(),'unique complete rheology names')
  DO i=1,SOURCE_TERM_COUNT
    term=source_term_descriptor_by_name(SOURCE_CATALOG(i)%name)
    CALL require(term%name==SOURCE_CATALOG(i)%name,'physical catalog lookup')
    CALL require(LEN_TRIM(term%target)>0,'physical target declared')
  END DO
  term=source_term_descriptor_by_name('SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS')
  CALL require(term%name=='','legacy mathematical addend removed')
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
      CALL registry%initialize_legacy(partition,ids,methods,MERGE(7,0,enabled(3)),active,error)
      CALL require(error=='' .AND. registry%is_initialized(), &
           'legacy combination normalizes: '//TRIM(error))
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
  CALL verify_native_parser()
  CALL verify_native_schedule()
  WRITE(*,'(A)') 'D17 physical catalog, native configuration and 864 legacy mappings verified'
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
    imex_velocity_group=3;imex_energy_group=1;imex_pore_pressure_group=2
    imex_group_method=[CHARACTER(LEN=16)::'ANALYTIC','ANALYTIC','NEWTON']
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

  SUBROUTINE verify_native_parser()
    INTEGER :: unit,i
    CHARACTER(LEN=160), PARAMETER :: lines(11)=[CHARACTER(LEN=160):: &
         '&source_groups', &
         ' source_group_count=1,', &
         ' source_group_id(1)=5,', &
         ' source_group_order(1)=1,', &
         " source_group_cadence(1)='STAGE',", &
         " source_group_method(1)='NEWTON',", &
         ' source_group_term_count(1)=1,', &
         " source_group_terms(1,1)='RHEOLOGY_LAVA',", &
         ' source_group_target_count(1)=1,', &
         " source_group_targets(1,1)='VELOCITY',", &
         '/']
    OPEN(NEWUNIT=unit,STATUS='scratch',ACTION='readwrite',FORM='formatted')
    DO i=1,SIZE(lines)
      WRITE(unit,'(A)') TRIM(lines(i))
    END DO
    CALL read_native_source_groups(unit,.FALSE.,error)
    CALL require(error=='' .AND. native_source_configuration_present, &
         'native SOURCE_GROUPS namelist parsed')
    CALL require(source_group_terms(1,1)=='RHEOLOGY_LAVA', &
         'native parser preserves user source spelling')
    REWIND(unit)
    CALL read_native_source_groups(unit,.TRUE.,error)
    CALL require(INDEX(error,'cannot be combined')>0, &
         'native and explicit legacy IMEX input conflict')
    CLOSE(unit)
    CALL clear_source_input()
    OPEN(NEWUNIT=unit,STATUS='scratch',ACTION='readwrite',FORM='formatted')
    WRITE(unit,'(A)') '&source_groups'
    WRITE(unit,'(A)') ' source_group_count=1'
    WRITE(unit,'(A)') '/'
    WRITE(unit,'(A)') '&source_groups'
    WRITE(unit,'(A)') ' source_group_count=1'
    WRITE(unit,'(A)') '/'
    CALL read_native_source_groups(unit,.FALSE.,error)
    CALL require(INDEX(error,'may appear only once')>0, &
         'duplicate native namelist rejected')
    CLOSE(unit)
    CALL clear_source_input()
  END SUBROUTINE verify_native_parser

  SUBROUTINE set_native_group(slot,id,cadence,method,order,name,target)
    INTEGER, INTENT(IN) :: slot,id,order
    CHARACTER(*), INTENT(IN) :: cadence,method,name,target
    source_group_id(slot)=id
    source_group_order(slot)=order
    source_group_cadence(slot)=cadence
    source_group_method(slot)=method
    source_group_term_count(slot)=1
    source_group_terms(1,slot)=name
    source_group_target_count(slot)=0
    IF (LEN_TRIM(target)>0) THEN
      source_group_target_count(slot)=1
      source_group_targets(1,slot)=target
    END IF
  END SUBROUTINE set_native_group

  SUBROUTINE verify_native_schedule()
    TYPE(source_group_type), ALLOCATABLE :: g(:),s(:)
    CALL configure_lava_case()
    CALL clear_source_input()
    native_source_configuration_present=.TRUE.
    source_group_count=5
    CALL set_native_group(1,10,'STAGE','ANALYTIC',1,'THERMAL_COOLING','TEMPERATURE')
    CALL set_native_group(2,11,'STAGE','ANALYTIC',2,'PORE_PRESSURE_DISSIPATION','PORE_PRESSURE')
    CALL set_native_group(3,12,'STAGE','NEWTON',3,'RHEOLOGY_LAVA','VELOCITY')
    CALL set_native_group(4,13,'STAGE','EXPLICIT',1,'BOUNDARY_INJECTION','')
    CALL set_native_group(5,14,'STEP','EXPLICIT',1,'MASS_EXCHANGE_DEPOSITION','')
    CALL initialize_native_source_configuration(registry,partition,error)
    CALL require(error=='' .AND. registry%is_initialized(), &
         'native names case-insensitive and complete: '//TRIM(error))
    g=registry%stage_groups();s=registry%step_groups()
    CALL require(SIZE(g)==4 .AND. SIZE(s)==1,'native schedule retains five groups')
    CALL require(group_contains_name(g(1),'thermal_cooling') .AND. &
         group_contains_name(g(3),'rheology_lava'),'native canonical physical names')
    CALL require(group_contains_name(s(1),'mass_exchange_deposition'), &
         'native post-RK physical group')
    source_group_terms(1,3)='SOURCE_TERM_RHEOLOGY_LAVA_VISCOUS'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
         'native mathematical addend rejected')
    source_group_terms(1,3)='RHEOLOGY_LAVA'
    source_group_targets(1,3)='UNKNOWN'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(INDEX(error,'unknown or duplicate group target')>0, &
         'native unsupported target rejected')
    source_group_targets(1,3)='VELOCITY'
    source_group_terms(1,5)='BOUNDARY_INJECTION'
    CALL initialize_native_source_configuration(invalid,partition,error)
    CALL require(error/='' .AND. .NOT.invalid%is_initialized(), &
         'native missing and duplicated physical source rejected')
    CALL clear_source_input()
  END SUBROUTINE verify_native_schedule

  SUBROUTINE require(valid,message)
    LOGICAL, INTENT(IN) :: valid
    CHARACTER(*), INTENT(IN) :: message
    IF (valid) RETURN
    WRITE(error_unit,'(A)') 'FAIL: '//message
    ERROR STOP 1
  END SUBROUTINE require
END PROGRAM test_source_registry
