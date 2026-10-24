PROGRAM test_source_registry
  USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
  USE source_catalog_2d
  USE source_registry_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_NEWTON, IM_GROUP_METHOD_ANALYTIC
  USE parameters_2d, ONLY : rheology_id_lava
  IMPLICIT NONE
  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: registry,invalid
  TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:),wrong(:)
  TYPE(source_term_type) :: term
  CHARACTER(LEN=64), ALLOCATABLE :: active(:)
  CHARACTER(LEN=256) :: error
  INTEGER :: i,j,mode,a,b,c,bits,ids(3),methods(3),membership(7),nm
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
  WRITE(*,'(A)') 'D17 physical source catalog and 864 legacy mappings verified'
CONTAINS
  SUBROUTINE require(valid,message)
    LOGICAL, INTENT(IN) :: valid
    CHARACTER(*), INTENT(IN) :: message
    IF (valid) RETURN
    WRITE(error_unit,'(A)') 'FAIL: '//message
    ERROR STOP 1
  END SUBROUTINE require
END PROGRAM test_source_registry
