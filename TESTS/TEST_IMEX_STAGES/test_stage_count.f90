!> \brief Exercise stage-count validation through real input and workspace entry points.
!> \details Arguments select input parsing or direct workspace initialization and N_RK.
!>          Invalid workspace requests must fail before allocation: dimensions are unset.
PROGRAM test_stage_count
  USE parameters_2d, ONLY : n_RK, n_vars, n_eqns, kt_trace_reconstruction
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE time_integration_2d, ONLY : time_integration_workspace_type
  USE inpout_2d, ONLY : init_param, read_param
  USE source_input_2d, ONLY : native_source_configuration_present, source_group_count
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE source_configuration_2d, ONLY : initialize_native_source_configuration
  USE equation_metadata_2d, ONLY : equation_partition_type
  USE source_catalog_2d, ONLY : SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS, &
       SOURCE_CADENCE_STEP, SOURCE_METHOD_NEWTON, SOURCE_METHOD_EXPLICIT
  IMPLICIT NONE
  TYPE(time_integration_workspace_type) :: integration
  TYPE(equation_partition_type) :: partition
  TYPE(source_registry_type) :: registry
  TYPE(source_group_type), ALLOCATABLE :: stage_groups(:),step_groups(:)
  CHARACTER(LEN=32) :: mode, argument
  CHARACTER(LEN=256) :: error
  INTEGER :: expected

  CALL get_command_argument(1,mode)
  CALL get_command_argument(2,argument)
  READ(argument,*) expected
  SELECT CASE (TRIM(mode))
  CASE ('input')
     CALL init_param
     CALL read_param
     IF (n_RK /= expected) ERROR STOP 'input stage count changed unexpectedly'
  CASE ('native-source')
     CALL init_param
     CALL read_param
     IF (.NOT.native_source_configuration_present) ERROR STOP 'native source namelist was not parsed'
     IF (source_group_count/=2) ERROR STOP 'native source group count changed unexpectedly'
     CALL initialize_native_source_configuration(registry,partition,error)
     IF (LEN_TRIM(error)>0) ERROR STOP 'native source schedule rejected: '//TRIM(error)
     stage_groups=registry%stage_groups()
     step_groups=registry%step_groups()
     IF (SIZE(stage_groups)/=1 .OR. SIZE(step_groups)/=1) &
          ERROR STOP 'native source groups were not built'
     IF (stage_groups(1)%method/=SOURCE_METHOD_NEWTON .OR. &
          SIZE(stage_groups(1)%equation_ids)/=2 .OR. &
          ANY(stage_groups(1)%equation_ids/=[2,3])) &
          ERROR STOP 'native velocity target was not retained'
     IF (step_groups(1)%cadence/=SOURCE_CADENCE_STEP .OR. &
          step_groups(1)%method/=SOURCE_METHOD_EXPLICIT .OR. &
          SIZE(step_groups(1)%term_ids)/=1 .OR. &
          step_groups(1)%term_ids(1)/=SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS) &
          ERROR STOP 'native active mass-exchange term was not retained'
  CASE ('workspace')
     kt_trace_reconstruction = 'froude_blend'
     n_RK = expected
     ! Valid requests need sizes; invalid requests must be rejected before using them.
     IF (expected >= 2 .AND. expected <= 4) THEN
        n_vars=5; n_eqns=5; comp_cells_x=2; comp_cells_y=2
     END IF
     CALL integration%initialize
     CALL integration%finalize
  CASE DEFAULT
     ERROR STOP 'unknown stage-count test mode'
  END SELECT
  IF (kt_trace_reconstruction /= 'froude_blend') ERROR STOP 'unexpected KT trace reconstruction'
  WRITE(*,*) 'PASS: stage count accepted:', n_RK
END PROGRAM test_stage_count
