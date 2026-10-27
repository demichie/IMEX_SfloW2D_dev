!********************************************************************************
!> \brief Validated source groups derived from physical source configuration.
!> \details Builds metadata once, outside cell/stage loops. Does not evaluate
!>          physics or change equation_partition_type, tableaux or solver calls.
!>          Private storage is exposed only through copies; failed initialization
!>          releases all storage. Joint Newton blocks remain coupled.
!********************************************************************************
MODULE source_registry_2d
  USE source_catalog_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_group_type, source_registry_type
  PUBLIC :: group_contains_name, group_source_names

  !> \brief Source group identified by (cadence,method,id); id=0 is valid.
  TYPE :: source_group_type
    INTEGER :: id=0 !< Public or normalized internal group identity.
    INTEGER :: cadence=0 !< Stage or post-RK step application.
    INTEGER :: method=0 !< Explicit, Newton or analytic.
    CHARACTER(LEN=64), ALLOCATABLE :: source_names(:) !< Unique schedulable source identities.
    INTEGER, ALLOCATABLE :: equation_ids(:) !< Union of target equations or explicit footprint.
  END TYPE source_group_type

  !> \brief Validated metadata owned and finalized by one simulation context.
  TYPE :: source_registry_type
    PRIVATE
    INTEGER :: n_equations=0 !< Conservative layout size.
    INTEGER :: pore_equation=0 !< Pore index, zero when absent.
    LOGICAL :: initialized=.FALSE. !< Cached successful-validation status.
    TYPE(source_group_type), ALLOCATABLE :: stage_entries(:), step_entries(:)
  CONTAINS
    PROCEDURE :: initialize => initialize_registry
    PROCEDURE :: finalize => finalize_registry
    PROCEDURE :: is_initialized => registry_is_initialized
    PROCEDURE :: stage_groups => copy_stage_groups
    PROCEDURE :: step_groups => copy_step_groups
    PROCEDURE :: validate => validate_registry
  END TYPE source_registry_type
CONTAINS
  PURE INTEGER FUNCTION group_stage_rank(group) RESULT(rank)
    TYPE(source_group_type), INTENT(IN) :: group
    rank=1
    IF (group%method==SOURCE_METHOD_EXPLICIT) rank=2
  END FUNCTION group_stage_rank

  !> Membership in the single, canonical physical-source name array.
  PURE LOGICAL FUNCTION group_contains_name(group,name) RESULT(found)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(*), INTENT(IN) :: name
    found=.FALSE.
    IF (.NOT.ALLOCATED(group%source_names)) RETURN
    found=ANY(group%source_names==TRIM(name))
  END FUNCTION group_contains_name

  !> Return the group-owned names without an additional canonicalization layer.
  FUNCTION group_source_names(group) RESULT(names)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(LEN=64), ALLOCATABLE :: names(:)
    ALLOCATE(names(0))
    IF (ALLOCATED(group%source_names)) names=group%source_names
  END FUNCTION group_source_names

  !> \brief Copy and validate complete group/solve metadata transactionally.
  !> \param[in,out] this Registry; reset even when validation fails.
  !> \param[in] n_equations Conservative equation count.
  !> \param[in] pore_equation Optional pore equation index, or zero when absent.
  !> \param[in] stage_groups Ordered stage descriptors.
  !> \param[in] step_groups Ordered post-RK descriptors.
  !> \param[out] error Empty on success, otherwise a diagnostic before simulation.
  SUBROUTINE initialize_registry(this,n_equations,pore_equation,stage_groups,step_groups,error)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: n_equations,pore_equation
    TYPE(source_group_type), INTENT(IN) :: stage_groups(:),step_groups(:)
    CHARACTER(*), INTENT(OUT) :: error
    CALL this%finalize()
    this%n_equations=n_equations
    this%pore_equation=pore_equation
    this%stage_entries=stage_groups
    this%step_entries=step_groups
    CALL this%validate(error)
    IF (LEN_TRIM(error)>0) THEN
      CALL this%finalize()
      RETURN
    END IF
    this%initialized=.TRUE.
  END SUBROUTINE initialize_registry

  !> \brief Validate native or normalized descriptors before making them visible.
  !> \param[in] this Candidate metadata; validation does not mutate it.
  !> \param[out] error Empty when all membership/capability/block checks pass.
  SUBROUTINE validate_registry(this,error)
    CLASS(source_registry_type), INTENT(IN) :: this
    CHARACTER(*), INTENT(OUT) :: error
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    TYPE(source_group_type) :: group
    TYPE(source_term_type) :: term
    INTEGER, ALLOCATABLE :: expected(:)
    CHARACTER(LEN=64) :: rheology_name
    CHARACTER(LEN=64), ALLOCATABLE :: seen_names(:)
    INTEGER :: i,j,id,nstage,mass_group
    error='D17: invalid or duplicate source catalog name'
    IF (.NOT.source_catalog_names_unique()) RETURN
    error='D17: incomplete source registry'
    IF (.NOT.ALLOCATED(this%stage_entries)) RETURN
    IF (.NOT.ALLOCATED(this%step_entries)) RETURN
    IF (this%n_equations<4) RETURN
    IF (this%pore_equation/=0) THEN
      IF (this%pore_equation<5 .OR. this%pore_equation>this%n_equations) RETURN
    END IF
    ALLOCATE(seen_names(0))
    mass_group=0
    nstage=SIZE(this%stage_entries)
    groups=[this%stage_entries,this%step_entries]
    DO i=1,SIZE(groups)
      group=groups(i)
      rheology_name='' ! Rheology is one schedulable physical identity.
      error='D17: invalid group identity, method, cadence or order'
      IF (group%id<0) RETURN
      IF (i<=nstage) THEN
        IF (group%cadence/=SOURCE_CADENCE_STAGE) RETURN
      ELSE
        IF (group%cadence/=SOURCE_CADENCE_STEP) RETURN
      END IF
      SELECT CASE(group%method)
      CASE(SOURCE_METHOD_EXPLICIT)
        CONTINUE
      CASE(SOURCE_METHOD_NEWTON,SOURCE_METHOD_ANALYTIC)
        IF (group%cadence/=SOURCE_CADENCE_STAGE .AND. &
             group%cadence/=SOURCE_CADENCE_STEP) RETURN
      CASE DEFAULT
        RETURN
      END SELECT
      DO j=1,i-1
        IF (groups(j)%cadence/=group%cadence) CYCLE
        IF (groups(j)%id==group%id .AND. groups(j)%method==group%method) RETURN
        IF (group_stage_rank(groups(j))>group_stage_rank(group)) RETURN
      END DO
      error='D17: unallocated group descriptor arrays'
      IF (.NOT.ALLOCATED(group%source_names)) RETURN
      IF (.NOT.ALLOCATED(group%equation_ids)) RETURN
      error='D17: invalid or repeated group equation'
      IF (.NOT.valid_equations(group%equation_ids,this%n_equations)) RETURN
      error='D17: repeated or empty canonical source name'
      DO j=1,SIZE(group%source_names)
        IF (LEN_TRIM(group%source_names(j))==0) RETURN
        IF (COUNT(group%source_names==group%source_names(j))/=1) RETURN
      END DO
      error='D17: physical term duplicated, unknown or dependent'
      DO j=1,SIZE(group%source_names)
        term=source_term_descriptor_by_name(group%source_names(j))
        ! Unknown names yield an empty descriptor. No numeric identity
        ! is needed to validate membership or scheduling capabilities.
        IF (LEN_TRIM(term%name)==0) RETURN
        ! Mathematical rheology addends are not physical schedule entries
        ! evaluation; native and normalized groups own complete laws.
                IF (ANY(seen_names==term%name)) RETURN
        seen_names=[CHARACTER(LEN=64) :: seen_names,term%name]
      error='D17: term does not support this method/cadence or target'
        IF (.NOT.supports_group(term,group)) RETURN
        IF (target_index(term%target)>0) THEN
          expected=target_equation_ids(target_index(term%target),this%pore_equation)
          IF (ANY(expected<1)) RETURN
          IF (.NOT.ALL(in_list(expected,group%equation_ids))) RETURN
        END IF
        IF (INDEX(TRIM(term%name),'rheology_')==1) THEN
          IF (LEN_TRIM(rheology_name)>0 .AND. rheology_name/=term%name) RETURN
          rheology_name=term%name
        END IF
        IF (INDEX(TRIM(term%name),'mass_exchange_')==1) THEN
          IF (mass_group/=0 .AND. mass_group/=i) RETURN
          mass_group=i
        END IF
        error='D17: physical term duplicated, unknown or dependent'
      END DO
      IF (LEN_TRIM(rheology_name)>0 .AND. group%method/=SOURCE_METHOD_EXPLICIT) THEN
        error='D17: rheology group requires the complete law'
        IF (.NOT.ANY(group%source_names==rheology_name)) RETURN
      END IF
      ! No independently scheduled sub-solves: validate group equation shape.
      error='D17: incomplete or incompatible group equation target'
      IF (group%method==SOURCE_METHOD_NEWTON) THEN
        IF (SIZE(group%equation_ids)==0) RETURN
      ELSE IF (group%method==SOURCE_METHOD_ANALYTIC) THEN
        ! A group may combine independent complete analytic operators.
        ! The momentum pair is indivisible; energy and pore are scalar.
        IF (SIZE(group%equation_ids)==0) RETURN
        IF (ANY(group%equation_ids==2) .NEQV. &
             ANY(group%equation_ids==3)) RETURN
        DO j=1,SIZE(group%equation_ids)
          id=group%equation_ids(j)
          IF (id/=2 .AND. id/=3 .AND. id/=4 .AND. &
               (this%pore_equation==0 .OR. id/=this%pore_equation)) RETURN
        END DO
        IF (ANY(group%equation_ids==2) .AND. LEN_TRIM(rheology_name)==0) RETURN
      END IF
    END DO
    error=''
  END SUBROUTINE validate_registry

  !> \brief Release every owned descriptor and reset initialization status.
  !> \param[in,out] this Registry to release; repeated calls are safe.
  SUBROUTINE finalize_registry(this)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    IF (ALLOCATED(this%stage_entries)) DEALLOCATE(this%stage_entries)
    IF (ALLOCATED(this%step_entries)) DEALLOCATE(this%step_entries)
    this%n_equations=0; this%pore_equation=0; this%initialized=.FALSE.
  END SUBROUTINE finalize_registry

  !> \brief Query cached initialization status without walking term arrays.
  !> \param[in] this Registry.
  !> \return True only after successful transactional validation.
  PURE LOGICAL FUNCTION registry_is_initialized(this)
    CLASS(source_registry_type), INTENT(IN) :: this
    registry_is_initialized=this%initialized
  END FUNCTION registry_is_initialized

  !> \brief Return a detached copy of ordered stage metadata.
  !> \param[in] this Registry.
  !> \return Stage groups, or an allocated empty list after finalization.
  FUNCTION copy_stage_groups(this) RESULT(groups)
    CLASS(source_registry_type), INTENT(IN) :: this
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    ALLOCATE(groups(0))
    IF (ALLOCATED(this%stage_entries)) groups=this%stage_entries
  END FUNCTION copy_stage_groups

  !> \brief Return a detached copy of post-RK metadata.
  !> \param[in] this Registry.
  !> \return Step groups, or an allocated empty list after finalization.
  FUNCTION copy_step_groups(this) RESULT(groups)
    CLASS(source_registry_type), INTENT(IN) :: this
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    ALLOCATE(groups(0))
    IF (ALLOCATED(this%step_entries)) groups=this%step_entries
  END FUNCTION copy_step_groups

  !> \brief Map one target to conservative equation indices.
  !> \param[in] target Velocity, energy or pore target ID.
  !> \param[in] pore Pore equation index, zero when absent.
  !> \return Ordered equations; an invalid target returns an empty list.
  FUNCTION target_equation_ids(target,pore) RESULT(ids)
    INTEGER, INTENT(IN) :: target,pore
    INTEGER, ALLOCATABLE :: ids(:)
    SELECT CASE(target)
    CASE(SOURCE_TARGET_VELOCITY); ids=[2,3]
    CASE(SOURCE_TARGET_TEMPERATURE); ids=[4]
    CASE(SOURCE_TARGET_PORE_PRESSURE); ids=[pore]
    CASE DEFAULT; ids=[INTEGER::]
    END SELECT
  END FUNCTION target_equation_ids

  INTEGER FUNCTION target_index(name) RESULT(target)
    CHARACTER(*), INTENT(IN) :: name
    target=0
    SELECT CASE(TRIM(name))
    CASE('velocity'); target=1
    CASE('temperature'); target=2
    CASE('pore_pressure'); target=3
    CASE('mass_exchange'); target=-1
    END SELECT
  END FUNCTION target_index
  !> \brief Check index bounds and uniqueness before any vector indexing.
  !> \param[in] ids Conservative equation IDs.
  !> \param[in] n Equation count.
  !> \return True for a unique in-range set, including an empty explicit set.
  LOGICAL FUNCTION valid_equations(ids,n) RESULT(valid)
    INTEGER, INTENT(IN) :: ids(:),n
    INTEGER :: i
    valid=.FALSE.
    IF (ANY(ids<1) .OR. ANY(ids>n)) RETURN
    DO i=1,SIZE(ids)
      IF (COUNT(ids==ids(i))/=1) RETURN
    END DO
    valid=.TRUE.
  END FUNCTION valid_equations

  !> \brief Elementwise membership, also for empty descriptor lists.
  !> \param[in] ids Requested IDs.
  !> \param[in] list Candidate set.
  !> \return One flag per requested ID.
  FUNCTION in_list(ids,list) RESULT(found)
    INTEGER, INTENT(IN) :: ids(:),list(:)
    LOGICAL :: found(SIZE(ids))
    INTEGER :: i
    DO i=1,SIZE(ids)
      found(i)=ANY(list==ids(i))
    END DO
  END FUNCTION in_list

  !> \brief Compare unique ID sets without changing their recorded order.
  !> \param[in] left First list.
  !> \param[in] right Second list.
  !> \return True for equal-sized unique lists with identical membership.
  LOGICAL FUNCTION same_set(left,right) RESULT(equal)
    INTEGER, INTENT(IN) :: left(:),right(:)
    INTEGER :: i
    equal=.FALSE.
    IF (SIZE(left)/=SIZE(right)) RETURN
    IF (.NOT.ALL(in_list(left,right))) RETURN
    DO i=1,SIZE(left)
      IF (COUNT(left==left(i))/=1 .OR. COUNT(right==left(i))/=1) RETURN
    END DO
    equal=.TRUE.
  END FUNCTION same_set

  !> \brief Check implemented temporal capabilities for one physical term.
  !> \param[in] term Valid schedulable catalog descriptor.
  !> \param[in] group Candidate method/cadence.
  !> \return True when the term's existing interface supports this group.
  LOGICAL FUNCTION supports_group(term,group) RESULT(supported)
    TYPE(source_term_type), INTENT(IN) :: term
    TYPE(source_group_type), INTENT(IN) :: group
    supported=.FALSE.
    IF (LEN_TRIM(term%name)==0) RETURN
    SELECT CASE(group%method)
    CASE(SOURCE_METHOD_EXPLICIT)
      supported=term%has_explicit
      IF (group%cadence==SOURCE_CADENCE_STAGE) THEN
        supported=supported .AND. INDEX(TRIM(term%name),'mass_exchange_')/=1
      ELSE IF (group%cadence==SOURCE_CADENCE_STEP) THEN
        supported=supported
      ELSE
        supported=.FALSE.
      END IF
    CASE(SOURCE_METHOD_NEWTON)
      supported=(group%cadence==SOURCE_CADENCE_STAGE .OR. &
           group%cadence==SOURCE_CADENCE_STEP) .AND. term%has_newton
    CASE(SOURCE_METHOD_ANALYTIC)
      supported=(group%cadence==SOURCE_CADENCE_STAGE .OR. &
           group%cadence==SOURCE_CADENCE_STEP) .AND. term%has_analytic
    END SELECT
  END FUNCTION supports_group
END MODULE source_registry_2d
