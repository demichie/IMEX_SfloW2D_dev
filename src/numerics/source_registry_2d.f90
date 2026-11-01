!********************************************************************************
!> \brief Validated source groups derived from physical source configuration.
!> \details Builds metadata once, outside cell/stage loops. Does not evaluate
!>          physics or change equation_partition_type, tableaux or solver calls.
!>          Private storage is exposed only through copies; failed initialization
!>          releases all storage. Joint Newton blocks remain coupled.
!********************************************************************************
MODULE source_registry_2d
  USE source_catalog_2d
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_group_type, source_registry_type
  PUBLIC :: group_contains_name, group_source_names

  !> \brief Runtime method subgroup; id=0 is valid.
  TYPE :: source_group_type
    INTEGER :: id=0 !< Public or normalized internal group identity.
    INTEGER :: method=0 !< Explicit, Newton or analytic.
    CHARACTER(LEN=64), ALLOCATABLE :: source_names(:) !< Unique schedulable source identities.
    INTEGER, ALLOCATABLE :: equation_ids(:) !< Implicit equation block or explicit footprint.
  END TYPE source_group_type

  !> \brief Validated metadata owned and finalized by one simulation context.
  TYPE :: source_registry_type
    PRIVATE
    INTEGER :: n_equations=0 !< Conservative layout size.
    INTEGER :: pore_equation=0 !< Pore index, zero when absent.
    LOGICAL :: initialized=.FALSE. !< Cached successful-validation status.
    LOGICAL :: implicit_sources_flag=.FALSE. !< Any validated NEWTON source group.
    TYPE(source_group_type), ALLOCATABLE :: lie_entries(:)
  CONTAINS
    PROCEDURE :: groups => copy_lie_groups
    PROCEDURE :: initialize => initialize_lie_registry
    PROCEDURE :: finalize => finalize_registry
    PROCEDURE :: is_initialized => registry_is_initialized
    PROCEDURE :: has_implicit_sources => registry_has_implicit_sources
  END TYPE source_registry_type
CONTAINS
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

  !> \brief Release every owned descriptor and reset initialization status.
  !> \param[in,out] this Registry to release; repeated calls are safe.
  SUBROUTINE finalize_registry(this)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    IF (ALLOCATED(this%lie_entries)) DEALLOCATE(this%lie_entries)
    this%n_equations=0; this%pore_equation=0; this%initialized=.FALSE.
    this%implicit_sources_flag=.FALSE.
  END SUBROUTINE finalize_registry

  !> \brief Query cached initialization status without walking term arrays.
  !> \param[in] this Registry.
  !> \return True only after successful transactional validation.
  !> Cached query, set once following successful registry validation.
  PURE LOGICAL FUNCTION registry_has_implicit_sources(this)
    CLASS(source_registry_type), INTENT(IN) :: this
    registry_has_implicit_sources = this%initialized .AND. this%implicit_sources_flag
  END FUNCTION registry_has_implicit_sources

  PURE LOGICAL FUNCTION registry_is_initialized(this)
    CLASS(source_registry_type), INTENT(IN) :: this
    registry_is_initialized=this%initialized
  END FUNCTION registry_is_initialized

  !> Initialize one ordered group list for the future Lie integrator.
  SUBROUTINE initialize_lie_registry(this,n_equations,pore_equation,groups,error)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: n_equations,pore_equation
    TYPE(source_group_type), INTENT(IN) :: groups(:)
    CHARACTER(*), INTENT(OUT) :: error
    TYPE(source_term_type) :: term
    CHARACTER(LEN=64), ALLOCATABLE :: seen(:)
    INTEGER :: i,j,mass_group
    CALL this%finalize()
    mass_group=-1
    ALLOCATE(seen(0))
    error='D19: invalid Lie group configuration'
    IF (n_equations<4 .OR. SIZE(groups)==0) RETURN
    IF (pore_equation/=0) THEN
      IF (pore_equation<5 .OR. pore_equation>n_equations) RETURN
    END IF
    IF (groups(1)%id/=0 .OR. groups(1)%method/=SOURCE_METHOD_EXPLICIT) RETURN
    DO i=1,SIZE(groups)
      IF (i>1) THEN
        IF (groups(i)%id<=groups(i-1)%id) RETURN
      END IF
      IF (.NOT.ALLOCATED(groups(i)%source_names)) RETURN
      IF (i>1 .AND. SIZE(groups(i)%source_names)==0) RETURN
      IF (.NOT.ALLOCATED(groups(i)%equation_ids)) RETURN
      IF (.NOT.valid_equations(groups(i)%equation_ids,n_equations)) RETURN
      SELECT CASE(groups(i)%method)
      CASE(SOURCE_METHOD_EXPLICIT,SOURCE_METHOD_NEWTON,SOURCE_METHOD_ANALYTIC)
        CONTINUE
      CASE DEFAULT
        RETURN
      END SELECT
      IF (groups(i)%method==SOURCE_METHOD_NEWTON .AND. SIZE(groups(i)%equation_ids)==0) RETURN
      IF (groups(i)%method==SOURCE_METHOD_ANALYTIC .AND. SIZE(groups(i)%source_names)/=1) RETURN
      DO j=1,SIZE(groups(i)%source_names)
        term=source_term_descriptor_by_name(groups(i)%source_names(j))
        IF (LEN_TRIM(term%name)==0) RETURN
        IF (INDEX(TRIM(term%name),'mass_exchange_')==1) THEN
          IF (mass_group>=0 .AND. mass_group/=groups(i)%id) RETURN
          mass_group=groups(i)%id
        END IF
        IF (ANY(seen==term%name)) RETURN
        seen=[CHARACTER(LEN=64) :: seen,term%name]
        SELECT CASE(groups(i)%method)
        CASE(SOURCE_METHOD_EXPLICIT)
          IF (.NOT.term%has_explicit) RETURN
        CASE(SOURCE_METHOD_NEWTON)
          IF (.NOT.term%has_newton) RETURN
        CASE(SOURCE_METHOD_ANALYTIC)
          IF (.NOT.term%has_analytic) RETURN
        END SELECT
      END DO
    END DO
    this%n_equations=n_equations
    this%pore_equation=pore_equation
    this%lie_entries=groups
    this%implicit_sources_flag=ANY(groups%method==SOURCE_METHOD_NEWTON)
    this%initialized=.TRUE.
    error=''
  END SUBROUTINE initialize_lie_registry

  !> \brief Return a detached copy of the single Lie group registry.
  FUNCTION copy_lie_groups(this) RESULT(groups)
    CLASS(source_registry_type), INTENT(IN) :: this
    TYPE(source_group_type), ALLOCATABLE :: groups(:)
    ALLOCATE(groups(0))
    IF (ALLOCATED(this%lie_entries)) groups=this%lie_entries
  END FUNCTION copy_lie_groups

  !> \brief Map one target to conservative equation indices.
  !> \param[in] target Velocity, energy or pore target ID.
  !> \param[in] pore Pore equation index, zero when absent.
  !> \return Ordered equations; an invalid target returns an empty list.
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

END MODULE source_registry_2d
