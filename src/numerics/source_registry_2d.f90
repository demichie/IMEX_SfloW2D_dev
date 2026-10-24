!********************************************************************************
!> \brief Validated, simulation-owned D17 source groups and legacy normalization.
!> \details Builds metadata once, outside cell/stage loops. Does not evaluate
!>          physics or change equation_partition_type, tableaux or solver calls.
!>          Private storage is exposed only through copies; failed initialization
!>          releases all storage. Legacy joint Newton blocks remain coupled.
!********************************************************************************
MODULE source_registry_2d
  USE source_catalog_2d
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_group_type, source_registry_type
  PUBLIC :: group_contains, group_contains_name, group_source_names

  !> \brief Homogeneous source group; identity is cadence plus id, not order.
  TYPE :: source_group_type
    INTEGER :: id=0 !< Public or normalized internal group identity.
    INTEGER :: legacy_id=0 !< Original IMEX input group, zero for non-IM work.
    INTEGER :: cadence=0 !< Stage or post-RK step application.
    INTEGER :: method=0 !< Explicit, Newton or analytic.
    INTEGER :: phase=0 !< Implicit sweep, requested EX evaluation, or post-RK.
    INTEGER :: execution_order=0 !< Explicit order within its phase.
    CHARACTER(LEN=64), ALLOCATABLE :: term_names(:) !< Legacy component names for physical evaluation.
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
    PROCEDURE :: initialize_legacy => initialize_legacy_registry
    PROCEDURE :: finalize => finalize_registry
    PROCEDURE :: is_initialized => registry_is_initialized
    PROCEDURE :: stage_groups => copy_stage_groups
    PROCEDURE :: step_groups => copy_step_groups
    PROCEDURE :: validate => validate_registry
  END TYPE source_registry_type
CONTAINS
  !> \brief Test physical membership without accessing unallocated storage.
  !> \param[in] group Group descriptor.
  !> \param[in] term_id Physical catalog ID.
  !> \return True exactly when the group contains the term.
  PURE LOGICAL FUNCTION group_contains(group,term_id)
    TYPE(source_group_type), INTENT(IN) :: group
    INTEGER, INTENT(IN) :: term_id
    TYPE(source_term_type) :: term
    group_contains=.FALSE.
    term=source_term_descriptor(term_id)
    IF (LEN_TRIM(term%name)==0) RETURN
    ! IDs are compatibility inputs; membership is keyed by canonical name.
    group_contains=group_contains_name(group,term%name)
  END FUNCTION group_contains

  !> Query group membership by source name. The legacy integer-backed group
  !> representation is converted to names at this boundary until migration.
  PURE LOGICAL FUNCTION group_contains_name(group,name) RESULT(found)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(*), INTENT(IN) :: name
    found=.FALSE.
    ! Complete rheological laws and other schedulable sources are addressed
    ! by their canonical identities when the group has been normalized.
    IF (ALLOCATED(group%source_names)) THEN
      found=ANY(group%source_names==TRIM(name))
      IF (found) RETURN
    END IF
    ! Transitional support for callers identifying mathematical addends.
    IF (ALLOCATED(group%term_names)) THEN
      found=ANY(group%term_names==TRIM(name))
    END IF
  END FUNCTION group_contains_name

  !> Return one name per schedulable physical source. Rheological
  !> component identities collapse to a single complete-law name; their
  !> term_names remains an internal compatibility view for physical assemblers.
  FUNCTION group_source_names(group) RESULT(names)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(LEN=64), ALLOCATABLE :: names(:)
    CHARACTER(LEN=64) :: name
    CHARACTER(LEN=64), ALLOCATABLE :: buffer(:)
    INTEGER :: i, n
    ALLOCATE(names(0))
    IF (.NOT.ALLOCATED(group%term_names)) RETURN
    DO i=1,SIZE(group%term_names)
      name=canonical_source_name(group%term_names(i))
      IF (LEN_TRIM(name)==0) CYCLE
      IF (ANY(names==name)) CYCLE
      n=SIZE(names)
      ALLOCATE(buffer(n+1))
      buffer(1:n)=names
      buffer(n+1)=name
      CALL MOVE_ALLOC(buffer,names)
    END DO
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
    CALL populate_registry_term_names(this)
    CALL this%validate(error)
    IF (LEN_TRIM(error)>0) THEN
      CALL this%finalize()
      RETURN
    END IF
    this%initialized=.TRUE.
  END SUBROUTINE initialize_registry

  !> \brief Normalize legacy groups while retaining the original partition.
  !> \param[in,out] this Destination registry, cleared on failure.
  !> \param[in] partition Already validated, compact legacy equation partition.
  !> \param[in] legacy_target_groups Original velocity/temperature/pressure IDs.
  !> \param[in] legacy_methods Three resolved legacy group methods.
  !> \param[in] pore_equation Pore equation index, zero when inactive.
  !> \param[in] active_names Names of potentially active physical addends (cell guards still apply).
  !> \param[out] error Empty on success, otherwise a configuration diagnostic.
  !> \note Mixed analytic groups expand in velocity/temperature/pressure order;
  !>       split IDs start at five. Newton groups retain their full equation block.
  SUBROUTINE initialize_legacy_registry(this,partition,legacy_target_groups,legacy_methods, &
       pore_equation,active_names,error)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    TYPE(equation_partition_type), INTENT(IN) :: partition
    INTEGER, INTENT(IN) :: legacy_target_groups(3),legacy_methods(3)
    INTEGER, INTENT(IN) :: pore_equation
    CHARACTER(*), INTENT(IN) :: active_names(:)
    CHARACTER(*), INTENT(OUT) :: error
    TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:)
    TYPE(source_group_type) :: group
    TYPE(source_term_type) :: term
    INTEGER, ALLOCATABLE :: targets(:),equations(:),target_equations(:)
    CHARACTER(LEN=64), ALLOCATABLE :: names(:)
    INTEGER :: legacy,target,split_id,i,method,order
    LOGICAL :: enabled(3)
    CALL this%finalize()
    error=''
    IF (.NOT.partition%is_initialized()) THEN
      error='D17: invalid legacy equation partition'; RETURN
    END IF
    IF (partition%n_equations<4) THEN
      error='D17: source layout requires at least four equations'; RETURN
    END IF
    IF (ANY(legacy_target_groups<1) .OR. ANY(legacy_target_groups>3)) THEN
      error='D17: legacy group IDs must be from one to three'; RETURN
    END IF
    IF (ANY(legacy_methods/=IM_GROUP_METHOD_NEWTON .AND. &
         legacy_methods/=IM_GROUP_METHOD_ANALYTIC)) THEN
      error='D17: invalid legacy method'; RETURN
    END IF
    IF (partition%implicit(2) .NEQV. partition%implicit(3)) THEN
      error='D17: velocity must remain a two-equation block'; RETURN
    END IF
    IF (pore_equation/=0) THEN
      IF (pore_equation<5 .OR. pore_equation>partition%n_equations) THEN
        error='D17: invalid pore equation index'; RETURN
      END IF
      IF (.NOT.partition%implicit(pore_equation)) THEN
        error='D17: pore target must be implicit in the legacy partition'; RETURN
      END IF
    END IF
    enabled=[partition%implicit(2),partition%implicit(4),pore_equation/=0]
    ! Normalize legacy catalog IDs before checking group membership/targets.
    DO i=1,SIZE(active_names)
      term=source_term_descriptor_by_name(active_names(i))
      IF (.NOT.term%schedulable) THEN
        error='D17: active ID is unknown or dependent'; RETURN
      END IF
      IF (term%target>0) THEN
        IF (.NOT.enabled(term%target)) THEN
          error='D17: active term has no enabled legacy target'; RETURN
        END IF
      END IF
    END DO
    ! Verify the original IDs map to exactly the partition's compact groups.
    order=0
    DO legacy=1,3
      IF (.NOT.ANY(enabled .AND. legacy_target_groups==legacy)) CYCLE
      order=order+1
      IF (order>partition%n_IM_groups) THEN
        error='D17: legacy input does not match equation partition'; RETURN
      END IF
      IF (partition%IM_group_method(order)/=legacy_methods(legacy)) THEN
        error='D17: legacy method does not match equation partition'; RETURN
      END IF
      DO target=1,3
        IF (.NOT.enabled(target) .OR. legacy_target_groups(target)/=legacy) CYCLE
        target_equations=target_equation_ids(target,pore_equation)
        IF (ANY(partition%IM_group(target_equations)/=order)) THEN
          error='D17: legacy target does not match equation partition'; RETURN
        END IF
      END DO
    END DO
    IF (order/=partition%n_IM_groups .OR. &
         partition%n_implicit/=2*MERGE(1,0,enabled(1))+COUNT(enabled(2:3))) THEN
      error='D17: unsupported or missing legacy implicit target'; RETURN
    END IF
    ALLOCATE(stages(0),steps(0))
    split_id=4
    DO legacy=1,3
      targets=PACK([1,2,3],enabled .AND. legacy_target_groups==legacy)
      IF (SIZE(targets)==0) CYCLE
      method=SOURCE_METHOD_NEWTON
      IF (legacy_methods(legacy)==IM_GROUP_METHOD_ANALYTIC) method=SOURCE_METHOD_ANALYTIC
      IF (method==SOURCE_METHOD_NEWTON) THEN
        names=names_for_targets(active_names,targets)
        equations=PACK([(i,i=1,partition%n_equations)],partition%IM_group== &
             partition%IM_group(target_equation_ids_first(targets(1),pore_equation)))
        group=make_group(legacy,legacy,SOURCE_CADENCE_STAGE,method, &
             SOURCE_PHASE_RK_IMPLICIT,SIZE(stages)+1,names,equations)
        stages=[stages,group]
      ELSE
        DO i=1,SIZE(targets)
          target=targets(i)
          names=names_for_targets(active_names,[target])
          equations=target_equation_ids(target,pore_equation)
          order=legacy
          IF (SIZE(targets)>1) THEN
            split_id=split_id+1
            order=split_id
          END IF
          group=make_group(order,legacy,SOURCE_CADENCE_STAGE,method, &
               SOURCE_PHASE_RK_IMPLICIT,SIZE(stages)+1,names,equations)
          stages=[stages,group]
        END DO
      END IF
    END DO
    equations=[(i,i=1,partition%n_equations)]
    names=names_for_targets(active_names,[0])
    group=make_group(4,0,SOURCE_CADENCE_STAGE,SOURCE_METHOD_EXPLICIT, &
         SOURCE_PHASE_RK_EXPLICIT,1,names,equations)
    stages=[stages,group]
    names=names_for_targets(active_names,[-1])
    group=make_group(1,0,SOURCE_CADENCE_STEP,SOURCE_METHOD_EXPLICIT, &
         SOURCE_PHASE_POST_RK,1,names,equations)
    steps=[group]
    CALL this%initialize(partition%n_equations,pore_equation,stages,steps,error)
  END SUBROUTINE initialize_legacy_registry


  !> Materialize source identities once at registry initialization, outside
  !> stage and cell loops. Existing ID-backed consumers remain supported.
  SUBROUTINE populate_registry_term_names(this)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    INTEGER :: i
    DO i=1,SIZE(this%stage_entries)
      this%stage_entries(i)%source_names=group_source_names(this%stage_entries(i))
    END DO
    DO i=1,SIZE(this%step_entries)
      this%step_entries(i)%source_names=group_source_names(this%step_entries(i))
    END DO
  END SUBROUTINE populate_registry_term_names

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
    CHARACTER(LEN=64), ALLOCATABLE :: law_names(:),batch_names(:)
    CHARACTER(LEN=64), ALLOCATABLE :: seen_names(:)
    INTEGER :: i,j,id,nstage,family,mass_group
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
    family=0; mass_group=0
    nstage=SIZE(this%stage_entries)
    groups=[this%stage_entries,this%step_entries]
    DO i=1,SIZE(groups)
      group=groups(i)
      family=0 ! Rheology completeness is checked per group, not across groups.
      error='D17: invalid group identity, method, phase or order'
      IF (group%id<1 .OR. group%execution_order<1 .OR. group%legacy_id<0) RETURN
      IF (i<=nstage) THEN
        IF (group%cadence/=SOURCE_CADENCE_STAGE) RETURN
      ELSE
        IF (group%cadence/=SOURCE_CADENCE_STEP) RETURN
      END IF
      SELECT CASE(group%method)
      CASE(SOURCE_METHOD_EXPLICIT)
        IF (group%cadence==SOURCE_CADENCE_STAGE) THEN
          IF (group%phase/=SOURCE_PHASE_RK_EXPLICIT) RETURN
        ELSE
          IF (group%phase/=SOURCE_PHASE_POST_RK) RETURN
        END IF
      CASE(SOURCE_METHOD_NEWTON,SOURCE_METHOD_ANALYTIC)
        IF (group%cadence/=SOURCE_CADENCE_STAGE .OR. group%phase/=SOURCE_PHASE_RK_IMPLICIT) RETURN
      CASE DEFAULT
        RETURN
      END SELECT
      DO j=1,i-1
        IF (groups(j)%cadence/=group%cadence) CYCLE
        IF (groups(j)%id==group%id) RETURN
        IF (groups(j)%phase>group%phase) RETURN
        IF (groups(j)%phase==group%phase .AND. &
             groups(j)%execution_order>=group%execution_order) RETURN
      END DO
      error='D17: unallocated group descriptor arrays'
      IF (.NOT.ALLOCATED(group%term_names)) RETURN
      IF (.NOT.ALLOCATED(group%equation_ids)) RETURN
      error='D17: invalid or repeated group equation'
      IF (.NOT.valid_equations(group%equation_ids,this%n_equations)) RETURN
      error='D17: repeated or empty canonical source name'
      DO j=1,SIZE(group%term_names)
        IF (LEN_TRIM(group%term_names(j))==0) RETURN
        IF (COUNT(group%term_names==group%term_names(j))/=1) RETURN
      END DO
      IF (ALLOCATED(group%source_names)) THEN
        error='D17: canonical group source identity does not match the complete laws'
        IF (SIZE(group%source_names)/=SIZE(group_source_names(group))) RETURN
        IF (ANY(group%source_names/=group_source_names(group))) RETURN
      END IF
      error='D17: physical term duplicated, unknown or dependent'
      DO j=1,SIZE(group%term_names)
        term=source_term_descriptor_by_name(group%term_names(j))
        ! Unknown names yield an empty descriptor. No numeric identity
        ! is needed to validate membership or scheduling capabilities.
        IF (LEN_TRIM(term%name)==0) RETURN
        IF (.NOT.term%schedulable) RETURN
        IF (ANY(seen_names==term%name)) RETURN
        seen_names=[CHARACTER(LEN=64) :: seen_names,term%name]
      error='D17: term does not support this method/cadence or target'
        IF (.NOT.supports_group(term,group)) RETURN
        IF (term%target/=0) THEN
          expected=target_equation_ids(term%target,this%pore_equation)
          IF (ANY(expected<1)) RETURN
          IF (.NOT.ALL(in_list(expected,group%equation_ids))) RETURN
        END IF
        IF (term%model_id/=0) THEN
          IF (family/=0 .AND. family/=term%model_id) RETURN
          family=term%model_id
        END IF
        IF (term%role==SOURCE_ROLE_MASS) THEN
          IF (mass_group/=0 .AND. mass_group/=i) RETURN
          mass_group=i
        END IF
        error='D17: physical term duplicated, unknown or dependent'
      END DO
      ! A rheological law is one physical source even while the catalog
      ! temporarily contains its legacy formula addends. Never schedule
      ! a proper subset in a Newton or analytic group.
      IF (family/=0 .AND. group%method/=SOURCE_METHOD_EXPLICIT) THEN
        law_names=rheology_component_names(family)
        error='D17: rheology group requires the complete law bundle'
        IF (.NOT.ALL(name_membership(law_names,group%term_names))) RETURN
      END IF
      ! No independently scheduled sub-solves: validate group equation shape.
      error='D17: incomplete or incompatible group equation target'
      IF (group%method==SOURCE_METHOD_NEWTON) THEN
        IF (SIZE(group%equation_ids)==0) RETURN
      ELSE IF (group%method==SOURCE_METHOD_ANALYTIC) THEN
        SELECT CASE(SIZE(group%equation_ids))
        CASE(2)
          IF (.NOT.same_set(group%equation_ids,[2,3])) RETURN
          law_names=rheology_component_names(family)
          IF (SIZE(law_names)/=SIZE(group%term_names) .OR. &
               .NOT.ALL(name_membership(law_names,group%term_names))) THEN
            error='D17: analytic velocity requires the complete law bundle'; RETURN
          END IF
        CASE(1)
          id=group%equation_ids(1)
          IF (id/=4 .AND. id/=this%pore_equation) RETURN
          IF (id==this%pore_equation .AND. this%pore_equation==0) RETURN
        CASE DEFAULT
          RETURN
        END SELECT
      END IF
    END DO
    ! Inspect the actual nominal memberships, not numeric catalog positions.
    ! Any atomic batch touched by a group must be wholly present in one group.
    error='D17: coupled evaluation batch split or incomplete'
    DO i=1,SIZE(groups)
      DO j=1,SIZE(groups(i)%term_names)
        term=source_term_descriptor_by_name(groups(i)%term_names(j))
        IF (term%bundle_id==0) CYCLE
        batch_names=source_bundle_names(term%bundle_id)
        IF (SIZE(batch_names)<2) RETURN
        IF (.NOT.same_batch(batch_names,groups)) RETURN
      END DO
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

  !> \brief Obtain a target's first equation for legacy compact-group lookup.
  !> \param[in] target Valid enabled target.
  !> \param[in] pore Valid pore index when the target is pressure.
  !> \return First conservative index.
  INTEGER FUNCTION target_equation_ids_first(target,pore) RESULT(id)
    INTEGER, INTENT(IN) :: target,pore
    INTEGER, ALLOCATABLE :: ids(:)
    ids=target_equation_ids(target,pore)
    id=ids(1)
  END FUNCTION target_equation_ids_first

  !> \brief Select source names for target blocks, EX (0) or mass exchange (-1).
  !> \param[in] active_names Names normalized at the legacy boundary.
  !> \param[in] targets Target IDs, zero for EX or minus one for step mass.
  !> \return Stable-order canonical addend names; invalid IDs reach validation.
  FUNCTION names_for_targets(active_names,targets) RESULT(names)
    CHARACTER(*), INTENT(IN) :: active_names(:)
    INTEGER, INTENT(IN) :: targets(:)
    CHARACTER(LEN=64), ALLOCATABLE :: names(:)
    TYPE(source_term_type) :: term
    INTEGER :: i,target,n
    ALLOCATE(names(SIZE(active_names)))
    n=0
    DO i=1,SIZE(active_names)
      term=source_term_descriptor_by_name(active_names(i)); target=term%target
      IF (term%role==SOURCE_ROLE_MASS) target=-1
      IF (ANY(targets==target)) THEN
        n=n+1
        names(n)=active_names(i)
      END IF
    END DO
    names=names(:n)
  END FUNCTION names_for_targets

  !> \brief Construct one block with fully allocated arrays.
  !> \brief Construct a fully allocated homogeneous group.
  !> \param[in] id Group identity, independent of order.
  !> \param[in] legacy Original input group, or zero for EX/step.
  !> \param[in] cadence Stage or step namespace.
  !> \param[in] method Group method.
  !> \param[in] phase Evaluation phase within its cadence.
  !> \param[in] order Explicit position within that phase.
  !> \param[in] names Canonical component names, converted at the legacy boundary.
  !> \param[in] equations Conservative equation footprint.
  !> \return Group directly owning the equation footprint.
  FUNCTION make_group(id,legacy,cadence,method,phase,order,names,equations) RESULT(group)
    INTEGER, INTENT(IN) :: id,legacy,cadence,method,phase,order,equations(:)
    CHARACTER(*), INTENT(IN) :: names(:)
    TYPE(source_group_type) :: group
    group%id=id; group%legacy_id=legacy; group%cadence=cadence
    group%method=method; group%phase=phase; group%execution_order=order
    ! Allocate the owned arrays explicitly to avoid compiler warnings about
    ! allocatable-component descriptor initialization on intrinsic assignment.
    ALLOCATE(group%term_names(SIZE(names)),group%equation_ids(SIZE(equations)))
    group%term_names=names
    group%equation_ids=equations
  END FUNCTION make_group

  !> Check catalogued mathematical names without converting group storage to IDs.
  PURE FUNCTION name_membership(required,available) RESULT(found)
    CHARACTER(*), INTENT(IN) :: required(:),available(:)
    LOGICAL :: found(SIZE(required))
    INTEGER :: i
    DO i=1,SIZE(required)
      found(i)=ANY(available==required(i))
    END DO
  END FUNCTION name_membership

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

  !> \brief Require all canonical batch members in the same group, or none.
  !> \param[in] names Canonical names of the atomic-batch addends.
  !> \param[in] groups All stage and step groups in execution order.
  !> \return True only for a wholly absent or unsplit batch.
  LOGICAL FUNCTION same_batch(names,groups) RESULT(valid)
    CHARACTER(*), INTENT(IN) :: names(:)
    TYPE(source_group_type), INTENT(IN) :: groups(:)
    INTEGER :: i,j,owner,found
    valid=.FALSE.
    owner=0
    found=0
    DO i=1,SIZE(names)
      DO j=1,SIZE(groups)
        IF (.NOT.ALLOCATED(groups(j)%term_names)) CYCLE
        IF (.NOT.ANY(groups(j)%term_names==names(i))) CYCLE
        found=found+1
        IF (owner==0) owner=j
        IF (owner/=j) RETURN
      END DO
    END DO
    valid=(found==0 .OR. found==SIZE(names))
  END FUNCTION same_batch

  !> \brief Check implemented temporal capabilities for one physical term.
  !> \param[in] term Valid schedulable catalog descriptor.
  !> \param[in] group Candidate method/cadence.
  !> \return True when the term's existing interface supports this group.
  LOGICAL FUNCTION supports_group(term,group) RESULT(supported)
    TYPE(source_term_type), INTENT(IN) :: term
    TYPE(source_group_type), INTENT(IN) :: group
    supported=.FALSE.
    IF (.NOT.term%schedulable) RETURN
    SELECT CASE(group%method)
    CASE(SOURCE_METHOD_EXPLICIT)
      supported=term%has_explicit
      IF (group%cadence==SOURCE_CADENCE_STAGE) THEN
        supported=supported .AND. term%role/=SOURCE_ROLE_MASS
      ELSE IF (group%cadence==SOURCE_CADENCE_STEP) THEN
        supported=supported .AND. term%role==SOURCE_ROLE_MASS
      ELSE
        supported=.FALSE.
      END IF
    CASE(SOURCE_METHOD_NEWTON)
      supported=group%cadence==SOURCE_CADENCE_STAGE .AND. term%has_newton
    CASE(SOURCE_METHOD_ANALYTIC)
      supported=group%cadence==SOURCE_CADENCE_STAGE .AND. term%has_analytic
    END SELECT
  END FUNCTION supports_group
END MODULE source_registry_2d
