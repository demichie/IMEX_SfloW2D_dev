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
  PUBLIC :: source_group_type, source_solver_type, source_registry_type
  PUBLIC :: group_contains, group_contains_name

  !> \brief One ordered numerical block, containing all of its physical terms.
  TYPE :: source_solver_type
    INTEGER :: id=0 !< Unique positive registry-local solver ID.
    INTEGER :: kind=0 !< Newton block or a specific analytic target solve.
    INTEGER :: execution_order=0 !< Position within its group.
    INTEGER, ALLOCATABLE :: target_ids(:) !< Scalar/vector target blocks.
    INTEGER, ALLOCATABLE :: equation_ids(:) !< Full conservative equation indices.
    INTEGER, ALLOCATABLE :: term_ids(:) !< Temporary compatibility representation.
    CHARACTER(LEN=64), ALLOCATABLE :: term_names(:) !< Canonical physical source names.
  END TYPE source_solver_type

  !> \brief Homogeneous source group; identity is cadence plus id, not order.
  TYPE :: source_group_type
    INTEGER :: id=0 !< Public or normalized internal group identity.
    INTEGER :: legacy_id=0 !< Original IMEX input group, zero for non-IM work.
    INTEGER :: cadence=0 !< Stage or post-RK step application.
    INTEGER :: method=0 !< Explicit, Newton or analytic.
    INTEGER :: phase=0 !< Implicit sweep, requested EX evaluation, or post-RK.
    INTEGER :: execution_order=0 !< Explicit order within its phase.
    INTEGER, ALLOCATABLE :: term_ids(:) !< Temporary compatibility representation.
    CHARACTER(LEN=64), ALLOCATABLE :: term_names(:) !< Canonical physical source names.
    INTEGER, ALLOCATABLE :: equation_ids(:) !< Union of target equations or explicit footprint.
    INTEGER, ALLOCATABLE :: implicit_solver_ids(:) !< Registry-local Newton solve IDs.
    INTEGER, ALLOCATABLE :: analytic_solver_ids(:) !< Exactly one for an analytic group.
  END TYPE source_group_type

  !> \brief Validated metadata owned and finalized by one simulation context.
  TYPE :: source_registry_type
    PRIVATE
    INTEGER :: n_equations=0 !< Conservative layout size.
    INTEGER :: pore_equation=0 !< Pore index, zero when absent.
    LOGICAL :: initialized=.FALSE. !< Cached successful-validation status.
    TYPE(source_group_type), ALLOCATABLE :: stage_entries(:), step_entries(:)
    TYPE(source_solver_type), ALLOCATABLE :: solver_entries(:)
  CONTAINS
    PROCEDURE :: initialize => initialize_registry
    PROCEDURE :: initialize_legacy => initialize_legacy_registry
    PROCEDURE :: finalize => finalize_registry
    PROCEDURE :: is_initialized => registry_is_initialized
    PROCEDURE :: stage_groups => copy_stage_groups
    PROCEDURE :: step_groups => copy_step_groups
    PROCEDURE :: solvers => copy_solvers
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
    group_contains=.FALSE.
    IF (ALLOCATED(group%term_ids)) group_contains=ANY(group%term_ids==term_id)
  END FUNCTION group_contains

  !> Query group membership by source name. The legacy integer-backed group
  !> representation is converted to names at this boundary until migration.
  PURE LOGICAL FUNCTION group_contains_name(group,name) RESULT(found)
    TYPE(source_group_type), INTENT(IN) :: group
    CHARACTER(*), INTENT(IN) :: name
    found=.FALSE.
    IF (ALLOCATED(group%term_names)) THEN
      found=ANY(group%term_names==TRIM(name))
    ELSE IF (ALLOCATED(group%term_ids)) THEN
      found=ANY(source_names_for_ids(group%term_ids)==TRIM(name))
    END IF
  END FUNCTION group_contains_name

  !> \brief Copy and validate complete group/solve metadata transactionally.
  !> \param[in,out] this Registry; reset even when validation fails.
  !> \param[in] n_equations Conservative equation count.
  !> \param[in] pore_equation Optional pore equation index, or zero when absent.
  !> \param[in] stage_groups Ordered stage descriptors.
  !> \param[in] step_groups Ordered post-RK descriptors.
  !> \param[in] solvers Ordered block descriptors.
  !> \param[out] error Empty on success, otherwise a diagnostic before simulation.
  SUBROUTINE initialize_registry(this,n_equations,pore_equation,stage_groups,step_groups,solvers,error)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: n_equations,pore_equation
    TYPE(source_group_type), INTENT(IN) :: stage_groups(:),step_groups(:)
    TYPE(source_solver_type), INTENT(IN) :: solvers(:)
    CHARACTER(*), INTENT(OUT) :: error
    CALL this%finalize()
    this%n_equations=n_equations
    this%pore_equation=pore_equation
    this%stage_entries=stage_groups
    this%step_entries=step_groups
    this%solver_entries=solvers
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
  !> \param[in] active_terms Potentially active physical IDs (cell guards still apply).
  !> \param[out] error Empty on success, otherwise a configuration diagnostic.
  !> \note Mixed analytic groups expand in velocity/temperature/pressure order;
  !>       split IDs start at five. Newton groups retain their full equation block.
  SUBROUTINE initialize_legacy_registry(this,partition,legacy_target_groups,legacy_methods, &
       pore_equation,active_terms,error)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    TYPE(equation_partition_type), INTENT(IN) :: partition
    INTEGER, INTENT(IN) :: legacy_target_groups(3),legacy_methods(3)
    INTEGER, INTENT(IN) :: pore_equation,active_terms(:)
    CHARACTER(*), INTENT(OUT) :: error
    TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:)
    TYPE(source_solver_type), ALLOCATABLE :: solves(:)
    TYPE(source_group_type) :: group
    TYPE(source_solver_type) :: solver
    TYPE(source_term_type) :: term
    INTEGER, ALLOCATABLE :: targets(:),terms(:),equations(:),target_equations(:)
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
    DO i=1,SIZE(active_terms)
      term=source_term_descriptor(active_terms(i))
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
    ALLOCATE(stages(0),steps(0),solves(0))
    split_id=4
    DO legacy=1,3
      targets=PACK([1,2,3],enabled .AND. legacy_target_groups==legacy)
      IF (SIZE(targets)==0) CYCLE
      method=SOURCE_METHOD_NEWTON
      IF (legacy_methods(legacy)==IM_GROUP_METHOD_ANALYTIC) method=SOURCE_METHOD_ANALYTIC
      IF (method==SOURCE_METHOD_NEWTON) THEN
        terms=terms_for_targets(active_terms,targets)
        equations=PACK([(i,i=1,partition%n_equations)],partition%IM_group== &
             partition%IM_group(target_equation_ids_first(targets(1),pore_equation)))
        solver=make_solver(SIZE(solves)+1,SOURCE_SOLVE_NEWTON_BLOCK,targets,equations,terms)
        solves=[solves,solver]
        group=make_group(legacy,legacy,SOURCE_CADENCE_STAGE,method, &
             SOURCE_PHASE_RK_IMPLICIT,SIZE(stages)+1,terms,equations,[solver%id])
        stages=[stages,group]
      ELSE
        DO i=1,SIZE(targets)
          target=targets(i)
          terms=terms_for_targets(active_terms,[target])
          equations=target_equation_ids(target,pore_equation)
          solver=make_solver(SIZE(solves)+1,SOURCE_SOLVE_ANALYTIC_VELOCITY+target-1, &
               [target],equations,terms)
          solves=[solves,solver]
          order=legacy
          IF (SIZE(targets)>1) THEN
            split_id=split_id+1
            order=split_id
          END IF
          group=make_group(order,legacy,SOURCE_CADENCE_STAGE,method, &
               SOURCE_PHASE_RK_IMPLICIT,SIZE(stages)+1,terms,equations,[solver%id])
          stages=[stages,group]
        END DO
      END IF
    END DO
    equations=[(i,i=1,partition%n_equations)]
    terms=terms_for_targets(active_terms,[0])
    group=make_group(4,0,SOURCE_CADENCE_STAGE,SOURCE_METHOD_EXPLICIT, &
         SOURCE_PHASE_RK_EXPLICIT,1,terms,equations,[INTEGER::])
    stages=[stages,group]
    terms=terms_for_targets(active_terms,[-1])
    group=make_group(1,0,SOURCE_CADENCE_STEP,SOURCE_METHOD_EXPLICIT, &
         SOURCE_PHASE_POST_RK,1,terms,equations,[INTEGER::])
    steps=[group]
    CALL this%initialize(partition%n_equations,pore_equation,stages,steps,solves,error)
  END SUBROUTINE initialize_legacy_registry

  !> Translate legacy catalog positions to their canonical source names.
  PURE FUNCTION source_names_for_ids(ids) RESULT(names)
    INTEGER, INTENT(IN) :: ids(:)
    CHARACTER(LEN=64) :: names(SIZE(ids))
    TYPE(source_term_type) :: term
    INTEGER :: i

    names=''
    DO i=1,SIZE(ids)
      term=source_term_descriptor(ids(i))
      names(i)=term%name
    END DO
  END FUNCTION source_names_for_ids

  !> Materialize source identities once at registry initialization, outside
  !> stage and cell loops. Existing ID-backed consumers remain supported.
  SUBROUTINE populate_registry_term_names(this)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    INTEGER :: i

    DO i=1,SIZE(this%stage_entries)
      IF (ALLOCATED(this%stage_entries(i)%term_ids)) &
           this%stage_entries(i)%term_names=source_names_for_ids(this%stage_entries(i)%term_ids)
    END DO
    DO i=1,SIZE(this%step_entries)
      IF (ALLOCATED(this%step_entries(i)%term_ids)) &
           this%step_entries(i)%term_names=source_names_for_ids(this%step_entries(i)%term_ids)
    END DO
    DO i=1,SIZE(this%solver_entries)
      IF (ALLOCATED(this%solver_entries(i)%term_ids)) &
           this%solver_entries(i)%term_names=source_names_for_ids(this%solver_entries(i)%term_ids)
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
    TYPE(source_solver_type) :: solver
    TYPE(source_term_type) :: term
    INTEGER, ALLOCATABLE :: expected(:),terms(:),bundle_terms(:)
    LOGICAL :: seen(SOURCE_TERM_COUNT)
    INTEGER :: i,j,k,id,si,nstage,family,mass_group,term_group(SOURCE_TERM_COUNT)
    LOGICAL, ALLOCATABLE :: used_solver(:)
    error='D17: invalid or duplicate source catalog name'
    IF (.NOT.source_catalog_names_unique()) RETURN
    error='D17: incomplete source registry'
    IF (.NOT.ALLOCATED(this%stage_entries)) RETURN
    IF (.NOT.ALLOCATED(this%step_entries)) RETURN
    IF (.NOT.ALLOCATED(this%solver_entries)) RETURN
    IF (this%n_equations<4) RETURN
    IF (this%pore_equation/=0) THEN
      IF (this%pore_equation<5 .OR. this%pore_equation>this%n_equations) RETURN
    END IF
    ALLOCATE(used_solver(SIZE(this%solver_entries)))
    used_solver=.FALSE.; seen=.FALSE.; term_group=0; family=0; mass_group=0
    nstage=SIZE(this%stage_entries)
    groups=[this%stage_entries,this%step_entries]
    DO i=1,SIZE(groups)
      group=groups(i)
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
      IF (.NOT.ALLOCATED(group%term_ids)) RETURN
      IF (.NOT.ALLOCATED(group%equation_ids)) RETURN
      IF (.NOT.ALLOCATED(group%implicit_solver_ids)) RETURN
      IF (.NOT.ALLOCATED(group%analytic_solver_ids)) RETURN
      error='D17: invalid or repeated group equation'
      IF (.NOT.valid_equations(group%equation_ids,this%n_equations)) RETURN
      ! Name arrays are materialized on registry initialization. Keep direct
      ! validation of pre-normalized legacy descriptors backward-compatible.
      IF (ALLOCATED(group%term_names)) THEN
        error='D17: group source name and legacy membership disagree'
        IF (SIZE(group%term_names)/=SIZE(group%term_ids)) RETURN
        IF (ANY(group%term_names/=source_names_for_ids(group%term_ids))) RETURN
      END IF
      error='D17: physical term duplicated, unknown or dependent'
      DO j=1,SIZE(group%term_ids)
        id=group%term_ids(j)
        ! Validated registry groups resolve their source through the canonical
        ! name. Raw legacy descriptors are still accepted by direct validate().
        IF (ALLOCATED(group%term_names)) THEN
          term=source_term_descriptor_by_name(group%term_names(j))
        ELSE
          term=source_term_descriptor(id)
        END IF
        IF (term%id/=id .OR. .NOT.term%schedulable) RETURN
        IF (seen(id)) RETURN
        seen(id)=.TRUE.; term_group(id)=i
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
      error='D17: invalid solver count for group method'
      SELECT CASE(group%method)
      CASE(SOURCE_METHOD_EXPLICIT)
        IF (SIZE(group%implicit_solver_ids)+SIZE(group%analytic_solver_ids)/=0) RETURN
      CASE(SOURCE_METHOD_NEWTON)
        IF (SIZE(group%implicit_solver_ids)/=1 .OR. SIZE(group%analytic_solver_ids)/=0) RETURN
      CASE(SOURCE_METHOD_ANALYTIC)
        IF (SIZE(group%analytic_solver_ids)/=1 .OR. SIZE(group%implicit_solver_ids)/=0) RETURN
      END SELECT
      ! Transitional descriptor: exactly one solve per non-explicit group.
      IF (group%method/=SOURCE_METHOD_EXPLICIT) THEN
        error='D17: missing, duplicate or inconsistent solver descriptor'
        IF (group%method==SOURCE_METHOD_NEWTON) THEN
          si=group%implicit_solver_ids(1)
        ELSE
          si=group%analytic_solver_ids(1)
        END IF
        IF (si<1 .OR. si>SIZE(this%solver_entries)) RETURN
        IF (used_solver(si)) RETURN
        used_solver(si)=.TRUE.; solver=this%solver_entries(si)
        IF (solver%id/=si .OR. solver%execution_order/=1) RETURN
        IF (.NOT.ALLOCATED(solver%target_ids)) RETURN
        IF (.NOT.ALLOCATED(solver%equation_ids)) RETURN
        IF (.NOT.ALLOCATED(solver%term_ids)) RETURN
        IF (SIZE(solver%target_ids)<1) RETURN
        IF (.NOT.valid_equations(solver%equation_ids,this%n_equations)) RETURN
        IF (SIZE(solver%equation_ids)==0) RETURN
        expected=[INTEGER::]
        DO k=1,SIZE(solver%target_ids)
          id=solver%target_ids(k)
          IF (id<1 .OR. id>3) RETURN
          IF (COUNT(solver%target_ids==id)/=1) RETURN
          expected=[expected,target_equation_ids(id,this%pore_equation)]
        END DO
        IF (.NOT.same_set(expected,solver%equation_ids)) RETURN
        ! With one solve per implicit group, its equation map is the group map.
        IF (.NOT.same_set(solver%equation_ids,group%equation_ids)) RETURN
        IF (group%method==SOURCE_METHOD_NEWTON) THEN
          IF (solver%kind/=SOURCE_SOLVE_NEWTON_BLOCK) RETURN
        ELSE
          IF (SIZE(solver%target_ids)/=1) RETURN
          IF (solver%kind/=SOURCE_SOLVE_ANALYTIC_VELOCITY+solver%target_ids(1)-1) RETURN
        END IF
        DO k=1,SIZE(solver%term_ids)
          term=source_term_descriptor(solver%term_ids(k))
          IF (.NOT.ANY(solver%target_ids==term%target)) RETURN
        END DO
        IF (group%method==SOURCE_METHOD_ANALYTIC) THEN
          IF (solver%target_ids(1)==SOURCE_TARGET_VELOCITY) THEN
            terms=rheology_source_terms(family)
            IF (.NOT.same_set(terms,solver%term_ids)) THEN
              error='D17: analytic velocity requires the complete law bundle'; RETURN
            END IF
          END IF
        END IF
      END IF
      error='D17: solver coverage does not match group equations/terms'
      IF (group%method/=SOURCE_METHOD_EXPLICIT) THEN
        ! A single solve must contain precisely the group's physical sources.
        IF (.NOT.same_set(solver%term_ids,group%term_ids)) RETURN
      END IF
    END DO
    error='D17: orphan solver descriptor'
    IF (.NOT.ALL(used_solver)) RETURN
    error='D17: coupled evaluation batch split or incomplete'
    DO id=1,SOURCE_TERM_COUNT
      term=source_term_descriptor(id)
      IF (term%bundle_id==0) CYCLE
      bundle_terms=source_bundle_terms(term%bundle_id)
      IF (SIZE(bundle_terms)<2) RETURN
      IF (id/=bundle_terms(1)) CYCLE
      IF (.NOT.same_batch(bundle_terms,term_group)) RETURN
    END DO
    error=''
  END SUBROUTINE validate_registry

  !> \brief Release every owned descriptor and reset initialization status.
  !> \param[in,out] this Registry to release; repeated calls are safe.
  SUBROUTINE finalize_registry(this)
    CLASS(source_registry_type), INTENT(INOUT) :: this
    IF (ALLOCATED(this%stage_entries)) DEALLOCATE(this%stage_entries)
    IF (ALLOCATED(this%step_entries)) DEALLOCATE(this%step_entries)
    IF (ALLOCATED(this%solver_entries)) DEALLOCATE(this%solver_entries)
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

  !> \brief Return a detached copy of numerical block descriptors.
  !> \param[in] this Registry.
  !> \return Solvers, or an allocated empty list after finalization.
  FUNCTION copy_solvers(this) RESULT(solvers)
    CLASS(source_registry_type), INTENT(IN) :: this
    TYPE(source_solver_type), ALLOCATABLE :: solvers(:)
    ALLOCATE(solvers(0))
    IF (ALLOCATED(this%solver_entries)) solvers=this%solver_entries
  END FUNCTION copy_solvers

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

  !> \brief Select terms for target blocks, EX (0) or mass exchange (-1).
  !> \param[in] active_terms Physical IDs proposed by the configuration bridge.
  !> \param[in] targets Target IDs, zero for EX or minus one for step mass.
  !> \return Stable-order subset; invalid terms reach validation rather than vanish.
  FUNCTION terms_for_targets(active_terms,targets) RESULT(ids)
    INTEGER, INTENT(IN) :: active_terms(:),targets(:)
    INTEGER, ALLOCATABLE :: ids(:)
    TYPE(source_term_type) :: term
    INTEGER :: i,target
    ids=[INTEGER::]
    DO i=1,SIZE(active_terms)
      term=source_term_descriptor(active_terms(i)); target=term%target
      IF (term%role==SOURCE_ROLE_MASS) target=-1
      IF (ANY(targets==target)) ids=[ids,active_terms(i)]
    END DO
  END FUNCTION terms_for_targets

  !> \brief Construct one block with fully allocated arrays.
  !> \param[in] id Unique solver ID.
  !> \param[in] kind Newton or analytic solve kind.
  !> \param[in] targets Declared target blocks.
  !> \param[in] equations Ordered equation footprint.
  !> \param[in] terms Physical term bundle.
  !> \return Descriptor with order one in its legacy-normalized group.
  FUNCTION make_solver(id,kind,targets,equations,terms) RESULT(solver)
    INTEGER, INTENT(IN) :: id,kind,targets(:),equations(:),terms(:)
    TYPE(source_solver_type) :: solver
    solver%id=id; solver%kind=kind; solver%execution_order=1
    solver%target_ids=targets; solver%equation_ids=equations; solver%term_ids=terms
  END FUNCTION make_solver

  !> \brief Construct a fully allocated homogeneous group.
  !> \param[in] id Group identity, independent of order.
  !> \param[in] legacy Original input group, or zero for EX/step.
  !> \param[in] cadence Stage or step namespace.
  !> \param[in] method Group method.
  !> \param[in] phase Evaluation phase within its cadence.
  !> \param[in] order Explicit position within that phase.
  !> \param[in] terms Physical term membership.
  !> \param[in] equations Conservative equation footprint.
  !> \param[in] solves Registry-local block IDs, empty for explicit groups.
  !> \return Group with exactly the solver list required by its method.
  FUNCTION make_group(id,legacy,cadence,method,phase,order,terms,equations,solves) RESULT(group)
    INTEGER, INTENT(IN) :: id,legacy,cadence,method,phase,order,terms(:),equations(:),solves(:)
    TYPE(source_group_type) :: group
    group%id=id; group%legacy_id=legacy; group%cadence=cadence
    group%method=method; group%phase=phase; group%execution_order=order
    group%term_ids=terms; group%equation_ids=equations
    group%implicit_solver_ids=[INTEGER::]; group%analytic_solver_ids=[INTEGER::]
    IF (method==SOURCE_METHOD_NEWTON) group%implicit_solver_ids=solves
    IF (method==SOURCE_METHOD_ANALYTIC) group%analytic_solver_ids=solves
  END FUNCTION make_group

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

  !> \brief Require every member of an atomic kernel batch in one group.
  !> \param[in] ids Required physical batch IDs.
  !> \param[in] term_group Global group position per physical term, zero if absent.
  !> \return True if the whole batch is absent or present in one group.
  LOGICAL FUNCTION same_batch(ids,term_group) RESULT(valid)
    INTEGER, INTENT(IN) :: ids(:),term_group(:)
    valid=.TRUE.
    IF (ALL(term_group(ids)==0)) RETURN
    valid=ALL(term_group(ids)>0) .AND. ALL(term_group(ids)==term_group(ids(1)))
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
