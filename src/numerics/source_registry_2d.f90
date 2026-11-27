!********************************************************************************
!> \brief Validate and own the ordered Lie groups selected by physical configuration.
!> \details source_configuration_2d assembles candidate groups from the selected rheology,
!>          thermal, pore-pressure, injection and mass-exchange settings. This module checks
!>          those candidates against source_catalog_2d and stores an immutable-in-practice
!>          schedule for time_integration_2d and lie_source_executor_2d.
!>
!>          Each group has a Lie ordering ID, a scope (L1, L2 or both layers), an
!>          integration method, canonical source names and a conservative equation set.
!>          The registry rejects unsupported methods, duplicate source identities,
!>          invalid group ordering and incompatible Newton equation footprints. Interlayer
!>          mass exchange and drag have their own pair-wide rules. No physical rates
!>          or evolving flow states are computed or stored here.
!>
!>          Validation is transactional: an unsuccessful initialization leaves the registry
!>          uninitialized. The operational executor receives a detached copy of the groups,
!>          so worker code does not mutate the registry's owned metadata.
!********************************************************************************
MODULE source_registry_2d
    USE source_catalog_2d
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: source_group_type, source_registry_type
    PUBLIC :: group_contains_name
    PUBLIC :: SOURCE_SCOPE_L1, SOURCE_SCOPE_L2, SOURCE_SCOPE_INTERLAYER

    INTEGER, PARAMETER :: SOURCE_SCOPE_L1=1, SOURCE_SCOPE_L2=2, SOURCE_SCOPE_INTERLAYER=3

    !> \brief Ordered, scoped description of one physical-source integration group.
    !> \details GROUP selects the position in the Lie source sequence, while method
    !>          selects explicit evaluation, a coupled Newton solve or an analytic substep.
    !>          source_names enumerate the configured laws; equation_ids specify the
    !>          conservative subset for Newton, the full set for explicit groups and the
    !>          empty set for analytic groups. scope distinguishes L1, L2 and paired updates.
    !>          Resolved IDs and the bit mask accelerate execution without re-parsing names.
    TYPE :: source_group_type
        INTEGER :: id=0 !< Public or normalized internal group identity.
        INTEGER :: method=0 !< Explicit, Newton or analytic.
        INTEGER :: scope=SOURCE_SCOPE_L1 !< Source ownership: lower, upper or interlayer.
        CHARACTER(LEN=64), ALLOCATABLE :: source_names(:) !< Unique schedulable source identities.
        INTEGER, ALLOCATABLE :: equation_ids(:) !< Full structurally affected Newton equation set or explicit footprint.
        LOGICAL :: rheology_implicit=.FALSE. !< Cached Newton implicit-rheology role.
        LOGICAL :: thermal_active=.FALSE. !< Cached thermal source presence.
        INTEGER :: rheology_id=RHEOLOGY_ID_NONE !< Cached numeric source dispatch.
        INTEGER :: source_mask=0 !< Immutable bitset of validated catalogue IDs.
    END TYPE source_group_type

    !> \brief Single owner of the validated group schedule for one simulation.
    !> \details Initialized during simulation setup and finalized with its context.
    !>          Its ordered descriptors are private; groups() supplies a detached
    !>          copy to the operational executor. The cached implicit-source flag
    !>          reports whether any validated group requires Newton integration.
    TYPE :: source_registry_type
    PRIVATE
        LOGICAL :: initialized=.FALSE. !< Cached successful-validation status.
        LOGICAL :: implicit_sources_flag=.FALSE. !< Any validated NEWTON source group.
        !> Validated and ordered Lie group metadata owned by this registry.
        TYPE(source_group_type), ALLOCATABLE :: lie_entries(:)
    CONTAINS

        PROCEDURE :: groups => copy_lie_groups
        PROCEDURE :: initialize => initialize_lie_registry
        PROCEDURE :: finalize => finalize_registry
        PROCEDURE :: is_initialized => registry_is_initialized
        PROCEDURE :: has_implicit_sources => registry_has_implicit_sources
    END TYPE source_registry_type

CONTAINS

    !> @brief Validate and store groups for the first-order Lie integrator.
    !> @details This setup-only validation starts by clearing any previous schedule.
    !>          It requires an explicit group zero, sorted group IDs/scopes and source
    !>          identities present in SOURCE_CATALOG. Method permissions, allowed
    !>          source combinations and per-layer uniqueness are checked before storing
    !>          descriptors. Newton equation_ids must equal the union of the structural
    !>          equation footprints of all their sources; analytic groups are exclusive.
    !>          Shared mass-exchange inventory sources must use the same explicit group.
    !>          Interlayer updates may only contain the supported paired mass or drag law.
    !>          A failure is returned through error and no initialized schedule is exposed.
    !> @param[in,out] this Registry receiving validated group metadata.
    !> @param[in] n_eqns Number of conservative equations.
    !> @param[in] idx_poreEqn Pressure-of-pore equation index, zero when absent.
    !> @param[in] groups Ordered source group candidates from input configuration.
    !> @param[out] error Validation diagnostic; empty on success.
    SUBROUTINE initialize_lie_registry(this,n_eqns,idx_poreEqn,groups,error)
        CLASS(source_registry_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: n_eqns,idx_poreEqn
        TYPE(source_group_type), INTENT(IN) :: groups(:)
        CHARACTER(*), INTENT(OUT) :: error
        TYPE(source_term_type) :: term
        CHARACTER(LEN=64), ALLOCATABLE :: seen(:)
        INTEGER, ALLOCATABLE :: seen_scope(:)
        INTEGER :: i,j,mass_group(2),k,eq,rheology_count(2),source_id,layer
        INTEGER, ALLOCATABLE :: expected(:),affected(:)
        CALL this%finalize()
        mass_group=-1
        rheology_count=0
        ALLOCATE(seen(0),seen_scope(0))
        error='D19: invalid Lie group configuration'
        IF (n_eqns<4 .OR. SIZE(groups)==0) RETURN

        IF (idx_poreEqn/=0) THEN
            IF (idx_poreEqn<5 .OR. idx_poreEqn>n_eqns) RETURN

        END IF

        IF (groups(1)%id/=0 .OR. groups(1)%method/=SOURCE_METHOD_EXPLICIT) RETURN

        DO i=1,SIZE(groups)

            IF (groups(i)%scope<SOURCE_SCOPE_L1 .OR. groups(i)%scope>SOURCE_SCOPE_INTERLAYER) RETURN
            ! Interlayer groups may use explicit mass transfer or analytic drag,
            ! each in an independently scheduled, single-source group.
            layer=groups(i)%scope
            IF (i>1) THEN
                IF (groups(i)%id<groups(i-1)%id) RETURN
                IF (groups(i)%id==groups(i-1)%id .AND. groups(i)%scope<=groups(i-1)%scope) RETURN
            END IF
            IF (groups(i)%id==0 .AND. groups(i)%method/=SOURCE_METHOD_EXPLICIT) RETURN

            IF (.NOT.ALLOCATED(groups(i)%source_names)) RETURN
            IF (groups(i)%id/=0 .AND. SIZE(groups(i)%source_names)==0) RETURN
            IF (.NOT.ALLOCATED(groups(i)%equation_ids)) RETURN
            IF (.NOT.valid_equations(groups(i)%equation_ids,n_eqns)) RETURN

            SELECT CASE(groups(i)%method)
            CASE(SOURCE_METHOD_EXPLICIT)
                IF (SIZE(groups(i)%equation_ids)/=n_eqns) RETURN
            CASE(SOURCE_METHOD_NEWTON)
                IF (SIZE(groups(i)%equation_ids)==0) RETURN
            CASE(SOURCE_METHOD_ANALYTIC)
                IF (SIZE(groups(i)%source_names)/=1) RETURN
                IF (SIZE(groups(i)%equation_ids)/=0) RETURN
            CASE DEFAULT
                RETURN

            END SELECT

            ALLOCATE(expected(0))

            DO j=1,SIZE(groups(i)%source_names)
                term=source_term_descriptor_by_name(groups(i)%source_names(j))
                IF (LEN_TRIM(term%name)==0) RETURN
                IF (groups(i)%scope==SOURCE_SCOPE_INTERLAYER) THEN
                    SELECT CASE(TRIM(term%name))
                    CASE('interlayer_mass_exchange')
                        IF (groups(i)%method/=SOURCE_METHOD_EXPLICIT) THEN
                            error='TL07: INTERLAYER mass exchange requires EXPLICIT'
                            RETURN
                        END IF
                    CASE('interlayer_drag')
                        IF (groups(i)%method/=SOURCE_METHOD_ANALYTIC .OR. groups(i)%id==0) THEN
                            error='TL09: interlayer drag requires a positive ANALYTIC group'
                            RETURN
                        END IF
                    CASE DEFAULT
                        error='TL07: only paired mass exchange or drag is supported in INTERLAYER'
                        RETURN
                    END SELECT
                ELSEIF (term%name=='interlayer_mass_exchange' .OR. term%name=='interlayer_drag') THEN
                    error='TL07: paired interlayer source requires INTERLAYER scope'
                    RETURN
                END IF

                IF (INDEX(TRIM(term%name),'mass_exchange_')==1) THEN
                    ! Inventory-coupled mass exchange currently supports only an
                    ! explicit Lie group; reject incompatible configurations here.
                    IF (groups(i)%method/=SOURCE_METHOD_EXPLICIT) RETURN
                    IF (mass_group(layer)>=0 .AND. mass_group(layer)/=groups(i)%id) RETURN
                    mass_group(layer)=groups(i)%id

                END IF

                IF (INDEX(TRIM(term%name),'rheology_')==1) THEN
                    rheology_count(layer)=rheology_count(layer)+1
                    IF (rheology_count(layer)>1) RETURN

                END IF

                IF (ANY((seen==term%name) .AND. (seen_scope==layer))) RETURN
                seen=[CHARACTER(LEN=64) :: seen,term%name]
                seen_scope=[seen_scope,layer]

                IF (groups(i)%method==SOURCE_METHOD_NEWTON) THEN
                    affected=source_newton_affected_equations(term%name)
                    IF (SIZE(affected)==0) RETURN

                    DO k=1,SIZE(affected)
                        eq=affected(k)
                        IF (.NOT.ANY(expected==eq)) expected=[expected,eq]

                        END DO

                END IF

                SELECT CASE(groups(i)%method)
                CASE(SOURCE_METHOD_EXPLICIT)
                    IF (.NOT.term%has_explicit) RETURN
                CASE(SOURCE_METHOD_NEWTON)
                    IF (.NOT.term%has_newton) RETURN
                CASE(SOURCE_METHOD_ANALYTIC)
                    IF (.NOT.term%has_analytic) RETURN

                END SELECT

            END DO

            IF (groups(i)%method==SOURCE_METHOD_NEWTON) THEN
                IF (.NOT.same_set(groups(i)%equation_ids,expected)) RETURN

            END IF

            DEALLOCATE(expected)

        END DO
        ! A vaporization-only group would silently evaluate a zero deposition rate.
        DO layer=SOURCE_SCOPE_L1,SOURCE_SCOPE_L2
            IF (ANY((seen=='mass_exchange_vaporization') .AND. (seen_scope==layer)) .AND. &
                .NOT.ANY((seen=='mass_exchange_deposition') .AND. (seen_scope==layer))) RETURN
        END DO
        this%lie_entries=groups
        ! Resolve source identities and capabilities once, outside all cell loops.
        DO i=1,SIZE(this%lie_entries)
            this%lie_entries(i)%rheology_id=RHEOLOGY_ID_NONE
            this%lie_entries(i)%rheology_implicit=.FALSE.
            this%lie_entries(i)%source_mask=0
            this%lie_entries(i)%thermal_active=.FALSE.

            DO j=1,SIZE(this%lie_entries(i)%source_names)
                source_id=source_term_id_by_name(this%lie_entries(i)%source_names(j))
                IF (source_id<=0) RETURN
                this%lie_entries(i)%source_mask=IOR(this%lie_entries(i)%source_mask,source_id_mask(source_id))
                IF (INDEX(TRIM(this%lie_entries(i)%source_names(j)),'rheology_')==1) THEN
                    this%lie_entries(i)%rheology_id=rheology_source_id(this%lie_entries(i)%source_names(j))
                    this%lie_entries(i)%rheology_implicit=                                                    &
                        rheology_name_has_role(this%lie_entries(i)%source_names(j),SOURCE_ROLE_IMPLICIT)
                ELSEIF (this%lie_entries(i)%source_names(j)=='thermal_cooling') THEN
                    this%lie_entries(i)%thermal_active=.TRUE.
                END IF
            END DO
        END DO
        this%implicit_sources_flag=ANY(groups%method==SOURCE_METHOD_NEWTON)
        this%initialized=.TRUE.
        error=''
    END SUBROUTINE initialize_lie_registry

    !> \brief Map one target to conservative equation indices.
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

    !> \brief Compare unique ID sets without changing their recorded order.
    !> \param[in] left First list.
    !> \param[in] right Second list.
    !> \return True for equal-sized unique lists with identical membership.
    LOGICAL FUNCTION same_set(left,right) RESULT(equal)
        INTEGER, INTENT(IN) :: left(:),right(:)
        INTEGER :: i
        equal=.FALSE.
        IF (SIZE(left)/=SIZE(right)) RETURN

        DO i=1,SIZE(left)
        IF (.NOT.ANY(right==left(i))) RETURN

        END DO

        DO i=1,SIZE(left)
        IF (COUNT(left==left(i))/=1 .OR. COUNT(right==left(i))/=1) RETURN

        END DO

        equal=.TRUE.
    END FUNCTION same_set

    !> \brief Release every owned descriptor and reset initialization status.
    !> \param[in,out] this Registry to release; repeated calls are safe.
    SUBROUTINE finalize_registry(this)
        CLASS(source_registry_type), INTENT(INOUT) :: this
        IF (ALLOCATED(this%lie_entries)) DEALLOCATE(this%lie_entries)
        this%initialized=.FALSE.
        this%implicit_sources_flag=.FALSE.
    END SUBROUTINE finalize_registry

    !> @brief Test whether a group schedules the named physical source.
    !> @param[in] group Lie group and its configured canonical identities.
    !> @param[in] name Canonical identity to search for.
    !> @return True if the name occurs in the allocated group source list.
    PURE LOGICAL FUNCTION group_contains_name(group,name) RESULT(found)
        TYPE(source_group_type), INTENT(IN) :: group
        CHARACTER(*), INTENT(IN) :: name
        found=.FALSE.
        IF (.NOT.ALLOCATED(group%source_names)) RETURN
        found=ANY(group%source_names==TRIM(name))
    END FUNCTION group_contains_name


    !> \brief Query cached initialization status without walking term arrays.
    !> \param[in] this Registry.
    !> \return True only after successful transactional validation.
    !> Cached query, set once following successful registry validation.
    PURE LOGICAL FUNCTION registry_has_implicit_sources(this)
        CLASS(source_registry_type), INTENT(IN) :: this
        registry_has_implicit_sources = this%initialized .AND. this%implicit_sources_flag
    END FUNCTION registry_has_implicit_sources

    !> @brief Query whether group metadata passed initialization and validation.
    !> @param[in] this Lie registry to examine.
    !> @return True only after successful initialization.
    PURE LOGICAL FUNCTION registry_is_initialized(this)
        CLASS(source_registry_type), INTENT(IN) :: this
        registry_is_initialized=this%initialized
    END FUNCTION registry_is_initialized

    !> @brief Return a detached copy of the ordered Lie source groups.
    !> @details The executor takes ownership of its copy, including allocatable name and
    !>          equation-index arrays. This prevents runtime bookkeeping from changing
    !>          the registry's validated description of the simulation.
    !> @param[in] this Initialized source registry.
    !> @return Independently allocated group descriptors for the caller.
    FUNCTION copy_lie_groups(this) RESULT(groups)
        CLASS(source_registry_type), INTENT(IN) :: this
        TYPE(source_group_type), ALLOCATABLE :: groups(:)
        ALLOCATE(groups(0))
        IF (ALLOCATED(this%lie_entries)) groups=this%lie_entries
    END FUNCTION copy_lie_groups

END MODULE source_registry_2d
