!********************************************************************************
!> \brief Validated source groups derived from physical source configuration.
!> \details Builds metadata once, outside cell/stage loops. Does not evaluate
!>          physics or depend on legacy global equation partitions.
!>          Private storage is exposed only through copies; failed initialization
!>          releases all storage. Joint Newton blocks remain coupled.
!********************************************************************************
MODULE source_registry_2d
    USE source_catalog_2d
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: source_group_type, source_registry_type
    PUBLIC :: group_contains_name, group_source_names

    !> \brief Ordered Lie source group; id=0 is the explicit spatial/source group.
    TYPE :: source_group_type
        INTEGER :: id=0 !< Public or normalized internal group identity.
        INTEGER :: method=0 !< Explicit, Newton or analytic.
        CHARACTER(LEN=64), ALLOCATABLE :: source_names(:) !< Unique schedulable source identities.
        INTEGER, ALLOCATABLE :: equation_ids(:) !< Full structurally affected Newton equation set or explicit footprint.
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

    !> @brief Validate and store groups for the first-order Lie integrator.
    !> @details Rejects inconsistent methods, duplicated identities and invalid
    !>          equation footprints; updates registry ownership on successful validation.
    !> @param[in,out] this Registry receiving validated group metadata.
    !> @param[in] n_equations Number of conservative equations.
    !> @param[in] pore_equation Pressure-of-pore equation index, zero when absent.
    !> @param[in] groups Ordered source group candidates from input configuration.
    !> @param[out] error Validation diagnostic; empty on success.
    SUBROUTINE initialize_lie_registry(this,n_equations,pore_equation,groups,error)
        CLASS(source_registry_type), INTENT(INOUT) :: this
        INTEGER, INTENT(IN) :: n_equations,pore_equation
        TYPE(source_group_type), INTENT(IN) :: groups(:)
        CHARACTER(*), INTENT(OUT) :: error
        TYPE(source_term_type) :: term
        CHARACTER(LEN=64), ALLOCATABLE :: seen(:)
        INTEGER :: i,j,mass_group,k,eq,rheology_count
        INTEGER, ALLOCATABLE :: expected(:),affected(:)
        CALL this%finalize()
        mass_group=-1
        rheology_count=0
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
            CASE(SOURCE_METHOD_EXPLICIT)
                IF (SIZE(groups(i)%equation_ids)/=n_equations) RETURN
            CASE(SOURCE_METHOD_NEWTON)
                IF (SIZE(groups(i)%equation_ids)==0) RETURN
            CASE(SOURCE_METHOD_ANALYTIC)
                IF (SIZE(groups(i)%source_names)/=1) RETURN
                IF (SIZE(groups(i)%equation_ids)/=0) RETURN
            CASE DEFAULT
                RETURN

            END SELECT

            IF (groups(i)%method==SOURCE_METHOD_NEWTON .AND. SIZE(groups(i)%equation_ids)==0) RETURN
            IF (groups(i)%method==SOURCE_METHOD_ANALYTIC .AND. SIZE(groups(i)%source_names)/=1) RETURN
            ALLOCATE(expected(0))

            DO j=1,SIZE(groups(i)%source_names)
                term=source_term_descriptor_by_name(groups(i)%source_names(j))
                IF (LEN_TRIM(term%name)==0) RETURN

                IF (INDEX(TRIM(term%name),'mass_exchange_')==1) THEN
                    ! Inventory-coupled mass exchange currently supports only an
                    ! explicit Lie group; reject incompatible configurations here.
                    IF (groups(i)%method/=SOURCE_METHOD_EXPLICIT) RETURN
                    IF (mass_group>=0 .AND. mass_group/=groups(i)%id) RETURN
                    mass_group=groups(i)%id

                END IF

                IF (INDEX(TRIM(term%name),'rheology_')==1) THEN
                    rheology_count=rheology_count+1
                    IF (rheology_count>1) RETURN

                END IF

                IF (ANY(seen==term%name)) RETURN
                seen=[CHARACTER(LEN=64) :: seen,term%name]

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
        IF (ANY(seen=='mass_exchange_vaporization') .AND.                                                              &
            .NOT.ANY(seen=='mass_exchange_deposition')) RETURN
        this%n_equations=n_equations
        this%pore_equation=pore_equation
        this%lie_entries=groups
        this%implicit_sources_flag=ANY(groups%method==SOURCE_METHOD_NEWTON)
        this%initialized=.TRUE.
        error=''
    END SUBROUTINE initialize_lie_registry

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
        this%n_equations=0; this%pore_equation=0; this%initialized=.FALSE.
        this%implicit_sources_flag=.FALSE.
    END SUBROUTINE finalize_registry
        PROCEDURE :: groups => copy_lie_groups
        PROCEDURE :: initialize => initialize_lie_registry
        PROCEDURE :: finalize => finalize_registry
        PROCEDURE :: is_initialized => registry_is_initialized
        PROCEDURE :: has_implicit_sources => registry_has_implicit_sources
    END TYPE source_registry_type
CONTAINS
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

    !> @brief Copy the source identities owned by a Lie group.
    !> @param[in] group Lie group with an optional allocated source list.
    !> @return Independent source name array, empty if none are configured.
    FUNCTION group_source_names(group) RESULT(names)
        TYPE(source_group_type), INTENT(IN) :: group
        CHARACTER(LEN=64), ALLOCATABLE :: names(:)
        ALLOCATE(names(0))
        IF (ALLOCATED(group%source_names)) names=group%source_names
    END FUNCTION group_source_names

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
    !> @param[in] this Initialized source registry.
    !> @return Independently allocated group descriptors for the caller.
    FUNCTION copy_lie_groups(this) RESULT(groups)
        CLASS(source_registry_type), INTENT(IN) :: this
        TYPE(source_group_type), ALLOCATABLE :: groups(:)
        ALLOCATE(groups(0))
        IF (ALLOCATED(this%lie_entries)) groups=this%lie_entries
    END FUNCTION copy_lie_groups

END MODULE source_registry_2d
