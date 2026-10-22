!********************************************************************************
!> \brief Immutable equation-partition metadata for one configured model
!********************************************************************************

MODULE equation_metadata_2d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: equation_partition_type
  PUBLIC :: IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON

  !> Analytic closed-form stage solve.
  INTEGER, PARAMETER :: IM_GROUP_METHOD_ANALYTIC = 1
  !> Newton stage solve.
  INTEGER, PARAMETER :: IM_GROUP_METHOD_NEWTON = 2

  !> \brief Equation partition and implicit-group metadata for local stage solves.
  !> \details Retains the full explicit/implicit mask and its compact global
  !>          index maps, and assigns each implicit equation to an ordered
  !>          group with a declared solve method and precomputed compact map.
  !>          Explicit equations have group zero; global equation indices and
  !>          the explicit/implicit partition remain unchanged by grouping.
  TYPE :: equation_partition_type

     !> Total number of equations in the configured model.
     INTEGER :: n_equations = 0

     !> Number of equations treated implicitly.
     INTEGER :: n_implicit = 0

     !> Number of equations treated explicitly.
     INTEGER :: n_explicit = 0

     !> Number of nonempty implicit groups, numbered from 1 through n_IM_groups.
     INTEGER :: n_IM_groups = 0

     !> True for equations belonging to the implicit partition; the source of truth for
     !>       explicit/implicit membership.
     LOGICAL, ALLOCATABLE :: implicit(:)

     !> Maps compact implicit indices to full-system equation indices.
     INTEGER, ALLOCATABLE :: implicit_map(:)

     !> Maps compact explicit indices to full-system equation indices.
     INTEGER, ALLOCATABLE :: explicit_map(:)

     !> One group number per equation; zero iff the equation is explicit.
     INTEGER, ALLOCATABLE :: IM_group(:)

     !> Solve method for each group, indexed from 1 through n_IM_groups.
     !> Values are IM_GROUP_METHOD_ANALYTIC or IM_GROUP_METHOD_NEWTON.
     INTEGER, ALLOCATABLE :: IM_group_method(:)

     !> Number of implicit equations in each group.
     INTEGER, ALLOCATABLE :: n_IM_group_equations(:)

     !> Compact global equation indices by group; unused trailing entries are zero.
     INTEGER, ALLOCATABLE :: IM_group_map(:,:)

   CONTAINS

     PROCEDURE :: initialize => initialize_equation_partition
     PROCEDURE :: finalize => finalize_equation_partition
     PROCEDURE :: is_initialized => equation_partition_is_initialized

  END TYPE equation_partition_type

CONTAINS

  !******************************************************************************
  !> \brief Initialize the equation partition and implicit-group metadata.
  !> \details Builds the global explicit/implicit mask and compact index maps.
  !>          If both optional group arguments are omitted, all implicit
  !>          equations are assigned to one NEWTON group (or no group when
  !>          every equation is explicit), preserving the legacy monolithic
  !>          solve contract. If supplied, group numbers must be dense and
  !>          positive for implicit equations, zero for explicit equations,
  !>          and each group's method must be supported. Also precomputes each
  !>          group's compact global equation map.
  !>
  !> \param[in,out] this Partition descriptor to finalize and populate.
  !> \param[in] implicit_mask One flag per full-system equation; true selects implicit treatment.
  !> \param[in] IM_group Optional group number per equation; zero selects the explicit partition.
  !> \param[in] IM_group_method Optional method identifier per group, indexed by group number.
  !> \note IM_group and IM_group_method must either both be present or both be omitted.
  !>       Invalid dimensions, membership, group numbering, or methods terminate with ERROR STOP.
  !******************************************************************************

  SUBROUTINE initialize_equation_partition( this, implicit_mask, IM_group, IM_group_method )

    CLASS(equation_partition_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: implicit_mask(:)
    INTEGER, INTENT(IN), OPTIONAL :: IM_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: IM_group_method(:)

    INTEGER :: equation
    INTEGER :: group
    INTEGER :: i_explicit
    INTEGER :: i_implicit
    INTEGER :: i_group

    CALL this%finalize

    IF ( PRESENT(IM_group) .NEQV. PRESENT(IM_group_method) ) THEN
       ERROR STOP 'initialize_equation_partition: both group arguments are required together'
    END IF

    this%n_equations = SIZE(implicit_mask)
    this%n_implicit = COUNT(implicit_mask)
    this%n_explicit = this%n_equations - this%n_implicit

    ALLOCATE( this%implicit(this%n_equations) )
    ALLOCATE( this%implicit_map(this%n_implicit) )
    ALLOCATE( this%explicit_map(this%n_explicit) )
    ALLOCATE( this%IM_group(this%n_equations) )

    this%implicit = implicit_mask

    i_implicit = 0
    i_explicit = 0

    DO equation = 1, this%n_equations

       IF ( this%implicit(equation) ) THEN
          i_implicit = i_implicit + 1
          this%implicit_map(i_implicit) = equation
       ELSE
          i_explicit = i_explicit + 1
          this%explicit_map(i_explicit) = equation
       END IF

    END DO

    IF ( PRESENT(IM_group) ) THEN

       IF ( SIZE(IM_group) .NE. this%n_equations ) THEN
          ERROR STOP 'initialize_equation_partition: group count must match equation count'
       END IF

       this%n_IM_groups = SIZE(IM_group_method)
       ALLOCATE( this%IM_group_method(this%n_IM_groups) )
       this%IM_group = IM_group
       this%IM_group_method = IM_group_method

    ELSEIF ( this%n_implicit .GT. 0 ) THEN

       this%n_IM_groups = 1
       ALLOCATE( this%IM_group_method(this%n_IM_groups) )
       this%IM_group = 0
       this%IM_group(this%implicit_map) = 1
       this%IM_group_method(1) = IM_GROUP_METHOD_NEWTON

    ELSE

       this%n_IM_groups = 0
       ALLOCATE( this%IM_group_method(0) )
       this%IM_group = 0

    END IF

    ALLOCATE( this%n_IM_group_equations(this%n_IM_groups) )
    ALLOCATE( this%IM_group_map(this%n_implicit,this%n_IM_groups) )
    this%n_IM_group_equations = 0
    this%IM_group_map = 0

    DO equation = 1, this%n_equations
       group = this%IM_group(equation)
       IF ( group .EQ. 0 ) CYCLE
       i_group = this%n_IM_group_equations(group) + 1
       this%n_IM_group_equations(group) = i_group
       this%IM_group_map(i_group,group) = equation
    END DO

    IF ( .NOT. this%is_initialized() ) THEN
       ERROR STOP 'initialize_equation_partition: invalid equation-group metadata'
    END IF

  END SUBROUTINE initialize_equation_partition


  !******************************************************************************
  !> \brief Release and reset all equation-partition and group metadata.
  !>
  !> \param[in,out] this Descriptor whose allocated masks, maps and group metadata are released.
  !******************************************************************************

  SUBROUTINE finalize_equation_partition( this )

    CLASS(equation_partition_type), INTENT(INOUT) :: this

    IF ( ALLOCATED(this%implicit) ) DEALLOCATE(this%implicit)
    IF ( ALLOCATED(this%implicit_map) ) DEALLOCATE(this%implicit_map)
    IF ( ALLOCATED(this%explicit_map) ) DEALLOCATE(this%explicit_map)
    IF ( ALLOCATED(this%IM_group) ) DEALLOCATE(this%IM_group)
    IF ( ALLOCATED(this%IM_group_method) ) DEALLOCATE(this%IM_group_method)
    IF ( ALLOCATED(this%n_IM_group_equations) )                             &
         DEALLOCATE(this%n_IM_group_equations)
    IF ( ALLOCATED(this%IM_group_map) ) DEALLOCATE(this%IM_group_map)

    this%n_equations = 0
    this%n_implicit = 0
    this%n_explicit = 0
    this%n_IM_groups = 0

  END SUBROUTINE finalize_equation_partition


  !******************************************************************************
  !> \brief Check the consistency of the equation partition and implicit groups.
  !>
  !> \param[in] this Partition descriptor to validate.
  !> \return True when counts, masks, compact maps, group membership and methods are consistent.
  !******************************************************************************

  LOGICAL FUNCTION equation_partition_is_initialized( this )

    CLASS(equation_partition_type), INTENT(IN) :: this

    INTEGER :: equation
    INTEGER :: group
    INTEGER :: i_implicit
    INTEGER :: i_explicit
    INTEGER :: i_group

    equation_partition_is_initialized = .FALSE.

    IF ( .NOT. ALLOCATED(this%implicit) ) RETURN
    IF ( .NOT. ALLOCATED(this%implicit_map) ) RETURN
    IF ( .NOT. ALLOCATED(this%explicit_map) ) RETURN
    IF ( .NOT. ALLOCATED(this%IM_group) ) RETURN
    IF ( .NOT. ALLOCATED(this%IM_group_method) ) RETURN
    IF ( .NOT. ALLOCATED(this%n_IM_group_equations) ) RETURN
    IF ( .NOT. ALLOCATED(this%IM_group_map) ) RETURN

    IF ( this%n_equations .LT. 0 .OR. this%n_implicit .LT. 0 .OR.          &
         this%n_explicit .LT. 0 .OR. this%n_IM_groups .LT. 0 ) RETURN
    IF ( this%n_implicit + this%n_explicit .NE. this%n_equations ) RETURN

    IF ( this%n_equations .NE. SIZE(this%implicit) .OR.                     &
         this%n_implicit .NE. SIZE(this%implicit_map) .OR.                   &
         this%n_explicit .NE. SIZE(this%explicit_map) .OR.                   &
         this%n_equations .NE. SIZE(this%IM_group) .OR.                      &
         this%n_IM_groups .NE. SIZE(this%IM_group_method) .OR.               &
         this%n_IM_groups .NE. SIZE(this%n_IM_group_equations) .OR.          &
         SIZE(this%IM_group_map,1) .NE. this%n_implicit .OR.                  &
         SIZE(this%IM_group_map,2) .NE. this%n_IM_groups ) RETURN
    IF ( COUNT(this%implicit) .NE. this%n_implicit ) RETURN

    IF ( ANY(this%implicit .NEQV. (this%IM_group .GT. 0)) ) RETURN
    IF ( ANY(this%IM_group .LT. 0) .OR.                                      &
         ANY(this%IM_group .GT. this%n_IM_groups) ) RETURN
    IF ( ANY(this%IM_group_method .NE. IM_GROUP_METHOD_ANALYTIC .AND.        &
         this%IM_group_method .NE. IM_GROUP_METHOD_NEWTON) ) RETURN

    i_implicit = 0
    i_explicit = 0
    DO equation = 1, this%n_equations
       IF ( this%implicit(equation) ) THEN
          i_implicit = i_implicit + 1
          IF ( this%implicit_map(i_implicit) .NE. equation ) RETURN
       ELSE
          i_explicit = i_explicit + 1
          IF ( this%explicit_map(i_explicit) .NE. equation ) RETURN
       END IF
    END DO

    DO group = 1, this%n_IM_groups
       IF ( .NOT. ANY(this%IM_group .EQ. group) ) RETURN
       IF ( this%n_IM_group_equations(group) .LT. 1 .OR.                    &
            this%n_IM_group_equations(group) .GT. this%n_implicit ) RETURN

       i_group = 0
       DO i_implicit = 1, this%n_implicit
          equation = this%implicit_map(i_implicit)
          IF ( this%IM_group(equation) .NE. group ) CYCLE
          i_group = i_group + 1
          IF ( i_group .GT. this%n_IM_group_equations(group) ) RETURN
          IF ( this%IM_group_map(i_group,group) .NE. equation ) RETURN
       END DO

       IF ( i_group .NE. this%n_IM_group_equations(group) ) RETURN
       IF ( ANY(this%IM_group_map(i_group+1:this%n_implicit,group) .NE. 0) ) &
            RETURN
    END DO

    equation_partition_is_initialized = .TRUE.

  END FUNCTION equation_partition_is_initialized

END MODULE equation_metadata_2d
