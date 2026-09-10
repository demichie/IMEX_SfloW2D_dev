!********************************************************************************
!> \brief Immutable equation-partition metadata for one configured model
!********************************************************************************
MODULE equation_metadata_2d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: equation_partition_type

  TYPE :: equation_partition_type

     !> Total number of equations in the configured model
     INTEGER :: n_equations = 0

     !> Number of equations treated implicitly
     INTEGER :: n_implicit = 0

     !> Number of equations treated explicitly
     INTEGER :: n_explicit = 0

     !> True for equations belonging to the implicit partition
     LOGICAL, ALLOCATABLE :: implicit(:)

     !> Maps compact implicit indices to full-system equation indices
     INTEGER, ALLOCATABLE :: implicit_map(:)

     !> Maps compact explicit indices to full-system equation indices
     INTEGER, ALLOCATABLE :: explicit_map(:)

   CONTAINS

     PROCEDURE :: initialize => initialize_equation_partition
     PROCEDURE :: finalize => finalize_equation_partition
     PROCEDURE :: is_initialized => equation_partition_is_initialized

  END TYPE equation_partition_type

CONTAINS

  !******************************************************************************
  !> \brief Configure the immutable explicit/implicit equation partition
  !******************************************************************************
  SUBROUTINE initialize_equation_partition( this, implicit_mask )

    CLASS(equation_partition_type), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: implicit_mask(:)

    INTEGER :: equation
    INTEGER :: i_explicit
    INTEGER :: i_implicit

    CALL this%finalize

    this%n_equations = SIZE(implicit_mask)
    this%n_implicit = COUNT(implicit_mask)
    this%n_explicit = this%n_equations - this%n_implicit

    ALLOCATE( this%implicit(this%n_equations) )
    ALLOCATE( this%implicit_map(this%n_implicit) )
    ALLOCATE( this%explicit_map(this%n_explicit) )

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

  END SUBROUTINE initialize_equation_partition


  !******************************************************************************
  !> \brief Release equation-partition storage
  !******************************************************************************
  SUBROUTINE finalize_equation_partition( this )

    CLASS(equation_partition_type), INTENT(INOUT) :: this

    IF ( ALLOCATED(this%implicit) ) DEALLOCATE(this%implicit)
    IF ( ALLOCATED(this%implicit_map) ) DEALLOCATE(this%implicit_map)
    IF ( ALLOCATED(this%explicit_map) ) DEALLOCATE(this%explicit_map)

    this%n_equations = 0
    this%n_implicit = 0
    this%n_explicit = 0

  END SUBROUTINE finalize_equation_partition


  !******************************************************************************
  !> \brief Report whether the partition has been configured
  !******************************************************************************
  LOGICAL FUNCTION equation_partition_is_initialized( this )

    CLASS(equation_partition_type), INTENT(IN) :: this

    equation_partition_is_initialized = .FALSE.

    IF ( .NOT. ALLOCATED(this%implicit) ) RETURN
    IF ( .NOT. ALLOCATED(this%implicit_map) ) RETURN
    IF ( .NOT. ALLOCATED(this%explicit_map) ) RETURN

    equation_partition_is_initialized =                                     &
         this%n_equations .EQ. SIZE(this%implicit)                           &
         .AND. this%n_implicit .EQ. SIZE(this%implicit_map)                  &
         .AND. this%n_explicit .EQ. SIZE(this%explicit_map)

  END FUNCTION equation_partition_is_initialized

END MODULE equation_metadata_2d
