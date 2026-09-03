!********************************************************************************
!> \brief Active computational-domain management
!>
!> This module owns the masks and compact index lists that identify cells and
!> interfaces requiring numerical work.
!********************************************************************************
MODULE domain_2d

  USE parameters_2d, ONLY : wp, n_RK
  USE parameters_2d, ONLY : radial_source_flag, bottom_radial_source_flag

  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, comp_cells_xy
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : source_cell, cell_source_fractions

  IMPLICIT NONE

  PRIVATE

  TYPE, PUBLIC :: domain_type

     REAL(wp), ALLOCATABLE :: solve_mask_time(:,:)

     LOGICAL, ALLOCATABLE :: solve_mask(:,:)
     LOGICAL, ALLOCATABLE :: solve_mask_temp(:,:)
     LOGICAL, ALLOCATABLE :: solve_mask_x(:,:)
     LOGICAL, ALLOCATABLE :: solve_mask_y(:,:)

     INTEGER :: solve_cells
     INTEGER :: solve_interfaces_x
     INTEGER :: solve_interfaces_y

     INTEGER, ALLOCATABLE :: j_cent(:)
     INTEGER, ALLOCATABLE :: k_cent(:)
     INTEGER, ALLOCATABLE :: j_stag_x(:)
     INTEGER, ALLOCATABLE :: k_stag_x(:)
     INTEGER, ALLOCATABLE :: j_stag_y(:)
     INTEGER, ALLOCATABLE :: k_stag_y(:)

   CONTAINS

     PROCEDURE :: initialize => initialize_domain
     PROCEDURE :: finalize => finalize_domain
     PROCEDURE :: check_solve

  END TYPE domain_type

  TYPE(domain_type), PUBLIC, TARGET :: domain

CONTAINS

  SUBROUTINE initialize_domain(this)

    CLASS(domain_type), INTENT(INOUT) :: this

    ALLOCATE( this%solve_mask_time(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%solve_mask(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%solve_mask_temp(comp_cells_x,comp_cells_y) )
    ALLOCATE( this%solve_mask_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( this%solve_mask_y(comp_cells_x,comp_interfaces_y) )

    this%solve_mask_time = 0.0_wp

    this%solve_mask(1,1:comp_cells_y) = .TRUE.
    this%solve_mask(comp_cells_x,1:comp_cells_y) = .TRUE.
    this%solve_mask(1:comp_cells_x,1) = .TRUE.
    this%solve_mask(1:comp_cells_x,comp_cells_y) = .TRUE.

    comp_cells_xy = comp_cells_x * comp_cells_y

    ALLOCATE( this%j_cent(comp_cells_xy) )
    ALLOCATE( this%k_cent(comp_cells_xy) )
    ALLOCATE( this%j_stag_x(comp_interfaces_x*comp_cells_y) )
    ALLOCATE( this%k_stag_x(comp_interfaces_x*comp_cells_y) )
    ALLOCATE( this%j_stag_y(comp_cells_x*comp_interfaces_y) )
    ALLOCATE( this%k_stag_y(comp_cells_x*comp_interfaces_y) )

  END SUBROUTINE initialize_domain

  SUBROUTINE finalize_domain(this)

    CLASS(domain_type), INTENT(INOUT) :: this

    DEALLOCATE( this%solve_mask_time )
    DEALLOCATE( this%solve_mask )
    DEALLOCATE( this%solve_mask_temp )
    DEALLOCATE( this%solve_mask_x )
    DEALLOCATE( this%solve_mask_y )

    DEALLOCATE( this%j_cent, this%k_cent )
    DEALLOCATE( this%j_stag_x, this%k_stag_x )
    DEALLOCATE( this%j_stag_y, this%k_stag_y )

  END SUBROUTINE finalize_domain

  SUBROUTINE check_solve(this, q, t, solve_all)

    IMPLICIT NONE

    CLASS(domain_type), INTENT(INOUT) :: this
    REAL(wp), INTENT(IN) :: q(:,:,:)
    REAL(wp), INTENT(IN) :: t
    LOGICAL, INTENT(IN) :: solve_all
    
    INTEGER :: i,j,k

    !$OMP PARALLEL
         
    IF ( solve_all ) THEN

       !$OMP WORKSHARE
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE.
       !$OMP END WORKSHARE

    ELSE
       
       !$OMP WORKSHARE
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .FALSE.
       !$OMP END WORKSHARE

    END IF
    !$OMP BARRIER
    
    !$OMP WORKSHARE
    WHERE ( ( q(1,2:comp_cells_x-1,2:comp_cells_y-1) .GT. 0.0_wp ) .AND.        &
         ( this%solve_mask_time(2:comp_cells_x-1,2:comp_cells_y-1) .LE. t ) )       &
         this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE.
    !$OMP END WORKSHARE
    
    !$OMP BARRIER

    IF ( bottom_radial_source_flag ) THEN

       !$OMP WORKSHARE
       WHERE ( cell_source_fractions .GT. 0.0_wp ) this%solve_mask = .TRUE.
       !$OMP END WORKSHARE
       
       !$OMP BARRIER

    END IF

    IF ( radial_source_flag ) THEN
             
       !$OMP DO private(j,k)
    
       DO k = 2,comp_cells_y-1
   
          DO j = 2,comp_cells_x-1

             IF ( source_cell(j,k) .EQ. 2 ) THEN
                
                this%solve_mask(j,k) = .TRUE.
                
             END IF

          END DO
          
       END DO

       !$OMP END DO

    END IF

    !$OMP BARRIER
    !$OMP MASTER

    DO i = 1,n_RK

       this%solve_mask_temp = this%solve_mask
       
       ! solution domain is extended to neighbours of positive-mass cells
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            this%solve_mask_temp(1:comp_cells_x-2,2:comp_cells_y-1)
       
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            this%solve_mask_temp(3:comp_cells_x,2:comp_cells_y-1)
       
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            this%solve_mask_temp(2:comp_cells_x-1,1:comp_cells_y-2)
       
       this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            this%solve_mask_temp(2:comp_cells_x-1,3:comp_cells_y)
       
    END DO

    !$OMP END MASTER
    !$OMP BARRIER

    !$OMP DO private(j,k)
    
    DO k = 1,comp_cells_y
       
       DO j = 1,comp_cells_x
          
          IF ( radial_source_flag ) THEN
             
             IF ( source_cell(j,k) .EQ. 1 ) this%solve_mask(j,k) = .FALSE.
             
          END IF
          
       END DO
       
    END DO
    
    !$OMP END DO
    !$OMP WORKSHARE

    this%solve_mask_x(1:comp_interfaces_x,1:comp_cells_y) = .FALSE.
    this%solve_mask_y(1:comp_cells_x,1:comp_interfaces_y) = .FALSE.

    !$OMP END WORKSHARE

    !$OMP END PARALLEL

    !----- check for cells where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          IF ( this%solve_mask(j,k) ) THEN

             i = i+1
             this%j_cent(i) = j
             this%k_cent(i) = k

             this%solve_mask_x(j,k) = .TRUE.
             this%solve_mask_x(j+1,k) = .TRUE.
             this%solve_mask_y(j,k) = .TRUE.
             this%solve_mask_y(j,k+1) = .TRUE.

          END IF

       END DO

    END DO

    this%solve_cells = i
    
    !----- check for y-interfaces where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_interfaces_x

          IF ( this%solve_mask_x(j,k) ) THEN

             i = i+1
             this%j_stag_x(i) = j
             this%k_stag_x(i) = k

          END IF

       END DO

    END DO

    this%solve_interfaces_x = i
   
    !----- check for y-interfaces where computation is needed

    i = 0
        
    DO k = 1,comp_interfaces_y

       DO j = 1,comp_cells_x

          IF ( this%solve_mask_y(j,k) ) THEN

             i = i+1
             this%j_stag_y(i) = j
             this%k_stag_y(i) = k

          END IF

       END DO

    END DO

    this%solve_interfaces_y = i

    RETURN

  END SUBROUTINE check_solve

END MODULE domain_2d
