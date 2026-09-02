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

  PUBLIC :: initialize_domain
  PUBLIC :: finalize_domain
  PUBLIC :: check_solve

  PUBLIC :: solve_mask_time
  PUBLIC :: solve_mask, solve_mask_temp, solve_mask_x, solve_mask_y
  PUBLIC :: solve_cells, solve_interfaces_x, solve_interfaces_y
  PUBLIC :: j_cent, k_cent
  PUBLIC :: j_stag_x, k_stag_x, j_stag_y, k_stag_y

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

  SUBROUTINE initialize_domain

    ALLOCATE( solve_mask_time(comp_cells_x,comp_cells_y) )
    ALLOCATE( solve_mask(comp_cells_x,comp_cells_y) )
    ALLOCATE( solve_mask_temp(comp_cells_x,comp_cells_y) )
    ALLOCATE( solve_mask_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( solve_mask_y(comp_cells_x,comp_interfaces_y) )

    solve_mask_time = 0.0_wp

    solve_mask(1,1:comp_cells_y) = .TRUE.
    solve_mask(comp_cells_x,1:comp_cells_y) = .TRUE.
    solve_mask(1:comp_cells_x,1) = .TRUE.
    solve_mask(1:comp_cells_x,comp_cells_y) = .TRUE.

    comp_cells_xy = comp_cells_x * comp_cells_y

    ALLOCATE( j_cent(comp_cells_xy) )
    ALLOCATE( k_cent(comp_cells_xy) )
    ALLOCATE( j_stag_x(comp_interfaces_x*comp_cells_y) )
    ALLOCATE( k_stag_x(comp_interfaces_x*comp_cells_y) )
    ALLOCATE( j_stag_y(comp_cells_x*comp_interfaces_y) )
    ALLOCATE( k_stag_y(comp_cells_x*comp_interfaces_y) )

  END SUBROUTINE initialize_domain

  SUBROUTINE finalize_domain

    DEALLOCATE( solve_mask_time )
    DEALLOCATE( solve_mask )
    DEALLOCATE( solve_mask_temp )
    DEALLOCATE( solve_mask_x )
    DEALLOCATE( solve_mask_y )

    DEALLOCATE( j_cent, k_cent )
    DEALLOCATE( j_stag_x, k_stag_x )
    DEALLOCATE( j_stag_y, k_stag_y )

  END SUBROUTINE finalize_domain

  SUBROUTINE check_solve(q, t, solve_all)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: q(:,:,:)
    REAL(wp), INTENT(IN) :: t
    LOGICAL, INTENT(IN) :: solve_all
    
    INTEGER :: i,j,k

    !$OMP PARALLEL
         
    IF ( solve_all ) THEN

       !$OMP WORKSHARE
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE. 
       !$OMP END WORKSHARE

    ELSE
       
       !$OMP WORKSHARE
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .FALSE.
       !$OMP END WORKSHARE

    END IF
    !$OMP BARRIER
    
    !$OMP WORKSHARE
    WHERE ( ( q(1,2:comp_cells_x-1,2:comp_cells_y-1) .GT. 0.0_wp ) .AND.        &
         ( solve_mask_time(2:comp_cells_x-1,2:comp_cells_y-1) .LE. t ) )       & 
         solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) = .TRUE.
    !$OMP END WORKSHARE
    
    !$OMP BARRIER

    IF ( bottom_radial_source_flag ) THEN

       !$OMP WORKSHARE
       WHERE ( cell_source_fractions .GT. 0.0_wp ) solve_mask = .TRUE.
       !$OMP END WORKSHARE
       
       !$OMP BARRIER

    END IF

    IF ( radial_source_flag ) THEN
             
       !$OMP DO private(j,k)
    
       DO k = 2,comp_cells_y-1
   
          DO j = 2,comp_cells_x-1

             IF ( source_cell(j,k) .EQ. 2 ) THEN
                
                solve_mask(j,k) = .TRUE.
                
             END IF

          END DO
          
       END DO

       !$OMP END DO

    END IF

    !$OMP BARRIER
    !$OMP MASTER

    DO i = 1,n_RK

       solve_mask_temp = solve_mask 
       
       ! solution domain is extended to neighbours of positive-mass cells
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(1:comp_cells_x-2,2:comp_cells_y-1)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(3:comp_cells_x,2:comp_cells_y-1)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(2:comp_cells_x-1,1:comp_cells_y-2)
       
       solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) =                          &
            solve_mask(2:comp_cells_x-1,2:comp_cells_y-1) .OR.                  &
            solve_mask_temp(2:comp_cells_x-1,3:comp_cells_y) 
       
    END DO

    !$OMP END MASTER
    !$OMP BARRIER

    !$OMP DO private(j,k)
    
    DO k = 1,comp_cells_y
       
       DO j = 1,comp_cells_x
          
          IF ( radial_source_flag ) THEN
             
             IF ( source_cell(j,k) .EQ. 1 ) solve_mask(j,k) = .FALSE.
             
          END IF
          
       END DO
       
    END DO
    
    !$OMP END DO
    !$OMP WORKSHARE

    solve_mask_x(1:comp_interfaces_x,1:comp_cells_y) = .FALSE.
    solve_mask_y(1:comp_cells_x,1:comp_interfaces_y) = .FALSE.

    !$OMP END WORKSHARE

    !$OMP END PARALLEL

    !----- check for cells where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          IF ( solve_mask(j,k) ) THEN

             i = i+1
             j_cent(i) = j
             k_cent(i) = k

             solve_mask_x(j,k) = .TRUE.
             solve_mask_x(j+1,k) = .TRUE.
             solve_mask_y(j,k) = .TRUE.
             solve_mask_y(j,k+1) = .TRUE.

          END IF

       END DO

    END DO

    solve_cells = i
    
    !----- check for y-interfaces where computation is needed
    i = 0
        
    DO k = 1,comp_cells_y

       DO j = 1,comp_interfaces_x

          IF ( solve_mask_x(j,k) ) THEN

             i = i+1
             j_stag_x(i) = j
             k_stag_x(i) = k

          END IF

       END DO

    END DO

    solve_interfaces_x = i
   
    !----- check for y-interfaces where computation is needed

    i = 0
        
    DO k = 1,comp_interfaces_y

       DO j = 1,comp_cells_x

          IF ( solve_mask_y(j,k) ) THEN

             i = i+1
             j_stag_y(i) = j
             k_stag_y(i) = k

          END IF

       END DO

    END DO

    solve_interfaces_y = i

    RETURN

  END SUBROUTINE check_solve

END MODULE domain_2d

