PROGRAM test_topography_faces

  USE parameters_2d, ONLY : wp, limiter, reconstr_coeff, theta
  USE geometry_2d, ONLY : B_cent, B_faceW, B_faceE, B_faceS, B_faceN
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, dx, dy, dx2, dy2
  USE geometry_2d, ONLY : limit, reconstruct_topography_faces

  IMPLICIT NONE

  REAL(wp), PARAMETER :: H0 = 20.0_wp

  REAL(wp) :: max_equilibrium_error
  REAL(wp) :: max_center_error
  REAL(wp) :: tolerance

  INTEGER :: limiter_id
  INTEGER :: j, k

  comp_cells_x = 8
  comp_cells_y = 7
  dx = 0.75_wp
  dy = 1.25_wp
  dx2 = 0.5_wp * dx
  dy2 = 0.5_wp * dy

  reconstr_coeff = 1.0_wp
  theta = 1.3_wp

  ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceW(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceE(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceS(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceN(comp_cells_x,comp_cells_y) )

  DO k = 1, comp_cells_y
     DO j = 1, comp_cells_x
        B_cent(j,k) = 2.0_wp + 0.17_wp * REAL(j,wp)                           &
             - 0.11_wp * REAL(k,wp) + 0.013_wp * REAL(j*j,wp)                &
             + 0.007_wp * REAL(j*k,wp) - 0.009_wp * REAL(k*k,wp)
     END DO
  END DO

  tolerance = 256.0_wp * EPSILON(1.0_wp) * H0

  DO limiter_id = 0, 7

     limiter(1) = limiter_id
     CALL reconstruct_topography_faces
     CALL evaluate_errors( max_equilibrium_error, max_center_error )

     IF ( max_equilibrium_error .GT. tolerance ) THEN
        WRITE(*,*) 'FAIL: limiter, equilibrium error = ',                    &
             limiter_id, max_equilibrium_error
        ERROR STOP 1
     END IF

     IF ( max_center_error .GT. tolerance ) THEN
        WRITE(*,*) 'FAIL: limiter, center error = ',                         &
             limiter_id, max_center_error
        ERROR STOP 1
     END IF

  END DO

  WRITE(*,*) 'PASS: topography face reconstruction preserves h+B'

CONTAINS

  SUBROUTINE evaluate_errors( equilibrium_error, center_error )

    REAL(wp), INTENT(OUT) :: equilibrium_error
    REAL(wp), INTENT(OUT) :: center_error

    REAL(wp) :: h_stencil(3)
    REAL(wp) :: coord_stencil(3)
    REAL(wp) :: hW, hE, hS, hN
    REAL(wp) :: slope

    equilibrium_error = 0.0_wp
    center_error = 0.0_wp

    DO k = 1, comp_cells_y
       DO j = 1, comp_cells_x

          CALL x_stencil_at_cell( j, k, h_stencil )
          coord_stencil = [ -dx, 0.0_wp, dx ]
          CALL limit( h_stencil, coord_stencil, limiter(1), slope )
          hW = H0 - B_cent(j,k) - reconstr_coeff * dx2 * slope
          hE = H0 - B_cent(j,k) + reconstr_coeff * dx2 * slope

          CALL y_stencil_at_cell( j, k, h_stencil )
          coord_stencil = [ -dy, 0.0_wp, dy ]
          CALL limit( h_stencil, coord_stencil, limiter(1), slope )
          hS = H0 - B_cent(j,k) - reconstr_coeff * dy2 * slope
          hN = H0 - B_cent(j,k) + reconstr_coeff * dy2 * slope

          equilibrium_error = MAX( equilibrium_error,                       &
               ABS(hW + B_faceW(j,k) - H0),                                 &
               ABS(hE + B_faceE(j,k) - H0),                                 &
               ABS(hS + B_faceS(j,k) - H0),                                 &
               ABS(hN + B_faceN(j,k) - H0) )

          center_error = MAX( center_error,                                 &
               ABS(0.5_wp * (B_faceW(j,k) + B_faceE(j,k)) - B_cent(j,k)),   &
               ABS(0.5_wp * (B_faceS(j,k) + B_faceN(j,k)) - B_cent(j,k)) )

       END DO
    END DO

  END SUBROUTINE evaluate_errors

  SUBROUTINE x_stencil_at_cell( j_cell, k_cell, values )

    INTEGER, INTENT(IN) :: j_cell, k_cell
    REAL(wp), INTENT(OUT) :: values(3)

    IF ( j_cell .EQ. 1 ) THEN
       values(1) = H0 - (2.0_wp * B_cent(1,k_cell) - B_cent(2,k_cell))
       values(2:3) = H0 - B_cent(1:2,k_cell)
    ELSEIF ( j_cell .EQ. comp_cells_x ) THEN
       values(1:2) = H0 - B_cent(comp_cells_x-1:comp_cells_x,k_cell)
       values(3) = H0 - (2.0_wp * B_cent(comp_cells_x,k_cell)                &
            - B_cent(comp_cells_x-1,k_cell))
    ELSE
       values = H0 - B_cent(j_cell-1:j_cell+1,k_cell)
    END IF

  END SUBROUTINE x_stencil_at_cell

  SUBROUTINE y_stencil_at_cell( j_cell, k_cell, values )

    INTEGER, INTENT(IN) :: j_cell, k_cell
    REAL(wp), INTENT(OUT) :: values(3)

    IF ( k_cell .EQ. 1 ) THEN
       values(1) = H0 - (2.0_wp * B_cent(j_cell,1) - B_cent(j_cell,2))
       values(2:3) = H0 - B_cent(j_cell,1:2)
    ELSEIF ( k_cell .EQ. comp_cells_y ) THEN
       values(1:2) = H0 - B_cent(j_cell,comp_cells_y-1:comp_cells_y)
       values(3) = H0 - (2.0_wp * B_cent(j_cell,comp_cells_y)                &
            - B_cent(j_cell,comp_cells_y-1))
    ELSE
       values = H0 - B_cent(j_cell,k_cell-1:k_cell+1)
    END IF

  END SUBROUTINE y_stencil_at_cell

END PROGRAM test_topography_faces
