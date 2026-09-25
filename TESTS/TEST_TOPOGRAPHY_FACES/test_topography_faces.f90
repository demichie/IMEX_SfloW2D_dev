PROGRAM test_topography_faces

  USE parameters_2d, ONLY : wp, limiter, reconstr_coeff, theta
  USE geometry_2d, ONLY : B_vertex, B_face_x, B_face_y, B_cent
  USE geometry_2d, ONLY : B_faceW, B_faceE, B_faceS, B_faceN
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : dx, dy, dx2, dy2
  USE geometry_2d, ONLY : limit, derive_topography_from_vertices
  USE geometry_2d, ONLY : reconstruct_topography_faces
  USE geometry_2d, ONLY : project_cell_field_to_vertices

  IMPLICIT NONE

  REAL(wp), PARAMETER :: H0 = 20.0_wp

  REAL(wp), ALLOCATABLE :: cell_field(:,:), vertex_field(:,:)
  REAL(wp) :: tolerance

  comp_cells_x = 8
  comp_cells_y = 7
  comp_interfaces_x = comp_cells_x + 1
  comp_interfaces_y = comp_cells_y + 1
  dx = 0.75_wp
  dy = 1.25_wp
  dx2 = 0.5_wp * dx
  dy2 = 0.5_wp * dy

  reconstr_coeff = 1.0_wp
  theta = 1.3_wp
  tolerance = 512.0_wp * EPSILON(1.0_wp) * H0

  ALLOCATE( B_vertex(comp_interfaces_x,comp_interfaces_y) )
  ALLOCATE( B_face_x(comp_interfaces_x,comp_cells_y) )
  ALLOCATE( B_face_y(comp_cells_x,comp_interfaces_y) )
  ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceW(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceE(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceS(comp_cells_x,comp_cells_y) )
  ALLOCATE( B_faceN(comp_cells_x,comp_cells_y) )
  ALLOCATE( cell_field(comp_cells_x,comp_cells_y) )
  ALLOCATE( vertex_field(comp_interfaces_x,comp_interfaces_y) )

  CALL check_analytic_bed(.FALSE.)
  CALL check_analytic_bed(.TRUE.)
  CALL check_one_cell_excavation
  CALL check_cell_to_vertex_projection
  CALL check_legacy_lake_at_rest

  WRITE(*,*) 'PASS: continuous shared Q1 topography verified'

CONTAINS

  SUBROUTINE check_analytic_bed(curved)

    LOGICAL, INTENT(IN) :: curved
    REAL(wp) :: x, y, expected, max_error
    INTEGER :: j, k

    DO k = 1, comp_interfaces_y
       y = REAL(k-1,wp) * dy
       DO j = 1, comp_interfaces_x
          x = REAL(j-1,wp) * dx
          B_vertex(j,k) = analytic_bed(x,y,curved)
       END DO
    END DO

    CALL derive_topography_from_vertices
    CALL check_q1_identities

    max_error = 0.0_wp
    DO k = 1, comp_cells_y
       y = (REAL(k,wp)-0.5_wp) * dy
       DO j = 1, comp_cells_x
          x = (REAL(j,wp)-0.5_wp) * dx
          expected = analytic_bed(x,y,curved)
          IF (curved) expected = expected + 0.25_wp *                         &
               (0.013_wp*dx**2 - 0.009_wp*dy**2)
          max_error = MAX(max_error,ABS(B_cent(j,k)-expected))
       END DO
    END DO

    CALL assert_small('analytic Q1 center values',max_error,tolerance)

  END SUBROUTINE check_analytic_bed

  PURE FUNCTION analytic_bed(x,y,curved) RESULT(value)

    REAL(wp), INTENT(IN) :: x, y
    LOGICAL, INTENT(IN) :: curved
    REAL(wp) :: value

    value = 2.0_wp + 0.17_wp*x - 0.11_wp*y
    IF (curved) value = value + 0.013_wp*x*x + 0.007_wp*x*y - 0.009_wp*y*y

  END FUNCTION analytic_bed

  SUBROUTINE check_one_cell_excavation

    REAL(wp) :: x, y
    INTEGER :: j, k

    DO k = 1, comp_interfaces_y
       y = REAL(k-1,wp) * dy
       DO j = 1, comp_interfaces_x
          x = REAL(j-1,wp) * dx
          B_vertex(j,k) = 8.0_wp + 0.08_wp*x - 0.04_wp*y
          IF ((j .GE. 4) .AND. (j .LE. 6) .AND.                              &
              (k .GE. 3) .AND. (k .LE. 5)) B_vertex(j,k) = B_vertex(j,k)     &
               - 3.0_wp
       END DO
    END DO

    CALL derive_topography_from_vertices
    CALL check_q1_identities

    IF (.NOT. (MINVAL(B_cent(4:5,3:4)) .LT. MINVAL(B_cent(1:2,1:2)))) THEN
       WRITE(*,*) 'FAIL: Q1 excavation was not represented'
       ERROR STOP 1
    END IF

  END SUBROUTINE check_one_cell_excavation

  SUBROUTINE check_cell_to_vertex_projection

    REAL(wp) :: projection_error
    INTEGER :: j, k

    DO k = 1, comp_cells_y
       DO j = 1, comp_cells_x
          cell_field(j,k) = 0.2_wp + 0.03_wp*REAL(j,wp)                       &
               + 0.05_wp*REAL(k,wp) + 0.01_wp*REAL(j*k,wp)
       END DO
    END DO

    CALL project_cell_field_to_vertices(cell_field,vertex_field)
    B_vertex = vertex_field
    CALL derive_topography_from_vertices

    projection_error = ABS(SUM(B_cent)-SUM(cell_field))
    CALL assert_small('mass-lumped cell-to-node integral',projection_error,   &
         tolerance*REAL(comp_cells_x*comp_cells_y,wp))

    CALL assert_small('cell-to-node southwest boundary',                     &
         ABS(vertex_field(1,1)-cell_field(1,1)),tolerance)
    CALL assert_small('cell-to-node interior average',                        &
         ABS(vertex_field(3,3)-0.25_wp*SUM(cell_field(2:3,2:3))),tolerance)

  END SUBROUTINE check_cell_to_vertex_projection

  SUBROUTINE check_q1_identities

    REAL(wp) :: center_error, face_error

    face_error = MAX(                                                          &
         MAXVAL(ABS(B_face_x-0.5_wp*(B_vertex(:,1:comp_cells_y)              &
              + B_vertex(:,2:comp_interfaces_y)))),                           &
         MAXVAL(ABS(B_face_y-0.5_wp*(B_vertex(1:comp_cells_x,:)              &
              + B_vertex(2:comp_interfaces_x,:)))) )

    center_error = MAX(                                                        &
         MAXVAL(ABS(B_cent-0.5_wp*(B_face_x(1:comp_cells_x,:)                &
              + B_face_x(2:comp_interfaces_x,:)))),                           &
         MAXVAL(ABS(B_cent-0.5_wp*(B_face_y(:,1:comp_cells_y)                &
              + B_face_y(:,2:comp_interfaces_y)))) )

    CALL assert_small('Q1 shared-face identities',face_error,tolerance)
    CALL assert_small('Q1 center/face identities',center_error,tolerance)

  END SUBROUTINE check_q1_identities

  SUBROUTINE check_legacy_lake_at_rest

    REAL(wp) :: equilibrium_error, center_error
    INTEGER :: limiter_id

    DO limiter_id = 0, 7
       limiter(1) = limiter_id
       CALL reconstruct_topography_faces
       CALL evaluate_legacy_errors(equilibrium_error,center_error)
       CALL assert_small('legacy lake-at-rest reconstruction',                &
            equilibrium_error,tolerance)
       CALL assert_small('legacy reconstructed center identity',              &
            center_error,tolerance)
    END DO

  END SUBROUTINE check_legacy_lake_at_rest

  SUBROUTINE evaluate_legacy_errors(equilibrium_error,center_error)

    REAL(wp), INTENT(OUT) :: equilibrium_error, center_error
    REAL(wp) :: h_stencil(3), coord_stencil(3)
    REAL(wp) :: hW, hE, hS, hN, slope
    INTEGER :: j, k

    equilibrium_error = 0.0_wp
    center_error = 0.0_wp

    DO k = 1, comp_cells_y
       DO j = 1, comp_cells_x
          CALL x_stencil_at_cell(j,k,h_stencil)
          coord_stencil = [ -dx, 0.0_wp, dx ]
          CALL limit(h_stencil,coord_stencil,limiter(1),slope)
          hW = H0-B_cent(j,k)-reconstr_coeff*dx2*slope
          hE = H0-B_cent(j,k)+reconstr_coeff*dx2*slope

          CALL y_stencil_at_cell(j,k,h_stencil)
          coord_stencil = [ -dy, 0.0_wp, dy ]
          CALL limit(h_stencil,coord_stencil,limiter(1),slope)
          hS = H0-B_cent(j,k)-reconstr_coeff*dy2*slope
          hN = H0-B_cent(j,k)+reconstr_coeff*dy2*slope

          equilibrium_error = MAX(equilibrium_error,                          &
               ABS(hW+B_faceW(j,k)-H0),ABS(hE+B_faceE(j,k)-H0),               &
               ABS(hS+B_faceS(j,k)-H0),ABS(hN+B_faceN(j,k)-H0))
          center_error = MAX(center_error,                                    &
               ABS(0.5_wp*(B_faceW(j,k)+B_faceE(j,k))-B_cent(j,k)),           &
               ABS(0.5_wp*(B_faceS(j,k)+B_faceN(j,k))-B_cent(j,k)))
       END DO
    END DO

  END SUBROUTINE evaluate_legacy_errors

  SUBROUTINE x_stencil_at_cell(j_cell,k_cell,values)

    INTEGER, INTENT(IN) :: j_cell, k_cell
    REAL(wp), INTENT(OUT) :: values(3)

    IF (j_cell .EQ. 1) THEN
       values(1) = H0-(2.0_wp*B_cent(1,k_cell)-B_cent(2,k_cell))
       values(2:3) = H0-B_cent(1:2,k_cell)
    ELSEIF (j_cell .EQ. comp_cells_x) THEN
       values(1:2) = H0-B_cent(comp_cells_x-1:comp_cells_x,k_cell)
       values(3) = H0-(2.0_wp*B_cent(comp_cells_x,k_cell)                    &
            - B_cent(comp_cells_x-1,k_cell))
    ELSE
       values = H0-B_cent(j_cell-1:j_cell+1,k_cell)
    END IF

  END SUBROUTINE x_stencil_at_cell

  SUBROUTINE y_stencil_at_cell(j_cell,k_cell,values)

    INTEGER, INTENT(IN) :: j_cell, k_cell
    REAL(wp), INTENT(OUT) :: values(3)

    IF (k_cell .EQ. 1) THEN
       values(1) = H0-(2.0_wp*B_cent(j_cell,1)-B_cent(j_cell,2))
       values(2:3) = H0-B_cent(j_cell,1:2)
    ELSEIF (k_cell .EQ. comp_cells_y) THEN
       values(1:2) = H0-B_cent(j_cell,comp_cells_y-1:comp_cells_y)
       values(3) = H0-(2.0_wp*B_cent(j_cell,comp_cells_y)                    &
            - B_cent(j_cell,comp_cells_y-1))
    ELSE
       values = H0-B_cent(j_cell,k_cell-1:k_cell+1)
    END IF

  END SUBROUTINE y_stencil_at_cell

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF (value .GT. limit_value) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_topography_faces
