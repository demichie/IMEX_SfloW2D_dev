PROGRAM test_hydrostatic_path

  USE parameters_2d, ONLY : wp
  USE geometry_2d, ONLY : B_vertex, B_face_x, B_face_y, B_cent
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
  USE geometry_2d, ONLY : derive_topography_from_vertices
  USE nonconservative_2d, ONLY : PATH_DIR_X, PATH_DIR_Y
  USE pccu_2d, ONLY : eval_hydrostatic_path

  IMPLICIT NONE

  REAL(wp), PARAMETER :: H0 = 20.0_wp
  REAL(wp), PARAMETER :: gamma0 = 9.0E3_wp

  REAL(wp) :: path(6), reverse_path(6)
  REAL(wp) :: max_cell_error, tolerance
  REAL(wp) :: x, y, h_minus, h_plus
  INTEGER :: j, k

  comp_cells_x = 8
  comp_cells_y = 7
  comp_interfaces_x = comp_cells_x + 1
  comp_interfaces_y = comp_cells_y + 1

  ALLOCATE( B_vertex(comp_interfaces_x,comp_interfaces_y) )
  ALLOCATE( B_face_x(comp_interfaces_x,comp_cells_y) )
  ALLOCATE( B_face_y(comp_cells_x,comp_interfaces_y) )
  ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )

  DO k = 1, comp_interfaces_y
     y = REAL(k-1,wp)
     DO j = 1, comp_interfaces_x
        x = REAL(j-1,wp)
        B_vertex(j,k) = 2.0_wp + 0.17_wp*x - 0.11_wp*y                    &
             + 0.013_wp*x*x + 0.007_wp*x*y - 0.009_wp*y*y
     END DO
  END DO
  CALL derive_topography_from_vertices

  max_cell_error = 0.0_wp
  DO k = 1, comp_cells_y
     DO j = 1, comp_cells_x
        h_minus = H0-B_face_x(j,k)
        h_plus = H0-B_face_x(j+1,k)
        CALL eval_hydrostatic_path(PATH_DIR_X,h_minus,gamma0,H0,            &
             h_plus,gamma0,H0,path)
        max_cell_error = MAX(max_cell_error,MAXVAL(ABS(path)))

        h_minus = H0-B_face_y(j,k)
        h_plus = H0-B_face_y(j,k+1)
        CALL eval_hydrostatic_path(PATH_DIR_Y,h_minus,gamma0,H0,            &
             h_plus,gamma0,H0,path)
        max_cell_error = MAX(max_cell_error,MAXVAL(ABS(path)))
     END DO
  END DO

  tolerance = 4096.0_wp*EPSILON(1.0_wp)*gamma0*H0**2
  CALL assert_small('lake-at-rest cell paths',max_cell_error,tolerance)

  CALL eval_hydrostatic_path(PATH_DIR_X,1.2_wp,8.5E3_wp,2.4_wp,            &
       0.7_wp,9.1E3_wp,2.9_wp,path)
  CALL eval_hydrostatic_path(PATH_DIR_X,0.7_wp,9.1E3_wp,2.9_wp,            &
       1.2_wp,8.5E3_wp,2.4_wp,reverse_path)
  tolerance = 4096.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp,MAXVAL(ABS(path)))
  CALL assert_small('path antisymmetry',MAXVAL(ABS(path+reverse_path)),tolerance)
  CALL assert_small('thermal path component',ABS(path(4)),0.0_wp)

  WRITE(*,*) 'PASS: production hydrostatic path identities verified'

CONTAINS

  SUBROUTINE assert_small(label,value,limit_value)

    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: value, limit_value

    IF (value .GT. limit_value) THEN
       WRITE(*,*) 'FAIL: ',TRIM(label),value,' > ',limit_value
       ERROR STOP 1
    END IF

  END SUBROUTINE assert_small

END PROGRAM test_hydrostatic_path
