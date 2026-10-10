!> \brief Capture refreshed Q1 bed geometry on arbitrary test grid sizes.
PROGRAM test_geometry_1d
  USE parameters_2d, ONLY : wp, slope_correction_flag
  USE geometry_2d
  USE omp_lib
  IMPLICIT NONE
  CHARACTER(LEN=32) :: argument
  INTEGER :: threads, actual_threads, bed, slope, j, k, output_unit
  REAL(wp) :: x, y

  CALL get_command_argument(1,argument); READ(argument,*) threads
  CALL get_command_argument(2,argument); READ(argument,*) comp_cells_x
  CALL get_command_argument(3,argument); READ(argument,*) comp_cells_y
  CALL get_command_argument(4,argument); READ(argument,*) bed
  CALL get_command_argument(5,argument); READ(argument,*) slope
  !$OMP PARALLEL
  !$OMP SINGLE
  actual_threads=omp_get_num_threads()
  !$OMP END SINGLE
  !$OMP END PARALLEL
  IF (actual_threads/=threads) ERROR STOP 'incorrect actual OpenMP team'
  comp_interfaces_x=comp_cells_x+1; comp_interfaces_y=comp_cells_y+1
  dx=0.5_wp; dy=0.75_wp; dx2=0.5_wp*dx; dy2=0.5_wp*dy
  slope_correction_flag=slope/=0
  ALLOCATE(B_vertex(comp_interfaces_x,comp_interfaces_y))
  ALLOCATE(B_face_x(comp_interfaces_x,comp_cells_y),B_face_y(comp_cells_x,comp_interfaces_y))
  ALLOCATE(B_cent(comp_cells_x,comp_cells_y))
  ALLOCATE(B_prime_x_geom(comp_cells_x,comp_cells_y),B_prime_y_geom(comp_cells_x,comp_cells_y))
  ALLOCATE(B_second_xx_geom(comp_cells_x,comp_cells_y),B_second_yy_geom(comp_cells_x,comp_cells_y))
  ALLOCATE(B_second_xy_geom(comp_cells_x,comp_cells_y),grav_coeff(comp_cells_x,comp_cells_y))
  ALLOCATE(grav_coeff_stag_x(comp_interfaces_x,comp_cells_y))
  ALLOCATE(grav_coeff_stag_y(comp_cells_x,comp_interfaces_y))
  DO k=1,comp_interfaces_y
     y=REAL(k-1,wp)*dy
     IF (comp_cells_y==1) y=0.0_wp
     DO j=1,comp_interfaces_x
        x=REAL(j-1,wp)*dx
        IF (comp_cells_x==1) x=0.0_wp
        B_vertex(j,k)=bed_value(x,y,bed)
     END DO
  END DO
  ! Poison every cache: a refresh must set all values, including inactive axes.
  B_cent=-999.0_wp; B_face_x=-999.0_wp; B_face_y=-999.0_wp
  B_prime_x_geom=-999.0_wp; B_prime_y_geom=-999.0_wp
  B_second_xx_geom=-999.0_wp; B_second_yy_geom=-999.0_wp; B_second_xy_geom=-999.0_wp
  grav_coeff=-999.0_wp; grav_coeff_stag_x=-999.0_wp; grav_coeff_stag_y=-999.0_wp
  CALL refresh_topography_geometry
  OPEN(NEWUNIT=output_unit,FILE='geometry.bin',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='REPLACE')
  WRITE(output_unit) B_vertex,B_face_x,B_face_y,B_cent
  WRITE(output_unit) B_prime_x_geom,B_prime_y_geom,B_second_xx_geom,B_second_yy_geom,B_second_xy_geom
  WRITE(output_unit) grav_coeff,grav_coeff_stag_x,grav_coeff_stag_y
  CLOSE(output_unit)
  WRITE(*,*) 'PASS: captured actual OpenMP team',actual_threads
CONTAINS
  !> \brief Evaluate a fixed flat, linear, quadratic or rough nodal test bed.
  !> \param[in] x Active x coordinate [m], zero along an inactive axis.
  !> \param[in] y Active y coordinate [m], zero along an inactive axis.
  !> \param[in] kind Bed selector from zero through three.
  !> \return Nodal elevation [m].
  PURE FUNCTION bed_value(x,y,kind) RESULT(value)
    REAL(wp), INTENT(IN) :: x,y
    INTEGER, INTENT(IN) :: kind
    REAL(wp) :: value
    value=10.0_wp
    IF (kind>=1) value=value+0.15_wp*x-0.11_wp*y
    IF (kind>=2) value=value+0.02_wp*x*x-0.01_wp*y*y+0.013_wp*x*y
    IF (kind==3) value=value+0.07_wp*SIN(0.7_wp*x+0.3_wp*y)
  END FUNCTION bed_value
END PROGRAM test_geometry_1d
