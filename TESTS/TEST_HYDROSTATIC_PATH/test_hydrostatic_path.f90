PROGRAM test_hydrostatic_path

  USE parameters_2d, ONLY : wp, limiter, reconstr_coeff, theta
  USE geometry_2d, ONLY : B_cent, B_faceW, B_faceE, B_faceS, B_faceN
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, dx, dy, dx2, dy2
  USE geometry_2d, ONLY : reconstruct_topography_faces
  USE well_balanced_2d, ONLY : eval_hydrostatic_path_integral

  IMPLICIT NONE

  REAL(wp), PARAMETER :: H0 = 20.0_wp
  REAL(wp), PARAMETER :: gamma0 = 9.0E3_wp

  REAL(wp) :: energy_path
  REAL(wp) :: max_cell_error
  REAL(wp) :: momentum_path
  REAL(wp) :: pressureL, pressureR
  REAL(wp) :: reverse_energy_path
  REAL(wp) :: reverse_momentum_path
  REAL(wp) :: tolerance
  REAL(wp) :: gravW, gravE, gravS, gravN
  REAL(wp) :: hW, hE, hS, hN

  INTEGER :: j, k

  comp_cells_x = 8
  comp_cells_y = 7
  dx = 0.75_wp
  dy = 1.25_wp
  dx2 = 0.5_wp * dx
  dy2 = 0.5_wp * dy

  limiter(1) = 3
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

  CALL reconstruct_topography_faces

  max_cell_error = 0.0_wp

  DO k = 1, comp_cells_y
     DO j = 1, comp_cells_x

        gravW = 0.75_wp + 0.01_wp * REAL(j-1,wp) + 0.004_wp * REAL(k,wp)
        gravE = 0.75_wp + 0.01_wp * REAL(j,wp) + 0.004_wp * REAL(k,wp)
        gravS = 0.80_wp + 0.006_wp * REAL(j,wp) + 0.008_wp * REAL(k-1,wp)
        gravN = 0.80_wp + 0.006_wp * REAL(j,wp) + 0.008_wp * REAL(k,wp)

        hW = H0 - B_faceW(j,k)
        hE = H0 - B_faceE(j,k)
        hS = H0 - B_faceS(j,k)
        hN = H0 - B_faceN(j,k)

        CALL eval_hydrostatic_path_integral( hW, hE, gamma0, gamma0,          &
             gravW, gravE, B_faceW(j,k), B_faceE(j,k), 0.0_wp, 0.0_wp,       &
             momentum_path, energy_path )

        pressureL = 0.5_wp * gamma0 * gravW * hW**2
        pressureR = 0.5_wp * gamma0 * gravE * hE**2
        max_cell_error = MAX( max_cell_error,                                &
             ABS(momentum_path - (pressureR - pressureL)), ABS(energy_path) )

        CALL eval_hydrostatic_path_integral( hS, hN, gamma0, gamma0,          &
             gravS, gravN, B_faceS(j,k), B_faceN(j,k), 0.0_wp, 0.0_wp,       &
             momentum_path, energy_path )

        pressureL = 0.5_wp * gamma0 * gravS * hS**2
        pressureR = 0.5_wp * gamma0 * gravN * hN**2
        max_cell_error = MAX( max_cell_error,                                &
             ABS(momentum_path - (pressureR - pressureL)), ABS(energy_path) )

     END DO
  END DO

  tolerance = 4096.0_wp * EPSILON(1.0_wp) * gamma0 * H0**2

  IF ( max_cell_error .GT. tolerance ) THEN
     WRITE(*,*) 'FAIL: lake-at-rest cell path error = ', max_cell_error
     ERROR STOP 1
  END IF

  CALL eval_hydrostatic_path_integral( 1.2_wp, 0.7_wp, 8.5E3_wp, 9.1E3_wp,  &
       0.81_wp, 0.93_wp, 2.4_wp, 2.9_wp, -0.3_wp, 1.1_wp,                  &
       momentum_path, energy_path )
  CALL eval_hydrostatic_path_integral( 0.7_wp, 1.2_wp, 9.1E3_wp, 8.5E3_wp,  &
       0.93_wp, 0.81_wp, 2.9_wp, 2.4_wp, 1.1_wp, -0.3_wp,                  &
       reverse_momentum_path, reverse_energy_path )

  tolerance = 4096.0_wp * EPSILON(1.0_wp)                                  &
       * MAX(1.0_wp, ABS(momentum_path), ABS(energy_path))

  IF ( ABS(momentum_path + reverse_momentum_path) .GT. tolerance .OR.       &
       ABS(energy_path + reverse_energy_path) .GT. tolerance ) THEN
     WRITE(*,*) 'FAIL: reversed path is not antisymmetric'
     ERROR STOP 1
  END IF

  CALL eval_hydrostatic_path_integral( 1.2_wp, 0.7_wp, 8.5E3_wp, 9.1E3_wp,  &
       0.9_wp, 0.9_wp, 2.4_wp, 2.4_wp, -0.3_wp, 1.1_wp,                    &
       momentum_path, energy_path )

  IF ( momentum_path .NE. 0.0_wp .OR. energy_path .NE. 0.0_wp ) THEN
     WRITE(*,*) 'FAIL: constant geometry path is not zero'
     ERROR STOP 1
  END IF

  WRITE(*,*) 'PASS: hydrostatic path identities verified'

END PROGRAM test_hydrostatic_path
