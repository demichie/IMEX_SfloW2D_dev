PROGRAM test_workspace_lifecycle

  USE parameters_2d, ONLY : wp, n_vars, n_eqns, verbose_level,             &
       radial_source_flag, lateral_source_flag, bottom_radial_source_flag, &
       bottom_fissural_source_flag, n_fissures, slope_correction_flag,      &
       curvature_term_flag, liquid_vaporization_flag
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, cell_size, x0, y0,   &
       topography_profile, nodata_topo, init_grid, B_vertex,                &
       sourceE_vect_x, sourceW_vect_x, sourceS_vect_x, sourceN_vect_x
  USE init_2d, ONLY : thickness_init, erodible_init, release_initialization_fields
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type
  USE mass_exchange_2d, ONLY : release_topography_workspace

  IMPLICIT NONE

  TYPE(reconstruction_workspace_type) :: reconstruction
  TYPE(hyperbolic_workspace_type) :: hyperbolic
  INTEGER :: j, k

  n_vars = 5
  n_eqns = 5
  verbose_level = -1
  comp_cells_x = 5
  comp_cells_y = 4
  cell_size = 1.0_wp
  x0 = 0.0_wp
  y0 = 0.0_wp
  radial_source_flag = .FALSE.
  lateral_source_flag = .FALSE.
  bottom_radial_source_flag = .FALSE.
  bottom_fissural_source_flag = .FALSE.
  slope_correction_flag = .FALSE.
  curvature_term_flag = .FALSE.
  liquid_vaporization_flag = .FALSE.
  n_fissures = 0
  nodata_topo = -9999.0_wp

  ALLOCATE(topography_profile(3,9,9))
  DO k = 1, 9
     DO j = 1, 9
        topography_profile(:,j,k) = [REAL(j-1,wp), REAL(k-1,wp), 2.0_wp]
     END DO
  END DO
  CALL init_grid
  CALL assert_true('input DEM released', .NOT.ALLOCATED(topography_profile))
  CALL assert_true('nodal bed retained', ALL(B_vertex == 2.0_wp))
  CALL assert_true('unused inlet vectors absent',                         &
       .NOT.ALLOCATED(sourceE_vect_x) .AND. .NOT.ALLOCATED(sourceW_vect_x) &
       .AND. .NOT.ALLOCATED(sourceS_vect_x) .AND. .NOT.ALLOCATED(sourceN_vect_x))

  ALLOCATE(thickness_init(comp_cells_x,comp_cells_y))
  ALLOCATE(erodible_init(comp_cells_x,comp_cells_y))
  CALL release_initialization_fields
  CALL release_initialization_fields
  CALL assert_true('initial rasters released',                           &
       .NOT.ALLOCATED(thickness_init) .AND. .NOT.ALLOCATED(erodible_init))
  CALL assert_true('bed unaffected by raster release', ALL(B_vertex == 2.0_wp))

  DO j = 1, 2
     CALL reconstruction%initialize
     CALL hyperbolic%initialize
     CALL assert_true('one speed bound per x face',                       &
          SIZE(hyperbolic%a_interface_xNeg) == (comp_cells_x+1)*comp_cells_y)
     CALL assert_true('one speed bound per y face',                       &
          SIZE(hyperbolic%a_interface_yPos) == comp_cells_x*(comp_cells_y+1))
     CALL hyperbolic%finalize
     CALL reconstruction%finalize
     CALL assert_true('speed bounds released', .NOT.ALLOCATED(hyperbolic%a_interface_xNeg))
     CALL assert_true('HP scratch released', .NOT.ALLOCATED(reconstruction%hp_scratch))
  END DO
  CALL release_topography_workspace
  CALL release_topography_workspace

  PRINT *, 'PASS: initialization storage and HP-PCCU workspace lifecycle verified'

CONTAINS

  SUBROUTINE assert_true(label, condition)
    CHARACTER(LEN=*), INTENT(IN) :: label
    LOGICAL, INTENT(IN) :: condition
    IF (.NOT.condition) THEN
       PRINT *, 'FAIL: ', label
       ERROR STOP 1
    END IF
  END SUBROUTINE assert_true

END PROGRAM test_workspace_lifecycle
