PROGRAM test_workspace_lifecycle

  USE parameters_2d, ONLY : wp, n_vars, n_eqns, verbose_level,             &
       radial_source_flag, lateral_source_flag, bottom_radial_source_flag, &
       bottom_fissural_source_flag, n_fissures, slope_correction_flag,      &
       curvature_term_flag, liquid_vaporization_flag, lava_flag, n_solid
  USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, cell_size, x0, y0,   &
       topography_profile, nodata_topo, init_grid, B_vertex,                &
       sourceE_vect_x, sourceW_vect_x, sourceS_vect_x, sourceN_vect_x, &
       B_cent, B_face_x, B_face_y, B_prime_x_geom, grav_coeff, &
       begin_geometry_transaction, rollback_geometry_transaction, &
       commit_geometry_transaction
  USE init_2d, ONLY : thickness_init, erodible_init, release_initialization_fields
  USE state_fields, ONLY : state_type, substrate_snapshot_type
  USE substrate_parameters_2d, ONLY : substrate
  USE reconstruction_2d, ONLY : reconstruction_workspace_type
  USE hyperbolic_2d, ONLY : hyperbolic_workspace_type
  IMPLICIT NONE

  TYPE(state_type) :: state
  TYPE(substrate_snapshot_type) :: substrate_before
  TYPE(reconstruction_workspace_type) :: reconstruction
  TYPE(hyperbolic_workspace_type) :: hyperbolic
  INTEGER :: j, k
  REAL(wp), ALLOCATABLE :: bed_before(:,:), center_before(:,:), &
       face_before(:,:), slope_before(:,:), gravity_before(:,:)

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

  ! Verify that the authoritative physical inventories and their aggregate
  ! maps are restored exactly after a simulated rejected transaction.
  ALLOCATE(state%deposit(comp_cells_x,comp_cells_y,1))
  ALLOCATE(state%erosion(comp_cells_x,comp_cells_y,1))
  ALLOCATE(state%erodible(1,comp_cells_x,comp_cells_y))
  ALLOCATE(state%deposit_tot(comp_cells_x,comp_cells_y))
  ALLOCATE(state%erosion_tot(comp_cells_x,comp_cells_y))
  state%deposit=1.0_wp
  state%erosion=2.0_wp
  state%erodible=3.0_wp
  state%deposit_tot=5.0_wp
  state%erosion_tot=6.0_wp
  CALL state%snapshot_substrate(substrate_before)
  state%deposit=-11.0_wp
  state%erosion=-12.0_wp
  state%erodible=-13.0_wp
  state%deposit_tot=-15.0_wp
  state%erosion_tot=-16.0_wp
  CALL state%restore_substrate(substrate_before)
  CALL assert_true('rollback restores deposition', ALL(state%deposit==1.0_wp))
  CALL assert_true('rollback restores erosion', ALL(state%erosion==2.0_wp))
  CALL assert_true('rollback restores erodible inventory', ALL(state%erodible==3.0_wp))
  CALL assert_true('rollback restores deposit totals', ALL(state%deposit_tot==5.0_wp))
  CALL assert_true('rollback restores erosion totals', ALL(state%erosion_tot==6.0_wp))

  ! Canonical lava deposit must survive a refresh even with no solid classes,
  ! without consulting the (deliberately invalid) granular porosity.
  lava_flag=.TRUE.
  n_solid=0
  substrate%erodible_porosity=1.0_wp
  state%deposit_tot=2.5_wp
  CALL state%refresh_substrate_totals()
  CALL assert_true('lava deposit survives n_solid=0 refresh', ALL(state%deposit_tot==2.5_wp))
  CALL assert_true('lava has no granular erosion total', ALL(state%erosion_tot==0.0_wp))
  ! The canonical lava inventory must also be restored after a rejected
  ! timestep, independently of the granular class inventory and its totals.
  ! Reuse the same snapshot to check allocation-safe repeated transactions.
  CALL state%snapshot_substrate(substrate_before)
  state%deposit_tot(2,3)=9.75_wp
  state%deposit(2,3,1)=17.0_wp
  state%erosion_tot(2,3)=8.0_wp
  CALL state%restore_substrate(substrate_before)
  CALL assert_true('lava rollback restores canonical thickness', &
       ALL(state%deposit_tot==2.5_wp))
  CALL assert_true('lava rollback restores class inventory', &
       ALL(state%deposit==1.0_wp))
  CALL assert_true('lava rollback restores erosion total', &
       ALL(state%erosion_tot==0.0_wp))
  CALL state%refresh_substrate_totals()
  CALL assert_true('lava rollback survives subsequent refresh', &
       ALL(state%deposit_tot==2.5_wp))
  ! Switching back to a zero-class granular carrier must clear the
  ! aggregate instead of interpreting the lava thickness as sediment.
  lava_flag=.FALSE.
  n_solid=0
  CALL state%refresh_substrate_totals()
  CALL assert_true('granular zero-class refresh clears thickness', &
       ALL(state%deposit_tot==0.0_wp))
  CALL assert_true('granular zero-class refresh clears erosion', &
       ALL(state%erosion_tot==0.0_wp))
  lava_flag=.FALSE.
  n_solid=1
  substrate%erodible_porosity=0.25_wp
  CALL state%refresh_substrate_totals()
  CALL assert_true('granular deposit uses porosity', ALL(ABS(state%deposit_tot-4.0_wp/3.0_wp)<1.0E-12_wp))

  ! Exercise the production geometry transaction: all fields must recover
  ! exactly, including cached centre/face elevation and derived coefficients.
  bed_before=B_vertex
  center_before=B_cent
  face_before=B_face_x
  slope_before=B_prime_x_geom
  gravity_before=grav_coeff
  CALL begin_geometry_transaction
  B_vertex=B_vertex+1.0_wp
  B_cent=B_cent+2.0_wp
  B_face_x=B_face_x+3.0_wp
  B_prime_x_geom=B_prime_x_geom+4.0_wp
  grav_coeff=grav_coeff+5.0_wp
  CALL rollback_geometry_transaction
  CALL assert_true('rollback restores nodal bed exactly', ALL(B_vertex==bed_before))
  CALL assert_true('rollback restores cell geometry exactly', ALL(B_cent==center_before))
  CALL assert_true('rollback restores face geometry exactly', ALL(B_face_x==face_before))
  CALL assert_true('rollback restores slope cache exactly', ALL(B_prime_x_geom==slope_before))
  CALL assert_true('rollback restores gravity cache exactly', ALL(grav_coeff==gravity_before))
  CALL begin_geometry_transaction
  CALL commit_geometry_transaction
  CALL assert_true('commit leaves nodal bed unchanged', ALL(B_vertex==bed_before))

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
     CALL assert_true('HP scratch released', .NOT.ALLOCATED(reconstruction%hp_eta_mask))
  END DO

  PRINT *, 'PASS: workspace lifecycle and geometry transaction rollback verified'

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
