!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!>
!> Owns grid coordinates, source geometry and the authoritative nodal bed B_vertex. Shared Q1
!> elevations feed HP-PCCU; separately filtered slopes and curvatures feed geometric corrections and
!> local rheology.
!*********************************************************************

MODULE geometry_2d

  USE parameters_2d, ONLY : wp , sp, xinf
  USE parameters_2d, ONLY : verbose_level

  IMPLICIT NONE

  !> Location of the centers (x) of the control volume of the domain
  REAL(wp), ALLOCATABLE :: x_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL(wp), ALLOCATABLE :: x_stag(:)

  !> Location of the centers (y) of the control volume of the domain
  REAL(wp), ALLOCATABLE :: y_comp(:)

  !> Location of the boundaries (y) of the control volumes of the domain
  REAL(wp), ALLOCATABLE :: y_stag(:)

  !> Authoritative continuous topography at Cartesian grid vertices.
  REAL(wp), ALLOCATABLE :: B_vertex(:,:)

  !> Unique Q1 bed elevation at x-normal and y-normal Cartesian faces.
  REAL(wp), ALLOCATABLE :: B_face_x(:,:)
  REAL(wp), ALLOCATABLE :: B_face_y(:,:)

  !> Topography at cell centers, derived from B_vertex.
  REAL(wp), ALLOCATABLE :: B_cent(:,:)

  LOGICAL, ALLOCATABLE :: B_nodata(:,:)

  INTEGER, ALLOCATABLE :: B_zone(:,:)


  ! TERMS FOR SLOPE AND AND CURVATURE CORRECTIONS

  !> Topography slope (x direction) at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_prime_x_geom(:,:)

  !> Topography 2nd x-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_xx_geom(:,:)

  !> Topography slope (y direction) at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_prime_y_geom(:,:)

  !> Topography 2nd y-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_yy_geom(:,:)

  !> Topography 2nd xy-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_xy_geom(:,:)


  !> Solution in ascii grid format (ESRI)
  REAL(wp), ALLOCATABLE :: grid_output(:,:)

  !> Integer solution in ascii grid format (ESRI)
  INTEGER, ALLOCATABLE :: grid_output_int(:,:)

  !> gravity coefficient (accounting for slope) at cell centers
  REAL(wp), ALLOCATABLE :: grav_coeff(:,:)

  !> modified gravity at cell x-faces
  REAL(wp), ALLOCATABLE :: grav_coeff_stag_x(:,:)

  !> modified gravity at cell y-faces
  REAL(wp), ALLOCATABLE :: grav_coeff_stag_y(:,:)

  !> deposit for the different classes
  REAL(wp), ALLOCATABLE :: deposit(:,:,:)

  !> total deposit 
  REAL(wp), ALLOCATABLE :: deposit_tot(:,:)

  !> erosion for the different classes
  REAL(wp), ALLOCATABLE :: erosion(:,:,:)

  !> total erosion 
  REAL(wp), ALLOCATABLE :: erosion_tot(:,:)

  !> erosdible substrate for the different classes
  REAL(wp), ALLOCATABLE :: erodible(:,:,:)

  REAL(wp), ALLOCATABLE :: topography_profile(:,:,:)

  REAL(wp) :: nodata_topo

  INTEGER, ALLOCATABLE :: source_cell(:,:)
  LOGICAL, ALLOCATABLE :: sourceE(:,:)
  LOGICAL, ALLOCATABLE :: sourceW(:,:)
  LOGICAL, ALLOCATABLE :: sourceS(:,:)
  LOGICAL, ALLOCATABLE :: sourceN(:,:)


  REAL(wp), ALLOCATABLE :: sourceE_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceE_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceW_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceW_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceS_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceS_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceN_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceN_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: cell_source_fractions(:,:)

   ! Cell-coverage fractions kept separate for each fissure.
   REAL(wp), ALLOCATABLE :: cell_fissure_fractions(:,:,:)

  !> Per-cell effective source-arc length used by the conservative lateral
  !> radial-source volume injection. Zero outside active source-boundary cells.
  REAL(wp), ALLOCATABLE :: cell_arc_perim(:,:)

  !> Outward unit normal associated with the active source arc in each cell.
  REAL(wp), ALLOCATABLE :: cell_arc_n_x(:,:)
  REAL(wp), ALLOCATABLE :: cell_arc_n_y(:,:)

  REAL(wp) :: pi_g

  INTEGER :: n_topography_profile_x, n_topography_profile_y

  REAL(wp) :: dx                 !< Control volumes size
  REAL(wp) :: x0                 !< Left of the physical domain
  REAL(wp) :: dy                 !< Control volumes size
  REAL(wp) :: y0                 !< Bottom of the physical domain
  REAL(wp) :: dx2                !< Half x Control volumes size
  REAL(wp) :: dy2                !< Half y Control volumes size

  REAL(wp) :: one_by_dx
  REAL(wp) :: one_by_dy

  INTEGER :: comp_cells_x      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_interfaces_x !< Number of interfaces (comp_cells_x+1)
  INTEGER :: comp_cells_y      !< Number of control volumes y in the comp. domain
  INTEGER :: comp_interfaces_y !< Number of interfaces (comp_cells_y+1)
  REAL(wp) :: cell_size
  INTEGER :: comp_cells_xy

CONTAINS

  !******************************************************************************
  !> \brief Allocate grid geometry and interpolate the input DEM to the authoritative nodal bed.
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !>
  !> \note Reads the configured domain and input topography_profile. Writes grid/bed/source
  !>       geometry, refreshes derived elevations and frees the input DEM after interpolation.
  !******************************************************************************

  SUBROUTINE init_grid

    USE parameters_2d, ONLY: radial_source_flag, lateral_source_flag
    USE parameters_2d, ONLY: eps_sing , eps_sing4
    USE parameters_2d, ONLY : bottom_radial_source_flag
      USE parameters_2d, ONLY : bottom_fissural_source_flag, n_fissures
    USE parameters_2d, ONLY : x_source , y_source , r_source , r2_source ,      &
         angle_source
      USE parameters_2d, ONLY : x_fissures_end_points, y_fissures_end_points,     &
          width_fissures
    USE parameters_2d, ONLY : liquid_vaporization_flag

    IMPLICIT none

    INTEGER j,k      !> loop counter

    comp_interfaces_x = comp_cells_x+1
    comp_interfaces_y = comp_cells_y+1

    ALLOCATE( x_comp(comp_cells_x) )
    ALLOCATE( x_stag(comp_interfaces_x) )
    ALLOCATE( y_comp(comp_cells_y) )
    ALLOCATE( y_stag(comp_interfaces_y) )

    ALLOCATE( source_cell(comp_cells_x,comp_cells_y) )

    ! cell where are equations are solved
    source_cell(1:comp_cells_x,1:comp_cells_y) = 0


    ALLOCATE( sourceE(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceW(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceN(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceS(comp_cells_x,comp_cells_y) )


    ! These vectors are consumed only by radial/lateral inlet boundaries.
    IF (radial_source_flag .OR. lateral_source_flag) THEN
       ALLOCATE( sourceE_vect_x(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceE_vect_y(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceW_vect_x(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceW_vect_y(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceS_vect_x(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceS_vect_y(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceN_vect_x(comp_cells_x,comp_cells_y) )
       ALLOCATE( sourceN_vect_y(comp_cells_x,comp_cells_y) )
    END IF

    ALLOCATE( B_vertex(comp_interfaces_x,comp_interfaces_y) )
    ALLOCATE( B_face_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( B_face_y(comp_cells_x,comp_interfaces_y) )
    ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_nodata(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_prime_x_geom(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_y_geom(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_second_xx_geom(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_second_yy_geom(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_second_xy_geom(comp_cells_x,comp_cells_y) )


    ALLOCATE( grid_output(comp_cells_x,comp_cells_y) )
    ALLOCATE( grid_output_int(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_coeff(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_coeff_stag_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( grav_coeff_stag_y(comp_cells_x,comp_interfaces_y) )

    ALLOCATE( cell_source_fractions(comp_cells_x,comp_cells_y) )
      ALLOCATE( cell_fissure_fractions(comp_cells_x,comp_cells_y,n_fissures) )
      cell_source_fractions = 0.0_wp
      cell_fissure_fractions = 0.0_wp

    ALLOCATE( cell_arc_perim(comp_cells_x,comp_cells_y) )
    ALLOCATE( cell_arc_n_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( cell_arc_n_y(comp_cells_x,comp_cells_y) )
    cell_arc_perim(:,:) = 0.0_wp
    cell_arc_n_x(:,:)   = 0.0_wp
    cell_arc_n_y(:,:)   = 0.0_wp

    IF ( comp_cells_x .GT. 1 ) THEN

       dx = cell_size

    ELSE

       dx = 1.0_wp

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       dy = cell_size

    ELSE

       dy = 1.0_wp

    END IF


    dx2 = dx / 2.0_wp
    dy2 = dy / 2.0_wp

    one_by_dx = 1.0_wp / dx
    one_by_dy = 1.0_wp / dy


    IF ( wp .EQ. sp ) THEN

       eps_sing=MIN(MIN( dx ** 4.0_wp,dy ** 4.0_wp ),1.0E-6_wp)

    ELSE

       eps_sing=MIN(MIN( dx ** 4.0_wp,dy ** 4.0_wp ),1.0E-10_wp)

    END IF

    eps_sing4 = eps_sing**4

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'eps_sing = ',eps_sing

    DO j=1,comp_interfaces_x

       x_stag(j) = x0 + (j-1) * dx

    END DO

    DO k=1,comp_interfaces_y

       y_stag(k) = y0 + (k-1) * dy

    END DO

    DO j=1,comp_cells_x

       x_comp(j) = 0.5_wp * ( x_stag(j) + x_stag(j+1) )

    END DO

    DO k=1,comp_cells_y

       y_comp(k) = 0.5_wp * ( y_stag(k) + y_stag(k+1) )

    END DO

    DO k=1,comp_cells_y

       DO j=1,comp_cells_x

          CALL interp_2d_nodata( topography_profile(1,:,:) ,                    &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_comp(j), y_comp(k) , B_nodata(j,k) )

       END DO

    ENDDO

    topography_profile(3,:,:) = MAX(0.0_wp,topography_profile(3,:,:))

    ! Sample the input topography directly at computational vertices. All
    ! center and face elevations are derived from this continuous Q1 field.
    DO k=1,comp_interfaces_y

       DO j=1,comp_interfaces_x

          CALL interp_2d_scalar( topography_profile(1,:,:) ,                    &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_stag(j), y_stag(k) , B_vertex(j,k) )

       END DO

    ENDDO

    CALL refresh_topography_geometry

    ALLOCATE(  B_zone(comp_cells_x,comp_cells_y) )

    B_zone = 0

    IF ( liquid_vaporization_flag ) CALL topography_zones

    pi_g = 4.0_wp * ATAN(1.0_wp)

    IF ( bottom_radial_source_flag ) THEN

       CALL compute_cell_fract(x_source,y_source,r_source,r2_source,            &
            angle_source,cell_source_fractions)

    END IF

    IF ( bottom_fissural_source_flag ) THEN

       DO j = 1, n_fissures
          CALL compute_cell_fissure_fraction(                                  &
               x_fissures_end_points(:,j), y_fissures_end_points(:,j),         &
               width_fissures(j), cell_fissure_fractions(:,:,j))
          cell_source_fractions = MAX(cell_source_fractions,                   &
               cell_fissure_fractions(:,:,j))
       END DO

    END IF

    ! The input DEM is used only above. Subsequent geometry refreshes and
    ! restarts use the authoritative nodal bed B_vertex.
    DEALLOCATE(topography_profile)

  END SUBROUTINE init_grid

  !******************************************************************************
  !> \brief Derive shared Q1 face and cell elevations from B_vertex.
  !>
  !> B_vertex is the only authoritative bed elevation. The two cells adjacent
  !> to an internal Cartesian face access the same stored B_face_x/B_face_y
  !> value, so face continuity is an exact storage identity.
  !>
  !> \note Reads B_vertex and writes B_face_x, B_face_y and B_cent without smoothing the
  !>       authoritative bed.
  !******************************************************************************

  SUBROUTINE derive_topography_from_vertices

    IMPLICIT NONE

    B_face_x(:,:) = 0.5_wp * ( B_vertex(:,1:comp_cells_y)                     &
         + B_vertex(:,2:comp_interfaces_y) )
    B_face_y(:,:) = 0.5_wp * ( B_vertex(1:comp_cells_x,:)                     &
         + B_vertex(2:comp_interfaces_x,:) )

    B_cent(:,:) = 0.25_wp * ( B_vertex(1:comp_cells_x,1:comp_cells_y)         &
         + B_vertex(2:comp_interfaces_x,1:comp_cells_y)                       &
         + B_vertex(1:comp_cells_x,2:comp_interfaces_y)                       &
         + B_vertex(2:comp_interfaces_x,2:comp_interfaces_y) )

  END SUBROUTINE derive_topography_from_vertices

  !******************************************************************************
  !> \brief Average adjacent cell values onto the shared vertex grid.
  !>
  !> Each vertex receives the arithmetic mean of all adjacent physical cells.
  !> This is the uniform-grid Q1 mass-lumped projection. It preserves the
  !> domain integral when the projected nodal field is averaged back to cells.
  !>
  !> \param[in] cell_field Cell-centered scalar proposal to be projected to the nodal grid.
  !> \param[out] vertex_field Vertex values obtained from the adjacent cells, including one-sided
  !>                          boundaries.
  !******************************************************************************

  SUBROUTINE project_cell_field_to_vertices( cell_field, vertex_field )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: cell_field(comp_cells_x,comp_cells_y)
    REAL(wp), INTENT(OUT) :: vertex_field(comp_interfaces_x,comp_interfaces_y)

    INTEGER :: j, k
    INTEGER :: j_first, j_last, k_first, k_last
    INTEGER :: adjacent_cells

    ! Gather at vertices instead of scattering from cells. This is race-free,
    ! needs no temporary weights, and gives bitwise-identical results for any
    ! OpenMP thread count because each vertex is evaluated by one iteration.
    !$OMP PARALLEL DO COLLAPSE(2)                                             &
    !$OMP & PRIVATE(j_first,j_last,k_first,k_last,adjacent_cells)
    DO k = 1, comp_interfaces_y
       DO j = 1, comp_interfaces_x
          j_first = MAX(1,j-1)
          j_last = MIN(comp_cells_x,j)
          k_first = MAX(1,k-1)
          k_last = MIN(comp_cells_y,k)
          adjacent_cells = (j_last-j_first+1) * (k_last-k_first+1)
          vertex_field(j,k) = SUM(cell_field(j_first:j_last,k_first:k_last)) &
               / REAL(adjacent_cells,wp)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE project_cell_field_to_vertices

  !******************************************************************************
  !> \brief Recompute Q1 elevations and filtered geometric derivatives after a bed update.
  !>
  !> \note Reads the current B_vertex and grid spacing. Updates all derived elevations, fitted
  !>       slopes/curvatures and centre/face gravity coefficients.
  !******************************************************************************

  SUBROUTINE refresh_topography_geometry

    IMPLICIT NONE

    ! Refresh shared Q1 elevations first. The derivative filter acts on the
    ! derived bed only; it must never overwrite the authoritative B_vertex.
    CALL derive_topography_from_vertices
    CALL topography_reconstruction

  END SUBROUTINE refresh_topography_geometry

  !******************************************************************************
  !> \brief Identify the largest connected region at the configured water elevation.
  !
  !> This subroutine search for the connected zones where topography elevation
  !> corresponds to a fixed value assigned in the input file (water_level). 
  !> A value 1 to is assigned to variable B_zone in the cells belonging to the
  !> largest connected area. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2021/07/21
  !>
  !> \note Reads B_cent and water_level; writes the connected-region indicator B_zone.
  !******************************************************************************

  SUBROUTINE topography_zones

    USE parameters_2D, ONLY: water_level

    IMPLICIT NONE

    INTEGER :: i,j,k

    INTEGER :: zone_counter
    INTEGER :: equi_list(0:1000)
    INTEGER :: equi_list_new(0:1000)
    INTEGER :: zone , zone_max
    INTEGER :: zone_cells(1:1000)

    WHERE ( ABS( B_cent - water_level ) .LE. 1.0E-1_wp )

       B_zone = -1

    END WHERE

    DO i=0,1000

       equi_list(i) = i

    END DO

    zone_counter = 0

    y_loop:DO k = 1,comp_cells_y

       x_loop:DO j = 1,comp_cells_x

          IF ( B_zone(j,k) .EQ. -1 ) THEN

             IF ( ( k .GT. 1 ) .AND. ( j.GT.1) ) THEN

                B_zone(j,k) = MAX( B_zone(j-1,k) , B_zone(j,k-1) )

                IF ( ( B_zone(j-1,k) .NE. B_zone(j,k-1) ) .AND.  &
                     (MIN( B_zone(j-1,k) , B_zone(j,k-1) ) .GT. 0 ) ) THEN

                   B_zone(j,k) = MIN( B_zone(j-1,k) , B_zone(j,k-1) )

                   equi_list(MAX( B_zone(j-1,k) , B_zone(j,k-1) )) = &
                        MIN( B_zone(j-1,k) , B_zone(j,k-1) )

                END IF

             ELSEIF ( j .GT. 1 ) THEN

                B_zone(j,k) = B_zone(j-1,k)

             ELSEIF ( k .GT. 1 ) THEN

                B_zone(j,k) = B_zone(j,k-1)

             END IF

             IF ( B_zone(j,k) .LE. 0 ) THEN

                zone_counter = zone_counter + 1
                B_zone(j,k) = zone_counter

             END IF

          END IF

       END DO x_loop

    END DO y_loop

    ! Search for smallest value in the equivalence list
    DO i=zone_counter,1,-1

       DO WHILE ( equi_list(equi_list(i)) .NE. equi_list(i) )

          equi_list(i) = equi_list(equi_list(i))

       END DO

    END DO

    equi_list_new(0:zone_counter) = 0
    i = 0

    DO zone = 1, zone_counter

       IF ( equi_list(zone) == zone ) THEN

          i = i + 1
          equi_list_new(zone) = i

       END IF
    END DO

    zone_counter = i

    zone_cells = 0

    !  Replace the labels by consecutive labels.
    y_loop2:DO k = 1,comp_cells_y

       x_loop2:DO j = 1,comp_cells_x

          B_zone(j,k) = equi_list_new(equi_list(B_zone(j,k)))

          IF ( B_zone(j,k) .GT. 0 ) THEN

             zone_cells( B_zone(j,k) ) = zone_cells( B_zone(j,k) ) + 1

          END IF

       END DO x_loop2

    END DO y_loop2

    WRITE(*,*) 'Number of cells of each conneceted area:'
    WRITE(*,*) zone_cells(1:zone_counter)

    zone_max = MAXLOC(zone_cells,1)

    y_loop3:DO k = 1,comp_cells_y

       x_loop3:DO j = 1,comp_cells_x

          IF ( B_zone(j,k) .NE. zone_max ) THEN

             B_zone(j,k) = 0

          ELSE

             B_zone(j,k) = 1

          END IF

       END DO x_loop3

    END DO y_loop3

    RETURN

  END SUBROUTINE topography_zones

  !******************************************************************************
  !> \brief Compute filtered bed slopes, curvatures and large-slope gravity coefficients.
  !
  !> A five-point polynomial least-squares filter computes the first and
  !> second bed derivatives at cell centers. The two-cell boundary band uses
  !> zero-order extrapolation from the nearest filtered interior value. The
  !> resulting slopes define the large-slope gravity correction at centers
  !> and shared Cartesian faces.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/08
  !>
  !> \note Reads B_cent, grid spacing and slope_correction_flag. Writes filtered derivative arrays
  !>       and G at cell centres/shared faces; it does not change B_vertex.
  !******************************************************************************

  SUBROUTINE topography_reconstruction

    USE parameters_2d, ONLY : slope_correction_flag

    IMPLICIT NONE

    INTEGER :: j,k,kk,kj
    REAL(wp) :: weighted_sum

    ! Five-point polynomial-fit coefficients differentiate the bed without
    ! changing the Q1 elevations used by HP-PCCU. These filtered derivatives
    ! enter large-slope/curvature corrections and local rheology instead.
    ! 1D Coefficients for 1st derivative: [-2, -1, 0, 1, 2]
    REAL(wp), PARAMETER :: c1(5) = [ -2.0_wp, -1.0_wp, 0.0_wp, 1.0_wp, 2.0_wp ]
    REAL(wp) :: norm1_x
    REAL(wp) :: norm1_y

    ! 1D Coefficients for 2nd derivative: [2, -1, -2, -1, 2]
    REAL(wp), PARAMETER :: c2(5) = [ 2.0_wp, -1.0_wp, -2.0_wp, -1.0_wp, 2.0_wp ]    
    REAL(wp) :: norm2_x
    REAL(wp) :: norm2_y

    REAL(wp) :: norm_xy

    norm1_x = 10.0_wp * dx
    norm1_y = 10.0_wp * dy
    norm2_x = 7.0_wp * dx**2
    norm2_y = 7.0_wp * dy**2
    norm_xy = (10.0_wp * dx) * (10.0_wp * dy)

    IF ( comp_cells_y .EQ. 1 ) THEN

       k = 1

       !$OMP PARALLEL DO PRIVATE(j,kj)
       DO j = 3, comp_cells_x - 2

          B_prime_x_geom(j,k)   = 0.0_wp
          B_prime_y_geom(j,k)   = 0.0_wp
          B_second_xx_geom(j,k) = 0.0_wp
          B_second_yy_geom(j,k) = 0.0_wp

          DO kj = 1, 5 ! Stencil in x-direction

             B_prime_x_geom(j,k)   = B_prime_x_geom(j,k)   + c1(kj) * B_cent(j + kj - 3, k)
             B_second_xx_geom(j,k) = B_second_xx_geom(j,k) + c2(kj) * B_cent(j + kj - 3, k)

          END DO

          B_prime_x_geom(j,k)   = B_prime_x_geom(j,k)   / norm1_x
          B_prime_y_geom(j,k)   = B_prime_y_geom(j,k)   / norm1_y
          B_second_xx_geom(j,k) = B_second_xx_geom(j,k) / norm2_x
          B_second_yy_geom(j,k) = B_second_yy_geom(j,k) / norm2_y

          B_second_xy_geom(j,k) = 0.0_wp

       END DO
       !$OMP END PARALLEL DO

    END IF

    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(j,k,kj,kk,weighted_sum)
    DO k = 3, comp_cells_y - 2

       DO j = 3, comp_cells_x - 2

          ! --- Pure derivatives (xx and yy) ---
          B_prime_x_geom(j,k)   = 0.0_wp
          B_prime_y_geom(j,k)   = 0.0_wp
          B_second_xx_geom(j,k) = 0.0_wp
          B_second_yy_geom(j,k) = 0.0_wp
          DO kj = 1, 5 ! Stencil in x-direction
             B_prime_x_geom(j,k)   = B_prime_x_geom(j,k)   + c1(kj) * B_cent(j + kj - 3, k)
             B_second_xx_geom(j,k) = B_second_xx_geom(j,k) + c2(kj) * B_cent(j + kj - 3, k)

          END DO
          DO kk = 1, 5 ! Stencil in y-direction
             B_prime_y_geom(j,k)   = B_prime_y_geom(j,k)   + c1(kk) * B_cent(j, k + kk - 3)
             B_second_yy_geom(j,k) = B_second_yy_geom(j,k) + c2(kk) * B_cent(j, k + kk - 3)
          END DO

          B_prime_x_geom(j,k)   = B_prime_x_geom(j,k)   / norm1_x
          B_prime_y_geom(j,k)   = B_prime_y_geom(j,k)   / norm1_y
          B_second_xx_geom(j,k) = B_second_xx_geom(j,k) / norm2_x
          B_second_yy_geom(j,k) = B_second_yy_geom(j,k) / norm2_y

          ! --- Mixed derivative (xy) using a 2D kernel ---
          ! The kernel is the outer product of the 1D first derivative kernels.
          weighted_sum = 0.0_wp
          DO kk = 1, 5 ! Stencil in y-direction
             DO kj = 1, 5 ! Stencil in x-direction
                weighted_sum = weighted_sum + c1(kj) * c1(kk) * B_cent(j + kj - 3, k + kk - 3)
             END DO
          END DO
          B_second_xy_geom(j,k) = weighted_sum / norm_xy

       END DO
    END DO
    !$OMP END PARALLEL DO

    !=======================================================================
    !  2. HANDLE BOUNDARY CELLS
    !     The derivatives in the 2-cell-wide border cannot be computed with
    !     the full 5-point stencil. A simple and robust strategy is to
    !     extrapolate the values from the nearest valid interior cell.
    !     This is a zero-order extrapolation (copying).
    !=======================================================================

    ! --- Handle left and right boundaries (columns j=1, 2, comp_cells_x-1, comp_cells_x) ---
    !$OMP PARALLEL DO PRIVATE(k)
    DO k = 1, comp_cells_y ! Loop over all rows
       ! Left boundary
       B_prime_x_geom(1,k)   = B_prime_x_geom(3,k)
       B_prime_x_geom(2,k)   = B_prime_x_geom(3,k)
       B_prime_y_geom(1,k)   = B_prime_y_geom(3,k)
       B_prime_y_geom(2,k)   = B_prime_y_geom(3,k)
       B_second_xx_geom(1,k) = B_second_xx_geom(3,k)
       B_second_xx_geom(2,k) = B_second_xx_geom(3,k)
       B_second_yy_geom(1,k) = B_second_yy_geom(3,k)
       B_second_yy_geom(2,k) = B_second_yy_geom(3,k)
       B_second_xy_geom(1,k) = B_second_xy_geom(3,k)
       B_second_xy_geom(2,k) = B_second_xy_geom(3,k)

       ! Right boundary
       B_prime_x_geom(comp_cells_x-1,k) = B_prime_x_geom(comp_cells_x-2,k)
       B_prime_x_geom(comp_cells_x,k)   = B_prime_x_geom(comp_cells_x-2,k)
       B_prime_y_geom(comp_cells_x-1,k) = B_prime_y_geom(comp_cells_x-2,k)
       B_prime_y_geom(comp_cells_x,k)   = B_prime_y_geom(comp_cells_x-2,k)
       B_second_xx_geom(comp_cells_x-1,k) = B_second_xx_geom(comp_cells_x-2,k)
       B_second_xx_geom(comp_cells_x,k)   = B_second_xx_geom(comp_cells_x-2,k)
       B_second_yy_geom(comp_cells_x-1,k) = B_second_yy_geom(comp_cells_x-2,k)
       B_second_yy_geom(comp_cells_x,k)   = B_second_yy_geom(comp_cells_x-2,k)
       B_second_xy_geom(comp_cells_x-1,k) = B_second_xy_geom(comp_cells_x-2,k)
       B_second_xy_geom(comp_cells_x,k)   = B_second_xy_geom(comp_cells_x-2,k)
    END DO
    !$OMP END PARALLEL DO

    IF ( comp_cells_y .GT. 1 ) THEN

       ! --- Handle top and bottom boundaries (rows k=1, 2, comp_cells_y-1, comp_cells_y) ---
       !$OMP PARALLEL DO PRIVATE(j)
       DO j = 1, comp_cells_x ! Loop over all columns
          ! Top boundary
          B_prime_x_geom(j,1)   = B_prime_x_geom(j,3)
          B_prime_x_geom(j,2)   = B_prime_x_geom(j,3)
          B_prime_y_geom(j,1)   = B_prime_y_geom(j,3)
          B_prime_y_geom(j,2)   = B_prime_y_geom(j,3)
          B_second_xx_geom(j,1) = B_second_xx_geom(j,3)
          B_second_xx_geom(j,2) = B_second_xx_geom(j,3)
          B_second_yy_geom(j,1) = B_second_yy_geom(j,3)
          B_second_yy_geom(j,2) = B_second_yy_geom(j,3)
          B_second_xy_geom(j,1) = B_second_xy_geom(j,3)
          B_second_xy_geom(j,2) = B_second_xy_geom(j,3)

          ! Bottom boundary
          B_prime_x_geom(j,comp_cells_y-1) = B_prime_x_geom(j,comp_cells_y-2)
          B_prime_x_geom(j,comp_cells_y)   = B_prime_x_geom(j,comp_cells_y-2)
          B_prime_y_geom(j,comp_cells_y-1) = B_prime_y_geom(j,comp_cells_y-2)
          B_prime_y_geom(j,comp_cells_y)   = B_prime_y_geom(j,comp_cells_y-2)
          B_second_xx_geom(j,comp_cells_y-1) = B_second_xx_geom(j,comp_cells_y-2)
          B_second_xx_geom(j,comp_cells_y)   = B_second_xx_geom(j,comp_cells_y-2)
          B_second_yy_geom(j,comp_cells_y-1) = B_second_yy_geom(j,comp_cells_y-2)
          B_second_yy_geom(j,comp_cells_y)   = B_second_yy_geom(j,comp_cells_y-2)
          B_second_xy_geom(j,comp_cells_y-1) = B_second_xy_geom(j,comp_cells_y-2)
          B_second_xy_geom(j,comp_cells_y)   = B_second_xy_geom(j,comp_cells_y-2)
       END DO
       !$OMP END PARALLEL DO

    END IF

    IF ( slope_correction_flag ) THEN

       ! Calculate grav_coeff using the retained LS first derivatives.
       grav_coeff = 1.0_wp / ( 1.0_wp + B_prime_x_geom**2 + B_prime_y_geom**2 )

       ! Interpolate the cell coefficient to the shared Cartesian faces.
       grav_coeff_stag_x(1,:) = grav_coeff(1,:)
       grav_coeff_stag_x(2:comp_interfaces_x-1,:) = 0.5_wp *                    &
            ( grav_coeff(1:comp_cells_x-1,:) + grav_coeff(2:comp_cells_x,:) )
       grav_coeff_stag_x(comp_interfaces_x,:) = grav_coeff(comp_cells_x,:)

       grav_coeff_stag_y(:,1) = grav_coeff(:,1)
       grav_coeff_stag_y(:,2:comp_interfaces_y-1) = 0.5_wp *                    &
            ( grav_coeff(:,1:comp_cells_y-1) + grav_coeff(:,2:comp_cells_y) )
       grav_coeff_stag_y(:,comp_interfaces_y) = grav_coeff(:,comp_cells_y)

    ELSE

       grav_coeff = 1.0_wp

       grav_coeff_stag_x = 1.0_wp
       grav_coeff_stag_y = 1.0_wp

    END IF

    RETURN

  END SUBROUTINE topography_reconstruction

  !******************************************************************************
  !> \brief Identify inlet cells and build their radial/lateral emission directions.
  !
  !> In this subroutine the source of mass is initialized. The cells belonging
  !> to the source are are identified ( source_cell(j,k) = 2 ).
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2021/04/30
  !>
  !> \note Reads radial/lateral source parameters and coordinates. Writes source masks, inlet
  !>       direction vectors and the integrated source perimeter.
  !******************************************************************************

  SUBROUTINE init_source

    USE parameters_2d, ONLY : x_source , y_source , r_source

    USE parameters_2d, ONLY : source_side
    USE parameters_2d, ONLY : x1_source , x2_source
    USE parameters_2d, ONLY : y1_source , y2_source
    USE parameters_2d, ONLY : lateral_source_flag
    USE parameters_2d, ONLY : radial_source_flag
    USE parameters_2d, ONLY : azimuth_source , arc_width_source
    USE parameters_2d, ONLY : r2_source

    IMPLICIT NONE

    REAL(wp) :: source_distance

    INTEGER :: j,k

    REAL(wp) :: side_fract
    REAL(wp) :: total_source

    REAL(wp) :: face_len
    REAL(wp) :: vx, vy, vmag, nx, ny
    REAL(wp) :: ang_compass, ang_diff, half_arc
    LOGICAL  :: in_arc
    REAL(wp) :: arc_perim_total
    REAL(wp) :: source_length_analytical
    REAL(wp) :: h_ell, rescale

    ! cells where are equations are solved
    source_cell(1:comp_cells_x,1:comp_cells_y) = 0

    cell_arc_perim(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    cell_arc_n_x(1:comp_cells_x,1:comp_cells_y)   = 0.0_wp
    cell_arc_n_y(1:comp_cells_x,1:comp_cells_y)   = 0.0_wp

    sourceE(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceW(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceN(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceS(1:comp_cells_x,1:comp_cells_y) = .FALSE.


    total_source = 0.0_wp

    IF ( lateral_source_flag ) THEN

       IF ( source_side .EQ. 'W' ) THEN

          DO k = 2,comp_cells_y-1

             IF ( ( y_stag(k) .LE. y2_source ) .AND. ( y_stag(k+1) .GE.         &
                  y1_source ) ) THEN

                source_cell(1,k) = 2
                sourceW(1,k) = .TRUE.

                side_fract = MIN ( MIN( (y_stag(k+1) - y1_source) , ( y2_source -     &
                     y_stag(k) ) ) , ( y2_source - y1_source ) , dy ) / dy

                WRITE(*,*) 'side_fract',side_fract

                sourceW_vect_x(1,k) = side_fract
                sourceW_vect_y(1,k) = 0.0_wp

                total_source = total_source + dy * side_fract

             END IF

          END DO

       END IF

       IF ( source_side .EQ. 'E' ) THEN

          DO k = 2,comp_cells_y-1

             IF ( ( y_stag(k) .LE. y2_source ) .AND. ( y_stag(k+1) .GE.         &
                  y1_source ) ) THEN

                source_cell(comp_cells_x,k) = 2
                sourceE(comp_cells_x,k) = .TRUE.

                side_fract = MIN ( MIN( (y_stag(k+1) - y1_source) , ( y2_source -     &
                     y_stag(k) ) ) , ( y2_source - y1_source ) , dy ) / dy

                sourceE_vect_x(comp_cells_x,k) = -side_fract
                sourceE_vect_y(comp_cells_x,k) = 0.0_wp

                total_source = total_source + dy * side_fract

             END IF

          END DO

       END IF

       IF ( source_side .EQ. 'S' ) THEN

          DO j = 2,comp_cells_x-1

             IF ( ( x_stag(j) .LE. x2_source ) .AND. ( x_stag(j+1) .GE.         &
                  x1_source ) ) THEN

                source_cell(j,1) = 2
                sourceS(j,1) = .TRUE.

                side_fract = MIN ( MIN( (x_stag(j+1) - x1_source) , ( x2_source -     &
                     x_stag(j) ) ) , ( x2_source - x1_source ) , dx ) / dx

                sourceS_vect_x(j,1) = 0.0_wp
                sourceS_vect_y(j,1) = side_fract

                total_source = total_source + dy * side_fract

             END IF

          END DO

       END IF

       IF ( source_side .EQ. 'N' ) THEN

          DO j = 2,comp_cells_x-1

             IF ( ( x_stag(j) .LE. x2_source ) .AND. ( x_stag(j+1) .GE.         &
                  x1_source ) ) THEN

                source_cell(j,1) = 2
                sourceN(j,1) = .TRUE.

                side_fract = MIN ( MIN( (x_stag(j+1) - x1_source) , ( x2_source -     &
                     x_stag(j) ) ) , ( x2_source - x1_source ) , dx ) / dx

                sourceS_vect_x(j,1) = 0.0_wp
                sourceS_vect_y(j,1) = - side_fract

                total_source = total_source + dy * side_fract

             END IF

          END DO

       END IF

       WRITE(*,*) 'Total source length (m) = ',total_source

       RETURN

    END IF

    DO k = 2,comp_cells_y-1

       DO j = 2,comp_cells_x-1

          IF ( ( x_comp(j) - x_source )**2 + ( y_comp(k) - y_source )**2 .LE.   &
               r_source**2 ) THEN

             ! cells where equations are not solved
             source_cell(j,k) = 1 

             ! check on west cell
             IF ( ( x_comp(j-1) - x_source )**2 + ( y_comp(k) - y_source )**2   &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j-1,k) = 2
                sourceE(j-1,k) = .TRUE.
                source_distance = SQRT( ( x_stag(j) - x_source )**2         &
                     + ( y_comp(k) - y_source )**2 )

                sourceE_vect_x(j-1,k) = ( x_stag(j) - x_source ) * r_source     &
                     / source_distance**2

                sourceE_vect_y(j-1,k) = ( y_comp(k) - y_source ) * r_source     &
                     / source_distance**2

                total_source = total_source + dx * ABS( sourceE_vect_x(j-1,k) )

             ELSEIF ( ( x_comp(j+1) - x_source )**2 + ( y_comp(k)-y_source )**2 &
                  .GE. r_source**2 ) THEN
                ! check on east cell

                ! cells where radial source boundary condition are applied
                source_cell(j+1,k) = 2
                sourceW(j+1,k) = .TRUE.
                source_distance = SQRT( ( x_stag(j+1) - x_source )**2       &
                     + ( y_comp(k) - y_source )**2 )

                sourceW_vect_x(j+1,k) = ( x_stag(j+1) - x_source ) * r_source   &
                     / source_distance**2

                sourceW_vect_y(j+1,k) = ( y_comp(k) - y_source ) * r_source     &
                     / source_distance**2

                total_source = total_source + dx * ABS( sourceW_vect_x(j+1,k) )

             END IF

             ! check on south cell
             IF ( ( x_comp(j) - x_source )**2 + ( y_comp(k-1) - y_source )**2   &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j,k-1) = 2
                sourceN(j,k-1) = .TRUE.
                source_distance = SQRT( ( x_comp(j) - x_source )**2         &
                     + ( y_stag(k) - y_source )**2 )

                sourceN_vect_x(j,k-1) = ( x_comp(j) - x_source ) * r_source     &
                     / source_distance**2

                sourceN_vect_y(j,k-1) = ( y_stag(k) - y_source ) * r_source     &
                     / source_distance**2

                total_source = total_source + dy * ABS( sourceN_vect_y(j,k-1) )

             ELSEIF ( ( x_comp(j)-x_source )**2 + ( y_comp(k+1) - y_source )**2 &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j,k+1) = 2
                sourceS(j,k+1) = .TRUE.
                source_distance = SQRT( ( x_comp(j) - x_source )**2         &
                     + ( y_stag(k+1) - y_source )**2 )

                sourceS_vect_x(j,k+1) = ( x_comp(j) - x_source ) * r_source     &
                     / source_distance**2

                sourceS_vect_y(j,k+1) = ( y_stag(k+1) - y_source ) * r_source   &
                     / source_distance**2

                total_source = total_source + dy * ABS( sourceS_vect_y(j,k+1) )

             END IF

          END IF

       END DO

    END DO

    ! ----------------------------------------------------------------------
    ! Volume-source weights for the lateral directed RADIAL_SOURCE.
    !
    ! The lateral source is reformulated as a volume source term in
    ! eval_expl_terms to guarantee that the integrated MFR equals what the
    ! user prescribes. For each ring-boundary cell (source_cell == 2),
    ! accumulate the in-arc face-length contribution and a face-length-
    ! weighted outward unit normal. cell_arc_perim is then rescaled to match
    ! the analytical arc length so that the per-cell injection rate
    ! rho_m * h * v * cell_arc_perim / cell_area integrates to MFR exactly.
    ! ----------------------------------------------------------------------
    IF ( radial_source_flag ) THEN

       half_arc = 0.5_wp * arc_width_source

       DO k = 1,comp_cells_y
          DO j = 1,comp_cells_x

             IF ( source_cell(j,k) .NE. 2 ) CYCLE

             ! East face
             IF ( sourceE(j,k) ) THEN
                vx = sourceE_vect_x(j,k)
                vy = sourceE_vect_y(j,k)
                vmag = SQRT( vx*vx + vy*vy )
                IF ( vmag .GT. 0.0_wp ) THEN
                   nx = vx / vmag
                   ny = vy / vmag
                   IF ( arc_width_source .GE. 360.0_wp ) THEN
                      in_arc = .TRUE.
                   ELSE
                      ang_compass = 90.0_wp - ATAN2(vy,vx) * 180.0_wp / pi_g
                      ang_compass = MOD(ang_compass + 360.0_wp, 360.0_wp)
                      ang_diff = MOD(ABS(ang_compass - azimuth_source)             &
                           + 360.0_wp, 360.0_wp)
                      IF ( ang_diff .GT. 180.0_wp ) ang_diff = 360.0_wp - ang_diff
                      in_arc = ( ang_diff .LE. half_arc )
                   END IF
                   IF ( in_arc ) THEN
                      face_len = dy
                      cell_arc_perim(j,k) = cell_arc_perim(j,k) + face_len
                      cell_arc_n_x(j,k)   = cell_arc_n_x(j,k)   + face_len * nx
                      cell_arc_n_y(j,k)   = cell_arc_n_y(j,k)   + face_len * ny
                   END IF
                END IF
             END IF

             ! West face
             IF ( sourceW(j,k) ) THEN
                vx = sourceW_vect_x(j,k)
                vy = sourceW_vect_y(j,k)
                vmag = SQRT( vx*vx + vy*vy )
                IF ( vmag .GT. 0.0_wp ) THEN
                   nx = vx / vmag
                   ny = vy / vmag
                   IF ( arc_width_source .GE. 360.0_wp ) THEN
                      in_arc = .TRUE.
                   ELSE
                      ang_compass = 90.0_wp - ATAN2(vy,vx) * 180.0_wp / pi_g
                      ang_compass = MOD(ang_compass + 360.0_wp, 360.0_wp)
                      ang_diff = MOD(ABS(ang_compass - azimuth_source)             &
                           + 360.0_wp, 360.0_wp)
                      IF ( ang_diff .GT. 180.0_wp ) ang_diff = 360.0_wp - ang_diff
                      in_arc = ( ang_diff .LE. half_arc )
                   END IF
                   IF ( in_arc ) THEN
                      face_len = dy
                      cell_arc_perim(j,k) = cell_arc_perim(j,k) + face_len
                      cell_arc_n_x(j,k)   = cell_arc_n_x(j,k)   + face_len * nx
                      cell_arc_n_y(j,k)   = cell_arc_n_y(j,k)   + face_len * ny
                   END IF
                END IF
             END IF

             ! North face
             IF ( sourceN(j,k) ) THEN
                vx = sourceN_vect_x(j,k)
                vy = sourceN_vect_y(j,k)
                vmag = SQRT( vx*vx + vy*vy )
                IF ( vmag .GT. 0.0_wp ) THEN
                   nx = vx / vmag
                   ny = vy / vmag
                   IF ( arc_width_source .GE. 360.0_wp ) THEN
                      in_arc = .TRUE.
                   ELSE
                      ang_compass = 90.0_wp - ATAN2(vy,vx) * 180.0_wp / pi_g
                      ang_compass = MOD(ang_compass + 360.0_wp, 360.0_wp)
                      ang_diff = MOD(ABS(ang_compass - azimuth_source)             &
                           + 360.0_wp, 360.0_wp)
                      IF ( ang_diff .GT. 180.0_wp ) ang_diff = 360.0_wp - ang_diff
                      in_arc = ( ang_diff .LE. half_arc )
                   END IF
                   IF ( in_arc ) THEN
                      face_len = dx
                      cell_arc_perim(j,k) = cell_arc_perim(j,k) + face_len
                      cell_arc_n_x(j,k)   = cell_arc_n_x(j,k)   + face_len * nx
                      cell_arc_n_y(j,k)   = cell_arc_n_y(j,k)   + face_len * ny
                   END IF
                END IF
             END IF

             ! South face
             IF ( sourceS(j,k) ) THEN
                vx = sourceS_vect_x(j,k)
                vy = sourceS_vect_y(j,k)
                vmag = SQRT( vx*vx + vy*vy )
                IF ( vmag .GT. 0.0_wp ) THEN
                   nx = vx / vmag
                   ny = vy / vmag
                   IF ( arc_width_source .GE. 360.0_wp ) THEN
                      in_arc = .TRUE.
                   ELSE
                      ang_compass = 90.0_wp - ATAN2(vy,vx) * 180.0_wp / pi_g
                      ang_compass = MOD(ang_compass + 360.0_wp, 360.0_wp)
                      ang_diff = MOD(ABS(ang_compass - azimuth_source)             &
                           + 360.0_wp, 360.0_wp)
                      IF ( ang_diff .GT. 180.0_wp ) ang_diff = 360.0_wp - ang_diff
                      in_arc = ( ang_diff .LE. half_arc )
                   END IF
                   IF ( in_arc ) THEN
                      face_len = dx
                      cell_arc_perim(j,k) = cell_arc_perim(j,k) + face_len
                      cell_arc_n_x(j,k)   = cell_arc_n_x(j,k)   + face_len * nx
                      cell_arc_n_y(j,k)   = cell_arc_n_y(j,k)   + face_len * ny
                   END IF
                END IF
             END IF

             ! Normalize the cell-aggregate outward unit normal.
             IF ( cell_arc_perim(j,k) .GT. 0.0_wp ) THEN
                cell_arc_n_x(j,k) = cell_arc_n_x(j,k) / cell_arc_perim(j,k)
                cell_arc_n_y(j,k) = cell_arc_n_y(j,k) / cell_arc_perim(j,k)
                vmag = SQRT( cell_arc_n_x(j,k)**2 + cell_arc_n_y(j,k)**2 )
                IF ( vmag .GT. 0.0_wp ) THEN
                   cell_arc_n_x(j,k) = cell_arc_n_x(j,k) / vmag
                   cell_arc_n_y(j,k) = cell_arc_n_y(j,k) / vmag
                END IF
             END IF

          END DO
       END DO

       arc_perim_total = SUM( cell_arc_perim(1:comp_cells_x,1:comp_cells_y) )
       WRITE(*,*) 'Radial source arc perimeter (discretized, m) = ',           &
            arc_perim_total

       ! Rescale cell_arc_perim so that its sum equals the analytical source
       ! length used by inpout to derive h_source and vel_source. This makes
       ! the per-cell injection rho_m * h * v * cell_arc_perim / cell_area
       ! integrate to the user-prescribed MFR exactly.
       h_ell = (r_source - r2_source)**2 / (r_source + r2_source)**2
       source_length_analytical = pi_g * (r_source + r2_source) *              &
            ( 1.0_wp + 0.25_wp * h_ell + h_ell**2 / 64.0_wp                    &
              + h_ell**3 / 256.0_wp + h_ell**4 * 25.0_wp / 16384.0_wp          &
              + h_ell**5 * 49.0_wp / 65536.0_wp )

       IF ( arc_width_source .LT. 360.0_wp ) THEN
          source_length_analytical = source_length_analytical                  &
               * ( arc_width_source / 360.0_wp )
       END IF

       IF ( arc_perim_total .GT. 0.0_wp ) THEN
          rescale = source_length_analytical / arc_perim_total
          cell_arc_perim(:,:) = cell_arc_perim(:,:) * rescale
          WRITE(*,*) 'Radial source analytical arc length (m) = ',             &
               source_length_analytical
          WRITE(*,*) 'Volume-source rescale factor              = ', rescale
       END IF

    END IF

    RETURN

  END SUBROUTINE init_source

  !-----------------------------------------------------------------------------
  !> \brief Bilinearly interpolate a scalar on a grid supplied as coordinate matrices.
  !>
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !>
  !> \param[in] x1 X coordinates of the original interpolation grid.
  !> \param[in] y1 Y coordinates of the original interpolation grid.
  !> \param[in] f1 Scalar values on the original grid.
  !> \param[in] x2 X coordinate of the requested interpolation point [m].
  !> \param[in] y2 Y coordinate of the requested interpolation point [m].
  !> \param[out] f2 Interpolated scalar value at (x2,y2).
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalar(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL(wp), INTENT(IN) :: x2, y2
    REAL(wp), INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF

    IF ( size(x1,1) .EQ. 1 ) THEN

       f2 = alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1)

    ELSEIF ( size(x1,2) .EQ. 1 ) THEN

       f2 = alfa_x * f1(ix,iy)  + ( 1.0_wp - alfa_x ) * f1(ix+1,iy)

    ELSE

       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1) ) &
            + ( 1.0_wp - alfa_x ) * ( alfa_y * f1(ix+1,iy) + ( 1.0_wp - alfa_y )&
            * f1(ix+1,iy+1) )

    END IF

  END SUBROUTINE interp_2d_scalar

  !-----------------------------------------------------------------------------
  !> \brief Check whether a bilinear interpolation stencil includes missing DEM data.
  !>
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !>
  !> \param[in] x1 X coordinates of the original interpolation grid.
  !> \param[in] y1 Y coordinates of the original interpolation grid.
  !> \param[in] f1 Scalar values on the original grid.
  !> \param[in] x2 X coordinate of the requested interpolation point [m].
  !> \param[in] y2 Y coordinate of the requested interpolation point [m].
  !> \param[out] f2 True if any required input stencil value equals nodata_topo.
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_nodata(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL(wp), INTENT(IN) :: x2, y2
    LOGICAL, INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF

    f2 = .FALSE.

    IF ( size(x1,1) .EQ. 1 ) THEN

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix,iy+1) .EQ. nodata_topo )

    ELSEIF ( size(x1,2) .EQ. 1 ) THEN

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix+1,iy) .EQ. nodata_topo ) 

    ELSE

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix,iy+1) .EQ. nodata_topo )&
            .OR. ( f1(ix+1,iy) .EQ. nodata_topo )                               &
            .OR. ( f1(ix+1,iy+1) .EQ. nodata_topo )

    END IF

  END SUBROUTINE interp_2d_nodata


  !-----------------------------------------------------------------------------
  !> \brief Bilinearly interpolate a scalar using one-dimensional coordinate arrays.
  !>
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !>
  !> \param[in] x1 X coordinates of the original interpolation grid.
  !> \param[in] y1 Y coordinates of the original interpolation grid.
  !> \param[in] f1 Scalar values on the original grid.
  !> \param[in] x2 X coordinate of the requested interpolation point [m].
  !> \param[in] y2 Y coordinate of the requested interpolation point [m].
  !> \param[out] f2 Interpolated scalar value at (x2,y2).
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalarB(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:) :: x1, y1
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: f1
    REAL(wp), INTENT(IN) :: x2, y2
    REAL(wp), INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1) ) / ( x1(2) - x1(1) ) ) + 1
       ix = MAX(0,MIN( ix , SIZE(x1)-1 ))
       alfa_x = ( x1(ix+1) - x2 ) / (  x1(ix+1) - x1(ix) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(y1) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1) ) / ( y1(2) - y1(1) ) ) + 1
       iy = MAX(1,MIN( iy , SIZE(y1)-1 ))
       alfa_y = ( y1(iy+1) - y2 ) / (  y1(iy+1) - y1(iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF

    IF ( ( alfa_x .LT. 0.0_wp ) .OR. ( alfa_x .GT. 1.0_wp )                     &
         .OR. ( alfa_y .LT. 0.0_wp ) .OR. ( alfa_y .GT. 1.0_wp ) ) THEN

       f2 = 0.0_wp
       RETURN

    END IF


    IF ( size(x1) .EQ. 1 ) THEN

       f2 = alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1)

    ELSEIF ( size(y1) .EQ. 1 ) THEN

       f2 = alfa_x * f1(ix,iy)  + ( 1.0_wp - alfa_x ) * f1(ix+1,iy)

    ELSE

       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1) ) &
            + ( 1.0_wp - alfa_x ) * ( alfa_y * f1(ix+1,iy) + ( 1.0_wp - alfa_y )&
            * f1(ix+1,iy+1) )

    END IF

    RETURN

  END SUBROUTINE interp_2d_scalarB


  !-----------------------------------------------------------------------------
  !> \brief Average an input raster over one target control volume using overlap weights.
  !>
  !> Scalar regrid (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !>
  !> \param[in] xin X coordinates of the input raster grid [m].
  !> \param[in] yin Y coordinates of the input raster grid [m].
  !> \param[in] fin Input raster values to average over the target cell.
  !> \param[in] xl Left boundary of the target cell [m].
  !> \param[in] xr Right boundary of the target cell [m].
  !> \param[in] yl Lower boundary of the target cell [m].
  !> \param[in] yr Upper boundary of the target cell [m].
  !> \param[out] fout Overlap-weighted scalar average over the target cell.
  !-----------------------------------------------------------------------------

  SUBROUTINE regrid_scalar(xin, yin, fin, xl, xr , yl, yr, fout)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:) :: xin, yin
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: fin
    REAL(wp), INTENT(IN) :: xl, xr , yl , yr
    REAL(wp), INTENT(OUT) :: fout

    INTEGER :: ix , iy
    INTEGER :: ix1 , ix2 , iy1 , iy2
    REAL(wp) :: alfa_x , alfa_y
    REAL(wp) :: dXin , dYin

    INTEGER nXin,nYin

    nXin = size(xin)-1
    nYin = size(yin)

    dXin = xin(2) - xin(1)
    dYin = yin(2) - yin(1)

    ix1 = MAX(1,CEILING( ( xl - xin(1) ) / dXin ))
    ix2 = MIN(nXin,CEILING( ( xr -xin(1) ) / dXin )+1)

    iy1 = MAX(1,CEILING( ( yl - yin(1) ) / dYin ))
    iy2 = MIN(nYin,CEILING( ( yr - yin(1) ) / dYin ) + 1)

    fout = 0.0_wp

    DO ix=ix1,ix2-1

       alfa_x = ( MIN(xr,xin(ix+1)) - MAX(xl,xin(ix)) ) / ( xr - xl )

       DO iy=iy1,iy2-1

          alfa_y = ( MIN(yr,yin(iy+1)) - MAX(yl,yin(iy)) ) / ( yr - yl )

          fout = fout + alfa_x * alfa_y * fin(ix,iy)

       END DO

    END DO

  END SUBROUTINE regrid_scalar

  !******************************************************************************
  !> \brief Apply the configured slope limiter to a three-point stencil.
  !
  !> This subroutine limits the slope of the linear reconstruction of 
  !> the physical variables, accordingly to the parameter "solve_limiter":\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] v Scalar values at the three stencil points.
  !> \param[in] z Coordinates of the three stencil points.
  !> \param[in] limiter Limiter selector; the supported cases are defined by this routine.
  !> \param[out] slope_lim Limited derivative along the stencil coordinate.
  !******************************************************************************

  SUBROUTINE limit( v , z , limiter , slope_lim )

    USE parameters_2d, ONLY : theta

    IMPLICIT none

    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp), INTENT(IN) :: z(3)
    INTEGER, INTENT(IN) :: limiter

    REAL(wp), INTENT(OUT) :: slope_lim

    REAL(wp) :: a , b , c

    REAL(wp) :: sigma1 , sigma2

    a = ( v(3) - v(2) ) / ( z(3) - z(2) )
    b = ( v(2) - v(1) ) / ( z(2) - z(1) )
    c = ( v(3) - v(1) ) / ( z(3) - z(1) )

    SELECT CASE (limiter)

    CASE ( 0 )

       slope_lim = 0.0_wp

    CASE ( 1 )

       ! minmod
       slope_lim = minmod(a,b)

    CASE ( 2 )

       ! superbee
       sigma1 = minmod( a , 2.0_wp*b )
       sigma2 = minmod( 2.0_wp*a , b )
       slope_lim = maxmod( sigma1 , sigma2 )

    CASE ( 3 )

       ! generalized minmod
       slope_lim = minmod( c , theta * minmod( a , b ) )

    CASE ( 4 )

       ! monotonized central-difference (MC, LeVeque p.112)
       slope_lim = minmod( c , 2.0 * minmod( a , b ) )

    CASE ( 5 )

       ! centered
       slope_lim = c

    CASE (6) 

       ! backward
       slope_lim = a

    CASE (7)

       !forward
       slope_lim = b

    END SELECT

  END SUBROUTINE limit

  !******************************************************************************
  !> \brief Return the smaller same-sign slope, or zero for incompatible slopes.
  !
  !> This function compute the minmod between two real numbers
  !> \date 2021/04/30
  !> @author 
  !> Mattia de' Michieli Vitturi
  !>
  !> \param[in] a First candidate slope.
  !> \param[in] b Second candidate slope.
  !> \return Same-sign minimum-magnitude slope, or zero near zero/opposite signs.
  !******************************************************************************

  REAL(wp) FUNCTION minmod(a,b)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: b
    REAL(wp) :: sa , sb 

    IF ( MIN(ABS(a),ABS(b)) .LE. 0.0e-40_wp ) THEN

       minmod = 0.0_wp

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       minmod = 0.5_wp * ( sa+sb ) * MIN( ABS(a) , ABS(b) )

    END IF

  END FUNCTION minmod

  !> \brief Return the larger same-sign slope, or zero for incompatible slopes.
  !>
  !> \param[in] a First candidate slope.
  !> \param[in] b Second candidate slope.
  !> \return Same-sign maximum-magnitude slope, or zero near zero/opposite signs.

  REAL(wp) function maxmod(a,b)

    IMPLICIT none

    REAL(wp) :: a , b , sa , sb 

    IF ( ABS(a*b) .LE. 0.0e-30_wp ) THEN

       maxmod = 0.0_wp

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       maxmod = 0.5_wp * ( sa+sb ) * MAX( ABS(a) , ABS(b) )

    END IF

  END function maxmod

  ! Fissural rectangle/cell intersections follow the segment-and-width source
  ! geometry introduced by Elisa Biagioli in BiElisa/IMEX_LavaFlow.
  !> \brief Compute each cell overlap fraction with a finite-width fissure rectangle.
  !>
  !> \param[in] endpoints_x X coordinates of the two fissure segment endpoints [m].
  !> \param[in] endpoints_y Y coordinates of the two fissure segment endpoints [m].
  !> \param[in] width Full width of the fissure rectangle [m].
  !> \param[out] cell_fraction Rectangle overlap divided by cell area, one value per cell.

  SUBROUTINE compute_cell_fissure_fraction(endpoints_x, endpoints_y, width,   &
       cell_fraction)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: endpoints_x(2), endpoints_y(2), width
    REAL(wp), INTENT(OUT) :: cell_fraction(comp_cells_x,comp_cells_y)

    REAL(wp) :: polygon_x(16), polygon_y(16)
    REAL(wp) :: clipped_x(16), clipped_y(16)
    REAL(wp) :: segment_x, segment_y, segment_length
    REAL(wp) :: normal_x, normal_y, half_width
    REAL(wp) :: min_x, max_x, min_y, max_y
    REAL(wp) :: boundary, previous_coordinate, current_coordinate
    REAL(wp) :: fraction, intersection_x, intersection_y, area_twice
    REAL(wp) :: denominator
    INTEGER :: cell_x, cell_y, vertex, next_vertex, edge
    INTEGER :: vertex_count, clipped_count, previous_vertex, axis
    LOGICAL :: keep_greater, previous_inside, current_inside

    cell_fraction = 0.0_wp
    segment_x = endpoints_x(2) - endpoints_x(1)
    segment_y = endpoints_y(2) - endpoints_y(1)
    segment_length = SQRT(segment_x**2 + segment_y**2)
    IF (segment_length .LE. EPSILON(1.0_wp)) RETURN

    normal_x = -segment_y / segment_length
    normal_y = segment_x / segment_length
    half_width = 0.5_wp * width

    min_x = MINVAL(endpoints_x) - half_width * ABS(normal_x)
    max_x = MAXVAL(endpoints_x) + half_width * ABS(normal_x)
    min_y = MINVAL(endpoints_y) - half_width * ABS(normal_y)
    max_y = MAXVAL(endpoints_y) + half_width * ABS(normal_y)

    DO cell_y = 1, comp_cells_y
       IF (y_stag(cell_y+1) .LT. min_y .OR. y_stag(cell_y) .GT. max_y) CYCLE

       DO cell_x = 1, comp_cells_x
          IF (x_stag(cell_x+1) .LT. min_x .OR. x_stag(cell_x) .GT. max_x) CYCLE

          vertex_count = 4
          polygon_x(1) = endpoints_x(1) + normal_x * half_width
          polygon_y(1) = endpoints_y(1) + normal_y * half_width
          polygon_x(2) = endpoints_x(2) + normal_x * half_width
          polygon_y(2) = endpoints_y(2) + normal_y * half_width
          polygon_x(3) = endpoints_x(2) - normal_x * half_width
          polygon_y(3) = endpoints_y(2) - normal_y * half_width
          polygon_x(4) = endpoints_x(1) - normal_x * half_width
          polygon_y(4) = endpoints_y(1) - normal_y * half_width

          ! Clip the fissure polygon against the cell's four axis-aligned edges.
          DO edge = 1, 4
             SELECT CASE (edge)
             CASE (1)
                axis = 1
                boundary = x_stag(cell_x)
                keep_greater = .TRUE.
             CASE (2)
                axis = 1
                boundary = x_stag(cell_x+1)
                keep_greater = .FALSE.
             CASE (3)
                axis = 2
                boundary = y_stag(cell_y)
                keep_greater = .TRUE.
             CASE (4)
                axis = 2
                boundary = y_stag(cell_y+1)
                keep_greater = .FALSE.
             END SELECT

             clipped_count = 0
             previous_vertex = vertex_count
             IF (axis .EQ. 1) THEN
                previous_coordinate = polygon_x(previous_vertex)
             ELSE
                previous_coordinate = polygon_y(previous_vertex)
             END IF
             IF (keep_greater) THEN
                previous_inside = previous_coordinate .GE. boundary
             ELSE
                previous_inside = previous_coordinate .LE. boundary
             END IF

             DO vertex = 1, vertex_count
                IF (axis .EQ. 1) THEN
                   current_coordinate = polygon_x(vertex)
                ELSE
                   current_coordinate = polygon_y(vertex)
                END IF
                IF (keep_greater) THEN
                   current_inside = current_coordinate .GE. boundary
                ELSE
                   current_inside = current_coordinate .LE. boundary
                END IF

                IF (current_inside .NEQV. previous_inside) THEN
                   IF (axis .EQ. 1) THEN
                      denominator = polygon_x(vertex) - polygon_x(previous_vertex)
                      fraction = (boundary - polygon_x(previous_vertex)) / denominator
                      intersection_x = boundary
                      intersection_y = polygon_y(previous_vertex) + fraction * &
                           (polygon_y(vertex) - polygon_y(previous_vertex))
                   ELSE
                      denominator = polygon_y(vertex) - polygon_y(previous_vertex)
                      fraction = (boundary - polygon_y(previous_vertex)) / denominator
                      intersection_x = polygon_x(previous_vertex) + fraction * &
                           (polygon_x(vertex) - polygon_x(previous_vertex))
                      intersection_y = boundary
                   END IF
                   clipped_count = clipped_count + 1
                   clipped_x(clipped_count) = intersection_x
                   clipped_y(clipped_count) = intersection_y
                END IF

                IF (current_inside) THEN
                   clipped_count = clipped_count + 1
                   clipped_x(clipped_count) = polygon_x(vertex)
                   clipped_y(clipped_count) = polygon_y(vertex)
                END IF

                previous_vertex = vertex
                previous_inside = current_inside
             END DO

             vertex_count = clipped_count
             IF (vertex_count .EQ. 0) EXIT
             polygon_x(1:vertex_count) = clipped_x(1:vertex_count)
             polygon_y(1:vertex_count) = clipped_y(1:vertex_count)
          END DO

          IF (vertex_count .LT. 3) CYCLE

          area_twice = 0.0_wp
          DO vertex = 1, vertex_count
             next_vertex = MOD(vertex, vertex_count) + 1
             area_twice = area_twice + polygon_x(vertex) * polygon_y(next_vertex) &
                  - polygon_x(next_vertex) * polygon_y(vertex)
          END DO
          cell_fraction(cell_x,cell_y) = MIN(1.0_wp,                          &
               0.5_wp * ABS(area_twice) / (dx * dy))
       END DO
    END DO

  END SUBROUTINE compute_cell_fissure_fraction

  !> \brief Estimate cell coverage by an elliptical source using subcell sampling.
  !>
  !> \param[in] xs X coordinate of the ellipse centre [m].
  !> \param[in] ys Y coordinate of the ellipse centre [m].
  !> \param[in] rs First semi-axis of the elliptical source [m].
  !> \param[in] r2s Second semi-axis of the elliptical source [m].
  !> \param[in] angles Clockwise angle of the ellipse relative to the x axis [degrees].
  !> \param[out] cell_fract Source-covered fraction of each computational cell.

  SUBROUTINE compute_cell_fract(xs,ys,rs,r2s,angles,cell_fract)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: xs,ys,rs,r2s,angles

    REAL(wp), INTENT(OUT) :: cell_fract(comp_cells_x,comp_cells_y)

    REAL(wp), ALLOCATABLE :: x_subgrid(:) , y_subgrid(:)

    INTEGER, ALLOCATABLE :: check_subgrid(:)

    INTEGER n_points , n_points2

    REAL(wp) :: source_area

    INTEGER :: h ,j, k

    REAL(wp) :: ang_rad

    ang_rad = angles / 180.0_wp * ATAN(1.0_wp)*4.0_wp

    WRITE(*,*) 'ang_rad',ang_rad

    n_points = 200
    n_points2 = n_points**2

    ALLOCATE( x_subgrid(n_points2) )
    ALLOCATE( y_subgrid(n_points2) )
    ALLOCATE( check_subgrid(n_points2) )

    x_subgrid = 0.0_wp
    y_subgrid = 0.0_wp

    DO h = 1,n_points

       x_subgrid(h:n_points2:n_points) = DBLE(h)
       y_subgrid((h-1)*n_points+1:h*n_points) = DBLE(h)

    END DO

    x_subgrid = ( 2.0_wp * x_subgrid - 1.0_wp ) / ( 2.0_wp * DBLE(n_points) )
    y_subgrid = ( 2.0_wp * y_subgrid - 1.0_wp ) / ( 2.0_wp * DBLE(n_points) )

    x_subgrid = ( x_subgrid - 0.5_wp ) * dx
    y_subgrid = ( y_subgrid - 0.5_wp ) * dy

    DO j=1,comp_cells_x

       DO k=1,comp_cells_y

          IF ( ( x_stag(j+1) .LT. ( xs - MAX(rs,r2s) ) ) .OR.                   &
               ( x_stag(j) .GT. ( xs + MAX(rs,r2s) ) ) .OR.                     &
               ( y_stag(k+1) .LT. ( ys - MAX(rs,r2s) ) ) .OR.                   &
               ( y_stag(k) .GT. ( ys + MAX(rs,r2s) ) ) ) THEN

             cell_fract(j,k) = 0.0_wp 

          ELSE

             check_subgrid = 0

             WHERE ( ( ( ( x_comp(j) + x_subgrid - xs )*COS(ang_rad) +          &
                  ( y_comp(k) + y_subgrid - ys )*SIN(ang_rad) )**2 / rs**2 +   &
                  ( - ( x_comp(j) + x_subgrid - xs )*SIN(ang_rad) +            &
                  ( y_comp(k) + y_subgrid - ys )*COS(ang_rad) )**2 / r2s**2 )  &
                  .LE. 1.0_wp )

                check_subgrid = 1

             END WHERE

             cell_fract(j,k) = REAL(SUM(check_subgrid),wp) / &
                  REAL(n_points2,wp)

          END IF

       ENDDO

    ENDDO

    source_area = dx*dy*SUM(cell_fract)

    IF ( VERBOSE_LEVEL .GE. 0 ) THEN

       WRITE(*,*) 'Source area =',source_area,' Error =',ABS( 1.0_wp -          &
            dx*dy*SUM(cell_fract) / ( 4.0_wp*ATAN(1.0_wp)*rs*r2s ) )

    END IF

    DEALLOCATE( x_subgrid )
    DEALLOCATE( y_subgrid )
    DEALLOCATE( check_subgrid )

    RETURN

  END SUBROUTINE compute_cell_fract

END MODULE geometry_2d
