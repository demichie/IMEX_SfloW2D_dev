!********************************************************************************
!> \brief Active computational-domain management
!>
!> This module owns the masks and compact index lists that identify cells and
!> interfaces requiring numerical work.
!>
!> Builds compact cell/face worksets from wet cells, inlet cells and stencil halos. Spatial kernels
!> consume these lists rather than scanning the entire domain.
!********************************************************************************

MODULE domain_2d
    USE parameters_2d, ONLY : idx_totMassEqn

    USE parameters_2d, ONLY : wp
    USE parameters_2d, ONLY : radial_source_flag, bottom_radial_source_flag,                                           &
        bottom_fissural_source_flag

    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y, comp_cells_xy
    USE geometry_2d, ONLY : comp_interfaces_x, comp_interfaces_y
    USE geometry_2d, ONLY : source_cell, cell_source_fractions
    USE diagnostics_2d, ONLY : fatal_error

    IMPLICIT NONE

    PRIVATE

    ! Two spatial cell layers cover the face and three-point reconstruction stencil.
    ! This width is independent of any temporal Runge-Kutta stage count.
    INTEGER, PARAMETER :: active_halo_layers = 2

    !> \brief Active-cell/face indices and the wider reconstruction halo.
    !> \details Solve lists determine where equations are advanced; reconstruction
    !>          lists also contain neighbouring dry cells read by the HP stencil.
    TYPE, PUBLIC :: domain_type

        REAL(wp), ALLOCATABLE :: solve_mask_time(:,:)
        !> Next pending initialization release; refreshed only when its time is reached.
        REAL(wp) :: next_release_time = HUGE(0.0_wp)
        !> Accepted Lie published wetness for the previous solve worklist.
        LOGICAL :: committed_mask_ready = .FALSE.

        LOGICAL, ALLOCATABLE :: solve_mask(:,:)
        LOGICAL, ALLOCATABLE :: solve_mask_temp(:,:)
        LOGICAL, ALLOCATABLE :: reconstruction_mask(:,:)
        LOGICAL, ALLOCATABLE :: solve_mask_x(:,:)
        LOGICAL, ALLOCATABLE :: solve_mask_y(:,:)

        INTEGER :: solve_cells
        INTEGER :: reconstruction_cells
        INTEGER :: hp_eta_cells
        INTEGER :: solve_interfaces_x
        INTEGER :: solve_interfaces_y

        INTEGER, ALLOCATABLE :: j_cent(:)
        INTEGER, ALLOCATABLE :: k_cent(:)
        LOGICAL, ALLOCATABLE :: hp_eta_mask(:,:)
        INTEGER, ALLOCATABLE :: hp_eta_j(:), hp_eta_k(:)
        INTEGER, ALLOCATABLE :: j_reconstruction(:)
        INTEGER, ALLOCATABLE :: k_reconstruction(:)
        INTEGER, ALLOCATABLE :: j_stag_x(:)
        INTEGER, ALLOCATABLE :: k_stag_x(:)
        INTEGER, ALLOCATABLE :: j_stag_y(:)
        INTEGER, ALLOCATABLE :: k_stag_y(:)

    CONTAINS

        PROCEDURE :: initialize => initialize_domain
        PROCEDURE :: finalize => finalize_domain
        PROCEDURE :: check_solve

    END TYPE domain_type

CONTAINS

    !> \brief Allocate and reset active-cell, halo and face worksets.
    !>
    !> \param[in,out] this Active-domain masks and compact cell/face worksets.

    SUBROUTINE initialize_domain(this)

        CLASS(domain_type), INTENT(INOUT) :: this

        ALLOCATE( this%solve_mask_time(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%solve_mask(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%solve_mask_temp(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%reconstruction_mask(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%hp_eta_mask(comp_cells_x,comp_cells_y) )
        ALLOCATE( this%solve_mask_x(comp_interfaces_x,comp_cells_y) )
        ALLOCATE( this%solve_mask_y(comp_cells_x,comp_interfaces_y) )

        this%solve_mask_time = 0.0_wp
        this%next_release_time = HUGE(0.0_wp)
        this%committed_mask_ready = .FALSE.

        this%solve_mask = .FALSE.
        this%solve_mask_temp = .FALSE.
        this%reconstruction_mask = .FALSE.
        this%hp_eta_mask = .FALSE.
        this%solve_mask_x = .FALSE.
        this%solve_mask_y = .FALSE.

        this%solve_mask(1,1:comp_cells_y) = .TRUE.
        this%solve_mask(comp_cells_x,1:comp_cells_y) = .TRUE.
        this%solve_mask(1:comp_cells_x,1) = .TRUE.
        this%solve_mask(1:comp_cells_x,comp_cells_y) = .TRUE.

        comp_cells_xy = comp_cells_x * comp_cells_y

        ALLOCATE( this%j_cent(comp_cells_xy) )
        ALLOCATE( this%k_cent(comp_cells_xy) )
        ALLOCATE( this%hp_eta_j(comp_cells_xy), this%hp_eta_k(comp_cells_xy) )
        ALLOCATE( this%j_reconstruction(comp_cells_xy) )
        ALLOCATE( this%k_reconstruction(comp_cells_xy) )
        ALLOCATE( this%j_stag_x(comp_interfaces_x*comp_cells_y) )
        ALLOCATE( this%k_stag_x(comp_interfaces_x*comp_cells_y) )
        ALLOCATE( this%j_stag_y(comp_cells_x*comp_interfaces_y) )
        ALLOCATE( this%k_stag_y(comp_cells_x*comp_interfaces_y) )

        this%solve_cells = 0
        this%reconstruction_cells = 0
        this%hp_eta_cells = 0
        this%solve_interfaces_x = 0
        this%solve_interfaces_y = 0

        this%j_cent = 0
        this%k_cent = 0
        this%j_reconstruction = 0
        this%k_reconstruction = 0
        this%j_stag_x = 0
        this%k_stag_x = 0
        this%j_stag_y = 0
        this%k_stag_y = 0

    END SUBROUTINE initialize_domain

    !> \brief Release active-domain masks and compact index lists.
    !>
    !> \param[in,out] this Active-domain masks and compact cell/face worksets.

    SUBROUTINE finalize_domain(this)

        CLASS(domain_type), INTENT(INOUT) :: this

        DEALLOCATE( this%solve_mask_time )
        DEALLOCATE( this%solve_mask )
        DEALLOCATE( this%solve_mask_temp )
        DEALLOCATE( this%reconstruction_mask )
        DEALLOCATE( this%hp_eta_mask, this%hp_eta_j, this%hp_eta_k )
        DEALLOCATE( this%solve_mask_x )
        DEALLOCATE( this%solve_mask_y )

        DEALLOCATE( this%j_cent, this%k_cent )
        DEALLOCATE( this%j_reconstruction, this%k_reconstruction )
        DEALLOCATE( this%j_stag_x, this%k_stag_x )
        DEALLOCATE( this%j_stag_y, this%k_stag_y )

    END SUBROUTINE finalize_domain

    !> \brief Rebuild active-cell and face worksets from the current flow and source geometry.
    !>
    !> \param[in,out] this Active-domain masks and compact cell/face worksets.
    !> \param[in] q Cell-centered conservative states, indexed as (variable,x-cell,y-cell).
    !> \param[in] t Current simulation or stage time [s].
    !> \param[in] solve_all When true, include the complete domain instead of only flow/source cells
    !>                      and halos.
    !> @param[in] q_other Optional second-layer conservative state for a shared
    !>                  coupled workset. Omitted preserves the original domain.

    SUBROUTINE check_solve(this, q, t, solve_all, q_other)

    IMPLICIT NONE

        CLASS(domain_type), INTENT(INOUT) :: this
        REAL(wp), INTENT(IN) :: q(:,:,:)
        REAL(wp), INTENT(IN) :: t
        LOGICAL, INTENT(IN) :: solve_all
        REAL(wp), INTENT(IN), OPTIONAL :: q_other(:,:,:)

        INTEGER :: i,j,k,n,km,kp
        INTEGER :: row_count(comp_cells_y),row_offset(comp_cells_y+1)
        LOGICAL :: active,rebuild_from_state,release_due

        rebuild_from_state=.NOT.this%committed_mask_ready
        release_due=this%next_release_time<=t

        IF (PRESENT(q_other)) THEN
            IF (SIZE(q_other,1)<idx_totMassEqn .OR. &
                SIZE(q_other,2)/=comp_cells_x .OR. SIZE(q_other,3)/=comp_cells_y) &
                CALL fatal_error('Coupled workset: invalid second-layer state shape')
        END IF

        ! On a normal step, accepted Lie publication has already written wetness
        ! for the old worklist. The full-state path is needed only without that
        ! publication (initialization, restart, and direct callers), or when a
        ! deferred initial release reaches its activation time.
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,active,km,kp)

        IF (solve_all) THEN

            !$OMP WORKSHARE
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1)=.TRUE.
            !$OMP END WORKSHARE

        ELSEIF (rebuild_from_state) THEN

            !$OMP WORKSHARE
            this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1)=.FALSE.
            !$OMP END WORKSHARE

        ELSE

            ! A halo cell processed by Lie may still contain delayed mass.
            ! Apply the current time gate to the previous compact worklist.
            !$OMP DO SCHEDULE(STATIC)
            DO i=1,this%solve_cells
                j=this%j_cent(i)
                k=this%k_cent(i)
                IF (j<=1 .OR. j>=comp_cells_x .OR. k<=1 .OR. k>=comp_cells_y) CYCLE
                IF (this%solve_mask_time(j,k)>t) this%solve_mask(j,k)=.FALSE.
            END DO
            !$OMP END DO

        END IF

        ! Mass scans are skipped in ordinary Lie-to-Lie steps. When a release
        ! time is reached, recover any newly eligible wet cells outside the
        ! previous worklist. The release-time map is set at initialization.
        IF (.NOT.solve_all .AND. (rebuild_from_state .OR. release_due)) THEN

            !$OMP WORKSHARE
            WHERE ((q(idx_totMassEqn,2:comp_cells_x-1,2:comp_cells_y-1)>0.0_wp) .AND. &
                (this%solve_mask_time(2:comp_cells_x-1,2:comp_cells_y-1)<=t)) &
                this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1)=.TRUE.
            !$OMP END WORKSHARE

            IF (PRESENT(q_other)) THEN

                !$OMP WORKSHARE
                WHERE ((q_other(idx_totMassEqn,2:comp_cells_x-1,2:comp_cells_y-1)>0.0_wp) .AND. &
                    (this%solve_mask_time(2:comp_cells_x-1,2:comp_cells_y-1)<=t)) &
                    this%solve_mask(2:comp_cells_x-1,2:comp_cells_y-1)=.TRUE.
                !$OMP END WORKSHARE

            END IF

        END IF

        IF (bottom_radial_source_flag .OR. bottom_fissural_source_flag) THEN

            !$OMP WORKSHARE
            WHERE (cell_source_fractions>0.0_wp) this%solve_mask=.TRUE.
            !$OMP END WORKSHARE

        END IF

        IF (radial_source_flag) THEN

            !$OMP DO SCHEDULE(STATIC)
            DO k=2,comp_cells_y-1
                DO j=2,comp_cells_x-1
                    IF (source_cell(j,k)==2) this%solve_mask(j,k)=.TRUE.
                END DO
            END DO
            !$OMP END DO

        END IF

        ! Ping-pong dilation: preserve the unmodified perimeter once, then
        ! write all interior destinations into the existing temporary mask.
        !$OMP SINGLE
        this%solve_mask_temp(1,:)=this%solve_mask(1,:)
        this%solve_mask_temp(comp_cells_x,:)=this%solve_mask(comp_cells_x,:)
        this%solve_mask_temp(:,1)=this%solve_mask(:,1)
        this%solve_mask_temp(:,comp_cells_y)=this%solve_mask(:,comp_cells_y)
        !$OMP END SINGLE

        !$OMP DO SCHEDULE(STATIC)
        DO k=2,comp_cells_y-1
            DO j=2,comp_cells_x-1
                this%solve_mask_temp(j,k)=this%solve_mask(j,k) .OR. &
                    this%solve_mask(j-1,k) .OR. this%solve_mask(j+1,k) .OR. &
                    this%solve_mask(j,k-1) .OR. this%solve_mask(j,k+1)
            END DO
        END DO
        !$OMP END DO

        ! Second dilation reads only the completed temporary mask. Apply the
        ! radial exclusion after the halo, and count each entire finalized row.
        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=0
            IF (k==1 .OR. k==comp_cells_y .OR. comp_cells_x<=2) THEN
                DO j=1,comp_cells_x
                    active=this%solve_mask_temp(j,k)
                    IF (radial_source_flag) THEN
                        IF (source_cell(j,k)==1) active=.FALSE.
                    END IF
                    this%solve_mask(j,k)=active
                    IF (active) n=n+1
                END DO
            ELSE
                km=k-1
                kp=k+1
                active=this%solve_mask_temp(1,k)
                IF (radial_source_flag) THEN
                    IF (source_cell(1,k)==1) active=.FALSE.
                END IF
                this%solve_mask(1,k)=active
                IF (active) n=n+1

                DO j=2,comp_cells_x-1
                    active=this%solve_mask_temp(j,k) .OR. &
                        this%solve_mask_temp(j-1,k) .OR. this%solve_mask_temp(j+1,k) .OR. &
                        this%solve_mask_temp(j,km) .OR. this%solve_mask_temp(j,kp)
                    IF (radial_source_flag) THEN
                        IF (source_cell(j,k)==1) active=.FALSE.
                    END IF
                    this%solve_mask(j,k)=active
                    IF (active) n=n+1
                END DO

                active=this%solve_mask_temp(comp_cells_x,k)
                IF (radial_source_flag) THEN
                    IF (source_cell(comp_cells_x,k)==1) active=.FALSE.
                END IF
                this%solve_mask(comp_cells_x,k)=active
                IF (active) n=n+1
            END IF
            row_count(k)=n
        END DO
        !$OMP END DO

        !$OMP WORKSHARE
        this%solve_mask_x=.FALSE.
        this%solve_mask_y=.FALSE.
        !$OMP END WORKSHARE

        !$OMP END PARALLEL

        ! The final mask is compacted by row. Each row owns a disjoint
        ! output interval; the old lexicographic (k,j) ordering is preserved.
        row_offset(1)=0
        DO k=1,comp_cells_y
            row_offset(k+1)=row_offset(k)+row_count(k)
        END DO

        this%solve_cells=row_offset(comp_cells_y+1)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,n) SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=row_offset(k)
            DO j=1,comp_cells_x
                IF (.NOT.this%solve_mask(j,k)) CYCLE
                n=n+1
                this%j_cent(n)=j
                this%k_cent(n)=k
            END DO
        END DO
        !$OMP END PARALLEL DO

        ! Adjacent solve cells share faces. Keep face writes out of the
        ! parallel emitter to avoid overlapping OpenMP writes.
        DO i=1,this%solve_cells
            j=this%j_cent(i)
            k=this%k_cent(i)
            this%solve_mask_x(j,k)=.TRUE.
            this%solve_mask_x(j+1,k)=.TRUE.
            this%solve_mask_y(j,k)=.TRUE.
            this%solve_mask_y(j,k+1)=.TRUE.
        END DO

        IF (this%solve_cells/=COUNT(this%solve_mask)) &
            CALL fatal_error('Lie OpenMP: inconsistent unique solve worklist')

        ! A later refresh without new Lie publication must rebuild wetness
        ! from the caller's current conservative state, not from halo values.
        this%committed_mask_ready=.FALSE.
        IF (solve_all .OR. rebuild_from_state .OR. release_due) THEN
            this%next_release_time=MINVAL(this%solve_mask_time, &
                MASK=this%solve_mask_time>t)
        END IF

        ! Preserve the former HP workset order: active cells, then neighbours.
        DO i = 1, this%hp_eta_cells
            this%hp_eta_mask(this%hp_eta_j(i),this%hp_eta_k(i)) = .FALSE.

        END DO

        this%hp_eta_cells = 0

        DO i = 1, this%solve_cells
            CALL append_hp_cell(this%j_cent(i),this%k_cent(i))

        END DO

        DO i = 1, this%solve_cells
            j = this%j_cent(i)
            k = this%k_cent(i)
            IF (j .GT. 1) CALL append_hp_cell(j-1,k)
            IF (j .LT. comp_cells_x) CALL append_hp_cell(j+1,k)
            IF (k .GT. 1) CALL append_hp_cell(j,k-1)
            IF (k .LT. comp_cells_y) CALL append_hp_cell(j,k+1)

        END DO

        ! HP reconstruction reads a two-cell stencil around every evolved cell.
        ! Keep this read-only halo separate from solve_mask so no additional cells
        ! are advanced by the time integrator.
        this%reconstruction_mask = this%solve_mask

        DO i = 1,2
            this%solve_mask_temp = this%reconstruction_mask

            IF ( comp_cells_x .GT. 1 ) THEN
                this%reconstruction_mask(2:comp_cells_x,:) =                                                           &
                    this%reconstruction_mask(2:comp_cells_x,:) .OR.                                                    &
                    this%solve_mask_temp(1:comp_cells_x-1,:)
                this%reconstruction_mask(1:comp_cells_x-1,:) =                                                         &
                    this%reconstruction_mask(1:comp_cells_x-1,:) .OR.                                                  &
                    this%solve_mask_temp(2:comp_cells_x,:)

            END IF

            IF ( comp_cells_y .GT. 1 ) THEN
                this%reconstruction_mask(:,2:comp_cells_y) =                                                           &
                    this%reconstruction_mask(:,2:comp_cells_y) .OR.                                                    &
                    this%solve_mask_temp(:,1:comp_cells_y-1)
                this%reconstruction_mask(:,1:comp_cells_y-1) =                                                         &
                    this%reconstruction_mask(:,1:comp_cells_y-1) .OR.                                                  &
                    this%solve_mask_temp(:,2:comp_cells_y)

            END IF

        END DO

        i = 0

        DO k = 1,comp_cells_y

            DO j = 1,comp_cells_x

                IF ( this%reconstruction_mask(j,k) ) THEN
                    i = i+1
                    this%j_reconstruction(i) = j
                    this%k_reconstruction(i) = k

                END IF

            END DO

        END DO

        this%reconstruction_cells = i
        IF (this%reconstruction_cells/=COUNT(this%reconstruction_mask))                                                &
            CALL fatal_error('Lie OpenMP: inconsistent unique reconstruction worklist')

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

    CONTAINS
    SUBROUTINE append_hp_cell(j_cell,k_cell)
        INTEGER, INTENT(IN) :: j_cell,k_cell
        IF (this%hp_eta_mask(j_cell,k_cell)) RETURN
        this%hp_eta_mask(j_cell,k_cell) = .TRUE.
        this%hp_eta_cells = this%hp_eta_cells + 1
        this%hp_eta_j(this%hp_eta_cells) = j_cell
        this%hp_eta_k(this%hp_eta_cells) = k_cell
    END SUBROUTINE append_hp_cell
    END SUBROUTINE check_solve

END MODULE domain_2d
