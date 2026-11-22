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

        INTEGER :: i,j,k,n,km,kp,direction,j_neighbor,k_neighbor
        INTEGER :: row_count(comp_cells_y),row_offset(comp_cells_y+1)
        INTEGER :: hp_row_offset(comp_cells_y+1)
        INTEGER :: row_count_x(comp_cells_y),row_offset_x(comp_cells_y+1)
        INTEGER :: row_count_y(comp_interfaces_y),row_offset_y(comp_interfaces_y+1)
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

        !$OMP END PARALLEL

        ! The final mask is compacted by row. Each row owns a disjoint
        ! output interval; the old lexicographic (k,j) ordering is preserved.
        row_offset(1)=0
        DO k=1,comp_cells_y
            row_offset(k+1)=row_offset(k)+row_count(k)
        END DO

        this%solve_cells=row_offset(comp_cells_y+1)

        ! Build dependent worksets in one OpenMP team. Each source-row owns
        ! disjoint output intervals; serial prefix sums preserve (k,j) order.
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,n,active,direction,j_neighbor,k_neighbor)

        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=row_offset(k)
            DO j=1,comp_cells_x
                IF (.NOT.this%solve_mask(j,k)) CYCLE
                n=n+1
                this%j_cent(n)=j
                this%k_cent(n)=k
            END DO
            IF (n/=row_offset(k+1)) CALL fatal_error('Lie OpenMP: inconsistent solve worklist')
        END DO
        !$OMP END DO

        ! The previous HP list contains unique cells, so its old membership
        ! flags can be cleared independently before the new list is emitted.
        !$OMP DO SCHEDULE(STATIC)
        DO i=1,this%hp_eta_cells
            this%hp_eta_mask(this%hp_eta_j(i),this%hp_eta_k(i))=.FALSE.
        END DO
        !$OMP END DO

        ! Preserve the exact HP list prefix: every solve cell in (k,j) order.
        !$OMP DO SCHEDULE(STATIC)
        DO i=1,this%solve_cells
            j=this%j_cent(i)
            k=this%k_cent(i)
            this%hp_eta_j(i)=j
            this%hp_eta_k(i)=k
            this%hp_eta_mask(j,k)=.TRUE.
        END DO
        !$OMP END DO

        ! The first adjacent solve cell in (k,j) order owns each non-solve HP
        ! neighbour. Counts and offsets allow deterministic emission without
        ! races, atomics or temporary per-thread lists.
        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=0
            DO i=row_offset(k)+1,row_offset(k+1)
                j=this%j_cent(i)
                DO direction=1,4
                    IF (first_hp_neighbor(j,k,direction)) n=n+1
                END DO
            END DO
            row_count(k)=n
        END DO
        !$OMP END DO

        !$OMP SINGLE
        hp_row_offset(1)=0
        DO k=1,comp_cells_y
            hp_row_offset(k+1)=hp_row_offset(k)+row_count(k)
        END DO
        this%hp_eta_cells=this%solve_cells+hp_row_offset(comp_cells_y+1)
        IF (this%hp_eta_cells>SIZE(this%hp_eta_j)) &
            CALL fatal_error('Lie OpenMP: HP worklist exceeds domain capacity')
        !$OMP END SINGLE

        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=this%solve_cells+hp_row_offset(k)
            DO i=row_offset(k)+1,row_offset(k+1)
                j=this%j_cent(i)
                DO direction=1,4
                    IF (.NOT.first_hp_neighbor(j,k,direction)) CYCLE
                    j_neighbor=j
                    k_neighbor=k
                    SELECT CASE(direction)
                    CASE(1)
                        j_neighbor=j-1
                    CASE(2)
                        j_neighbor=j+1
                    CASE(3)
                        k_neighbor=k-1
                    CASE(4)
                        k_neighbor=k+1
                    END SELECT
                    n=n+1
                    this%hp_eta_j(n)=j_neighbor
                    this%hp_eta_k(n)=k_neighbor
                    this%hp_eta_mask(j_neighbor,k_neighbor)=.TRUE.
                END DO
            END DO
            IF (n/=this%solve_cells+hp_row_offset(k+1)) &
                CALL fatal_error('Lie OpenMP: inconsistent HP worklist')
        END DO
        !$OMP END DO

        ! Face membership is the union of its adjacent solve cells. Computing
        ! each destination exactly once avoids the overlapping serial writes
        ! previously used to mark faces from the compact solve list. Fuse the
        ! first reconstruction halo with the face-mask scans and row counts.
        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_interfaces_y
            IF (k<=comp_cells_y) THEN
                n=0
                DO j=1,comp_interfaces_x
                    active=.FALSE.
                    IF (j>1) active=this%solve_mask(j-1,k)
                    IF (j<=comp_cells_x) active=active .OR. this%solve_mask(j,k)
                    this%solve_mask_x(j,k)=active
                    IF (active) n=n+1
                END DO
                row_count_x(k)=n

                DO j=1,comp_cells_x
                    active=this%solve_mask(j,k)
                    IF (j>1) active=active .OR. this%solve_mask(j-1,k)
                    IF (j<comp_cells_x) active=active .OR. this%solve_mask(j+1,k)
                    IF (k>1) active=active .OR. this%solve_mask(j,k-1)
                    IF (k<comp_cells_y) active=active .OR. this%solve_mask(j,k+1)
                    this%reconstruction_mask(j,k)=active
                END DO
            END IF

            n=0
            DO j=1,comp_cells_x
                active=.FALSE.
                IF (k>1) active=this%solve_mask(j,k-1)
                IF (k<=comp_cells_y) active=active .OR. this%solve_mask(j,k)
                this%solve_mask_y(j,k)=active
                IF (active) n=n+1
            END DO
            row_count_y(k)=n
        END DO
        !$OMP END DO

        ! The read-only reconstruction halo has exactly two four-neighbour
        ! rings. Reuse the existing temporary mask as the second destination,
        ! rather than copying the entire reconstruction mask twice.
        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            DO j=1,comp_cells_x
                active=this%reconstruction_mask(j,k)
                IF (j>1) active=active .OR. this%reconstruction_mask(j-1,k)
                IF (j<comp_cells_x) active=active .OR. this%reconstruction_mask(j+1,k)
                IF (k>1) active=active .OR. this%reconstruction_mask(j,k-1)
                IF (k<comp_cells_y) active=active .OR. this%reconstruction_mask(j,k+1)
                this%solve_mask_temp(j,k)=active
            END DO
        END DO
        !$OMP END DO

        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_cells_y
            n=0
            DO j=1,comp_cells_x
                active=this%solve_mask_temp(j,k)
                this%reconstruction_mask(j,k)=active
                IF (active) n=n+1
            END DO
            row_count(k)=n
        END DO
        !$OMP END DO

        !$OMP SINGLE
        row_offset(1)=0
        row_offset_x(1)=0
        row_offset_y(1)=0
        DO k=1,comp_cells_y
            row_offset(k+1)=row_offset(k)+row_count(k)
            row_offset_x(k+1)=row_offset_x(k)+row_count_x(k)
            row_offset_y(k+1)=row_offset_y(k)+row_count_y(k)
        END DO
        row_offset_y(comp_interfaces_y+1)=row_offset_y(comp_cells_y+1)+row_count_y(comp_interfaces_y)
        this%reconstruction_cells=row_offset(comp_cells_y+1)
        this%solve_interfaces_x=row_offset_x(comp_cells_y+1)
        this%solve_interfaces_y=row_offset_y(comp_interfaces_y+1)
        !$OMP END SINGLE

        !$OMP DO SCHEDULE(STATIC)
        DO k=1,comp_interfaces_y
            IF (k<=comp_cells_y) THEN
                n=row_offset(k)
                DO j=1,comp_cells_x
                    IF (.NOT.this%reconstruction_mask(j,k)) CYCLE
                    n=n+1
                    this%j_reconstruction(n)=j
                    this%k_reconstruction(n)=k
                END DO
                IF (n/=row_offset(k+1)) CALL fatal_error('Lie OpenMP: inconsistent reconstruction worklist')

                n=row_offset_x(k)
                DO j=1,comp_interfaces_x
                    IF (.NOT.this%solve_mask_x(j,k)) CYCLE
                    n=n+1
                    this%j_stag_x(n)=j
                    this%k_stag_x(n)=k
                END DO
                IF (n/=row_offset_x(k+1)) CALL fatal_error('Lie OpenMP: inconsistent x-face worklist')
            END IF

            n=row_offset_y(k)
            DO j=1,comp_cells_x
                IF (.NOT.this%solve_mask_y(j,k)) CYCLE
                n=n+1
                this%j_stag_y(n)=j
                this%k_stag_y(n)=k
            END DO
            IF (n/=row_offset_y(k+1)) CALL fatal_error('Lie OpenMP: inconsistent y-face worklist')
        END DO
        !$OMP END DO

        !$OMP END PARALLEL

        ! Accepted Lie publication is consumed; future direct calls must
        ! rebuild wetness from the caller's conservative states.
        this%committed_mask_ready=.FALSE.
        IF (solve_all .OR. rebuild_from_state .OR. release_due) THEN
            this%next_release_time=MINVAL(this%solve_mask_time, &
                MASK=this%solve_mask_time>t)
        END IF

    CONTAINS

        !> @brief Decide which solve cell first contributes a particular HP neighbour.
        !> @details Source cells are ordered by (k,j), with neighbours W,E,S,N.
        !>          Every non-solve neighbour is emitted by its earliest adjacent
        !>          source, preserving the historical HP list without duplicates.
        !> @param[in] j Source x index.
        !> @param[in] k Source y index.
        !> @param[in] direction Neighbour order: 1=W, 2=E, 3=S, 4=N.
        !> @return True only for the first source that reaches this neighbour.
        PURE LOGICAL FUNCTION first_hp_neighbor(j,k,direction) RESULT(first)
            INTEGER, INTENT(IN) :: j,k,direction

            first=.FALSE.
            SELECT CASE(direction)
            CASE(1)
                IF (j<=1) RETURN
                IF (this%solve_mask(j-1,k)) RETURN
                IF (k>1) THEN
                    IF (this%solve_mask(j-1,k-1)) RETURN
                END IF
                IF (j>2) THEN
                    IF (this%solve_mask(j-2,k)) RETURN
                END IF
            CASE(2)
                IF (j>=comp_cells_x) RETURN
                IF (this%solve_mask(j+1,k)) RETURN
                IF (k>1) THEN
                    IF (this%solve_mask(j+1,k-1)) RETURN
                END IF
            CASE(3)
                IF (k<=1) RETURN
                IF (this%solve_mask(j,k-1)) RETURN
                IF (k>2) THEN
                    IF (this%solve_mask(j,k-2)) RETURN
                END IF
                IF (j>1) THEN
                    IF (this%solve_mask(j-1,k-1)) RETURN
                END IF
                IF (j<comp_cells_x) THEN
                    IF (this%solve_mask(j+1,k-1)) RETURN
                END IF
            CASE(4)
                IF (k>=comp_cells_y) RETURN
                IF (this%solve_mask(j,k+1)) RETURN
            CASE DEFAULT
                RETURN
            END SELECT
            first=.TRUE.
        END FUNCTION first_hp_neighbor
    END SUBROUTINE check_solve

END MODULE domain_2d
