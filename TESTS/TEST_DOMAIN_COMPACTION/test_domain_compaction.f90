!> @brief Differential tests and performance experiments for serial and OpenMP compaction.
PROGRAM test_domain_compaction
    USE omp_lib, ONLY: omp_set_dynamic,omp_set_num_threads,omp_get_wtime
    USE domain_compaction_test_kernels, ONLY: compact_reference,compact_global_parallel,compact_sparse_parallel
    IMPLICIT NONE
    INTEGER, PARAMETER :: nthreads(4)=[1,2,4,8]
    INTEGER :: arg_count,i
    CHARACTER(LEN=64) :: arg

    CALL omp_set_dynamic(.FALSE.)
    CALL run_correctness

    arg_count=COMMAND_ARGUMENT_COUNT()
    DO i=1,arg_count
        CALL GET_COMMAND_ARGUMENT(i,arg)
        SELECT CASE (TRIM(arg))
        CASE ('--benchmark')
            CALL run_benchmarks(500,500)
        CASE ('--benchmark-large')
            CALL run_benchmarks(2000,3000)
        CASE DEFAULT
            WRITE(*,'(A,A)') 'Unknown option: ',TRIM(arg)
            ERROR STOP 1
        END SELECT
    END DO

CONTAINS

    !> @brief Compare exact sorted outputs over geometry, topology and thread-count cases.
    SUBROUTINE run_correctness
        INTEGER, PARAMETER :: sizes(2,9)=RESHAPE([ &
            1,1, 1,7, 7,1, 2,2, 3,4, 7,13, 19,9, 64,32, 129,81 ],[2,9])
        INTEGER :: s,kind,t,nx,ny,checked
        INTEGER, ALLOCATABLE :: old(:),external(:),expected(:),actual(:)
        LOGICAL, ALLOCATABLE :: mask(:,:)

        checked=0
        DO s=1,SIZE(sizes,2)
            nx=sizes(1,s)
            ny=sizes(2,s)
            DO kind=1,15
                CALL make_case(nx,ny,kind,old,external,mask)
                CALL compact_reference(mask,expected)
                DO t=1,SIZE(nthreads)
                    CALL omp_set_num_threads(nthreads(t))
                    CALL compact_global_parallel(mask,actual)
                    CALL compare_results(expected,actual,mask,nx,ny,kind,nthreads(t))
                    checked=checked+1
                    CALL compact_sparse_parallel(old,external,mask,actual)
                    CALL compare_results(expected,actual,mask,nx,ny,kind,nthreads(t))
                    checked=checked+1
                END DO
            END DO
        END DO
        WRITE(*,'(A,I0,A)') 'PASS: ',checked,' exact differential compaction cases'
        WRITE(*,'(2A)') '  Global/sparse OpenMP: empty/full/edge/strip/islands/', &
            'external/exclusions, 1/2/4/8 threads.'
    END SUBROUTINE run_correctness

    !> @brief Generate an old workset, external seeds, and an independent target mask.
    !> @details Uses direct grid-cell halo expansion, independent of the tested merge.
    !> @param[in] nx Number of x-cells.
    !> @param[in] ny Number of y-cells.
    !> @param[in] kind Geometry and activation case selector.
    !> @param[out] old Allocated sorted previous workset.
    !> @param[out] external Allocated sorted newly activated seeds.
    !> @param[out] mask Allocated reference mask after dilation and optional exclusions.
    SUBROUTINE make_case(nx,ny,kind,old,external,mask)
        INTEGER, INTENT(IN) :: nx,ny,kind
        INTEGER, ALLOCATABLE, INTENT(OUT) :: old(:),external(:)
        LOGICAL, ALLOCATABLE, INTENT(OUT) :: mask(:,:)
        LOGICAL, ALLOCATABLE :: old_mask(:,:),event_mask(:,:)
        INTEGER :: j,k,p,cx,cy
        LOGICAL :: selected

        ALLOCATE(old_mask(nx,ny),event_mask(nx,ny),mask(nx,ny))
        old_mask=.FALSE.
        event_mask=.FALSE.
        mask=.FALSE.
        cx=(nx+1)/2
        cy=(ny+1)/2

        DO k=1,ny
            DO j=1,nx
                p=j+(k-1)*nx
                selected=.FALSE.
                SELECT CASE (kind)
                CASE (1,12)
                    selected=.FALSE.
                CASE (13)
                    selected=.TRUE.
                CASE (14)
                    selected=MOD(p,97)<13
                CASE (15)
                    selected=MOD(p,MAX(1,nx*ny/32))==0
                CASE (2)
                    selected=j==cx .AND. k==cy
                CASE (3)
                    selected=(j==1 .OR. j==nx) .AND. (k==1 .OR. k==ny)
                CASE (4)
                    selected=k==cy
                CASE (5)
                    selected=ABS(j-cx)<=MAX(1,nx/6) .AND. ABS(k-cy)<=MAX(1,ny/6)
                CASE (6)
                    selected=MOD(p,13)==0 .OR. MOD(p,19)==0
                CASE (7)
                    selected=.TRUE.
                CASE (8)
                    selected=MOD(p,1000)==0
                CASE (9)
                    selected=MOD(p,11)==0
                CASE (10)
                    selected=j==MAX(1,MIN(nx,k))
                CASE (11)
                    selected=j==1 .OR. j==nx .OR. k==1 .OR. k==ny
                END SELECT
                old_mask(j,k)=selected
            END DO
        END DO

        IF (kind==3 .OR. kind==5 .OR. kind==8 .OR. kind==12) THEN
            event_mask(1,ny)=.TRUE.
            event_mask(nx,1)=.TRUE.
            event_mask(cx,cy)=.TRUE.
        END IF
        IF (kind==14) event_mask(MIN(nx,2),MAX(1,ny-1))=.TRUE.
        IF (kind==11) THEN
            event_mask(1,:)=.TRUE.
            event_mask(nx,:)=.TRUE.
            event_mask(:,1)=.TRUE.
            event_mask(:,ny)=.TRUE.
        END IF

        CALL compact_reference(old_mask,old)
        CALL compact_reference(event_mask,external)

        DO k=1,ny
            DO j=1,nx
                p=j+(k-1)*nx
                IF (old_mask(j,k) .AND. kind/=13 .AND. (kind==7 .OR. MOD(p+kind,11)/=0)) &
                    CALL activate_halo(mask,j,k)
                IF (event_mask(j,k)) CALL activate_halo(mask,j,k)
            END DO
        END DO
        ! Simulate source-cell exclusions applied after the halo expansion.
        IF (kind==10) mask(cx,:)=.FALSE.
    END SUBROUTINE make_case

    !> @brief Activate a cell and its four orthogonal neighbours in a logical mask.
    !> @param[in,out] mask Target logical mask, indexed as (x,y).
    !> @param[in] j One-based x-cell index.
    !> @param[in] k One-based y-cell index.
    SUBROUTINE activate_halo(mask,j,k)
        LOGICAL, INTENT(INOUT) :: mask(:,:)
        INTEGER, INTENT(IN) :: j,k
        INTEGER :: nx,ny

        nx=SIZE(mask,1)
        ny=SIZE(mask,2)
        mask(j,k)=.TRUE.
        IF (j>1) mask(j-1,k)=.TRUE.
        IF (j<nx) mask(j+1,k)=.TRUE.
        IF (k>1) mask(j,k-1)=.TRUE.
        IF (k<ny) mask(j,k+1)=.TRUE.
    END SUBROUTINE activate_halo

    !> @brief Assert exact equality, strict ordering, validity and mask reconstruction.
    !> @param[in] expected Serial reference indices.
    !> @param[in] actual Parallel compaction indices.
    !> @param[in] mask Ground-truth active cell membership.
    !> @param[in] nx Number of x-cells.
    !> @param[in] ny Number of y-cells.
    !> @param[in] kind Test geometry identifier, used for diagnostic errors.
    !> @param[in] nthreads OpenMP thread count, used for diagnostic errors.
    SUBROUTINE compare_results(expected,actual,mask,nx,ny,kind,nthreads)
        INTEGER, INTENT(IN) :: expected(:),actual(:),nx,ny,kind,nthreads
        LOGICAL, INTENT(IN) :: mask(:,:)
        LOGICAL, ALLOCATABLE :: reconstructed(:,:)
        INTEGER :: l,j,k

        IF (SIZE(expected)/=SIZE(actual)) THEN
            WRITE(*,'(A,5(I0,1X))') 'FAIL: count (nx ny case threads actual): ', &
                nx,ny,kind,nthreads,SIZE(actual)
            ERROR STOP 1
        END IF
        IF (ANY(actual/=expected)) THEN
            WRITE(*,'(A,4(I0,1X))') 'FAIL: different indices (nx ny case threads): ',nx,ny,kind,nthreads
            ERROR STOP 1
        END IF
        ALLOCATE(reconstructed(nx,ny))
        reconstructed=.FALSE.
        DO l=2,SIZE(actual)
            IF (actual(l)<=actual(l-1)) ERROR STOP 'Duplicate or unsorted index'
        END DO
        DO l=1,SIZE(actual)
            IF (actual(l)<1 .OR. actual(l)>nx*ny) ERROR STOP 'Out-of-range index'
            j=MOD(actual(l)-1,nx)+1
            k=(actual(l)-1)/nx+1
            reconstructed(j,k)=.TRUE.
        END DO
        IF (.NOT.ALL(reconstructed .EQV. mask)) ERROR STOP 'List/mask inconsistency'
    END SUBROUTINE compare_results

    !> @brief Compare serial, global OpenMP, and sparse OpenMP compaction on two geometries.
    !> @details All three timed routines allocate their output list. Tests and mask setup are
    !> excluded from timings; each median is over seven repetitions after a warmup.
    !> @param[in] nx Number of x-cells used in the performance experiment.
    !> @param[in] ny Number of y-cells used in the performance experiment.
    SUBROUTINE run_benchmarks(nx,ny)
        INTEGER, INTENT(IN) :: nx,ny
        REAL(KIND=8), PARAMETER :: fractions(4)=[0.001D0,0.01D0,0.1D0,1.0D0]
        INTEGER, PARAMETER :: nrep=7
        INTEGER :: t,f,r,j,k,cx,cy,half_x,half_y,kind
        REAL(KIND=8) :: t0,serial_ms,global_ms,sparse_ms
        REAL(KIND=8) :: serial_times(nrep),global_times(nrep),sparse_times(nrep)
        INTEGER, ALLOCATABLE :: old(:),external(:),expected(:),actual(:)
        LOGICAL, ALLOCATABLE :: old_mask(:,:),mask(:,:)

        ALLOCATE(old_mask(nx,ny),mask(nx,ny),external(0))
        cx=(nx+1)/2
        cy=(ny+1)/2
        WRITE(*,'(A,I0,A,I0)') 'BENCHMARK grid: ',nx,' x ',ny
        WRITE(*,'(2A)') ' Shape    Fraction Threads    Old cells   New cells', &
            '  Serial ms  Global ms  Sparse ms  Global/Sparse'

        DO kind=1,2
            DO f=1,SIZE(fractions)
                old_mask=.FALSE.
                half_x=MAX(0,INT(0.5D0*DBLE(nx)*SQRT(fractions(f))))
                half_y=MAX(0,INT(0.5D0*DBLE(ny)*SQRT(fractions(f))))
                SELECT CASE (kind)
                CASE (1)
                    old_mask(MAX(1,cx-half_x):MIN(nx,cx+half_x), &
                        MAX(1,cy-half_y):MIN(ny,cy+half_y))=.TRUE.
                CASE (2)
                    ! Long narrow corridor: occupy a prescribed fraction of rows.
                    half_y=MAX(0,INT(0.5D0*DBLE(ny)*fractions(f)))
                    old_mask(:,MAX(1,cy-half_y):MIN(ny,cy+half_y))=.TRUE.
                END SELECT

                CALL compact_reference(old_mask,old)
                mask=.FALSE.
                DO k=1,ny
                    DO j=1,nx
                        IF (.NOT.old_mask(j,k)) CYCLE
                        IF (MOD(j+nx*(k-1),17)==0) CYCLE
                        CALL activate_halo(mask,j,k)
                    END DO
                END DO
                CALL omp_set_num_threads(1)
                CALL compact_reference(mask,expected)
                CALL compact_reference(mask,actual)
                CALL compare_results(expected,actual,mask,nx,ny,kind,1)
                DO r=1,nrep
                    t0=omp_get_wtime()
                    CALL compact_reference(mask,actual)
                    serial_times(r)=1000.0D0*(omp_get_wtime()-t0)
                END DO
                CALL sort_times(serial_times)
                serial_ms=serial_times((nrep+1)/2)

                DO t=1,SIZE(nthreads)
                    CALL omp_set_num_threads(nthreads(t))
                    CALL compact_global_parallel(mask,actual)
                    CALL compare_results(expected,actual,mask,nx,ny,kind,nthreads(t))
                    CALL compact_sparse_parallel(old,external,mask,actual)
                    CALL compare_results(expected,actual,mask,nx,ny,kind,nthreads(t))
                    DO r=1,nrep
                        t0=omp_get_wtime()
                        CALL compact_global_parallel(mask,actual)
                        global_times(r)=1000.0D0*(omp_get_wtime()-t0)
                        t0=omp_get_wtime()
                        CALL compact_sparse_parallel(old,external,mask,actual)
                        sparse_times(r)=1000.0D0*(omp_get_wtime()-t0)
                    END DO
                    CALL sort_times(global_times)
                    CALL sort_times(sparse_times)
                    global_ms=global_times((nrep+1)/2)
                    sparse_ms=sparse_times((nrep+1)/2)
                    WRITE(*,'(A8,1X,F8.4,1X,I7,1X,I12,1X,I11,3(1X,F10.4),1X,F13.3)') &
                        MERGE('compact ','corridor',kind==1),fractions(f),nthreads(t), &
                        SIZE(old),SIZE(expected),serial_ms,global_ms,sparse_ms, &
                        global_ms/MAX(sparse_ms,TINY(1.0D0))
                END DO
                DEALLOCATE(old,expected,actual)
            END DO
        END DO
    END SUBROUTINE run_benchmarks

    !> @brief Sort elapsed durations to allow selection of their robust median.
    !> @param[in,out] values Elapsed durations in milliseconds.
    SUBROUTINE sort_times(values)
        REAL(KIND=8), INTENT(INOUT) :: values(:)
        REAL(KIND=8) :: tmp
        INTEGER :: i,j

        DO i=1,SIZE(values)-1
            DO j=i+1,SIZE(values)
                IF (values(i)<=values(j)) CYCLE
                tmp=values(i)
                values(i)=values(j)
                values(j)=tmp
            END DO
        END DO
    END SUBROUTINE sort_times

END PROGRAM test_domain_compaction
