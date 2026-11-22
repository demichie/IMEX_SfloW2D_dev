! Standalone test kernels: no dependency on the IMEX-SfloW2D solver modules.
MODULE domain_compaction_test_kernels
    USE omp_lib, ONLY: omp_get_max_threads
    USE iso_fortran_env, ONLY: int64
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: compact_reference, compact_sparse_parallel

CONTAINS

    ! Reference: globally scan the logical mask in the existing (k,j) order.
    SUBROUTINE compact_reference(mask, indices)
        LOGICAL, INTENT(IN) :: mask(:,:)
        INTEGER, ALLOCATABLE, INTENT(OUT) :: indices(:)
        INTEGER :: j,k,n,nx,ny

        nx=SIZE(mask,1)
        ny=SIZE(mask,2)
        ALLOCATE(indices(COUNT(mask)))
        n=0
        DO k=1,ny
            DO j=1,nx
                IF (.NOT.mask(j,k)) CYCLE
                n=n+1
                indices(n)=j+(k-1)*nx
            END DO
        END DO
    END SUBROUTINE compact_reference

    ! Sparse, deterministic compaction. Neither a global scan nor per-thread output lists.
    ! old_indices and external_seeds are sorted and contain no duplicates internally.
    ! All true mask entries must be inside the one-cell halo of their union.
    SUBROUTINE compact_sparse_parallel(old_indices,external_seeds,mask,indices)
        INTEGER, INTENT(IN) :: old_indices(:),external_seeds(:)
        LOGICAL, INTENT(IN) :: mask(:,:)
        INTEGER, ALLOCATABLE, INTENT(OUT) :: indices(:)
        INTEGER, ALLOCATABLE :: counts(:),offsets(:)
        INTEGER :: blocks,b,lo,hi,nx,ny,nfound,total

        nx=SIZE(mask,1)
        ny=SIZE(mask,2)
        total=nx*ny
        ! Small destination-index intervals help balance nonuniform sparse regions.
        blocks=MAX(1,32*omp_get_max_threads())
        ALLOCATE(counts(blocks),offsets(blocks))

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lo,hi,nfound) SCHEDULE(STATIC)
        DO b=1,blocks
            CALL block_limits(total,blocks,b,lo,hi)
            CALL merge_interval(old_indices,external_seeds,mask,lo,hi,nfound)
            counts(b)=nfound
        END DO
        !$OMP END PARALLEL DO

        nfound=0
        DO b=1,blocks
            offsets(b)=nfound
            nfound=nfound+counts(b)
        END DO
        ALLOCATE(indices(nfound))

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lo,hi,nfound) SCHEDULE(STATIC)
        DO b=1,blocks
            CALL block_limits(total,blocks,b,lo,hi)
            CALL merge_interval(old_indices,external_seeds,mask,lo,hi,nfound,indices,offsets(b))
            IF (nfound/=counts(b)) ERROR STOP 'Compaction count changed between passes'
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE compact_sparse_parallel

    SUBROUTINE block_limits(total,blocks,b,lo,hi)
        INTEGER, INTENT(IN) :: total,blocks,b
        INTEGER, INTENT(OUT) :: lo,hi

        lo=1+INT(INT(total,int64)*INT(b-1,int64)/INT(blocks,int64))
        hi=INT(INT(total,int64)*INT(b,int64)/INT(blocks,int64))
    END SUBROUTINE block_limits

    ! Merge five shifted views of each input, filtered by a destination-index interval.
    ! The shifted views are virtual: no candidate arrays or duplicates are stored.
    SUBROUTINE merge_interval(old_indices,external_seeds,mask,lo,hi,nfound,indices,offset)
        INTEGER, INTENT(IN) :: old_indices(:),external_seeds(:),lo,hi
        LOGICAL, INTENT(IN) :: mask(:,:)
        INTEGER, INTENT(OUT) :: nfound
        INTEGER, INTENT(INOUT), OPTIONAL :: indices(:)
        INTEGER, INTENT(IN), OPTIONAL :: offset
        INTEGER :: nx,ny,d,s,p,j,k,previous,delta(5),cursor(5,2),next_p(5,2)

        nx=SIZE(mask,1)
        ny=SIZE(mask,2)
        delta=[0,-1,1,-nx,nx]
        nfound=0
        IF (hi<lo) RETURN

        DO s=1,2
            DO d=1,5
                IF (s==1) THEN
                    cursor(d,s)=lower_bound(old_indices,lo-delta(d))
                    CALL next_candidate(old_indices,delta(d),d,nx,ny,lo,hi,cursor(d,s),next_p(d,s))
                ELSE
                    cursor(d,s)=lower_bound(external_seeds,lo-delta(d))
                    CALL next_candidate(external_seeds,delta(d),d,nx,ny,lo,hi,cursor(d,s),next_p(d,s))
                END IF
            END DO
        END DO

        previous=0
        DO
            p=MINVAL(next_p)
            IF (p==HUGE(p)) EXIT

            IF (p/=previous) THEN
                j=MOD(p-1,nx)+1
                k=(p-1)/nx+1
                IF (mask(j,k)) THEN
                    nfound=nfound+1
                    IF (PRESENT(indices)) THEN
                        IF (.NOT.PRESENT(offset)) ERROR STOP 'Missing output offset'
                        indices(offset+nfound)=p
                    END IF
                END IF
                previous=p
            END IF

            DO s=1,2
                DO d=1,5
                    IF (next_p(d,s)/=p) CYCLE
                    cursor(d,s)=cursor(d,s)+1
                    IF (s==1) THEN
                        CALL next_candidate(old_indices,delta(d),d,nx,ny,lo,hi,cursor(d,s),next_p(d,s))
                    ELSE
                        CALL next_candidate(external_seeds,delta(d),d,nx,ny,lo,hi,cursor(d,s),next_p(d,s))
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE merge_interval

    ! First position whose source index is >= value, or SIZE(indices)+1.
    INTEGER FUNCTION lower_bound(indices,value) RESULT(pos)
        INTEGER, INTENT(IN) :: indices(:),value
        INTEGER :: left,right,middle

        left=1
        right=SIZE(indices)+1
        DO WHILE (left<right)
            middle=left+(right-left)/2
            IF (indices(middle)<value) THEN
                left=middle+1
            ELSE
                right=middle
            END IF
        END DO
        pos=left
    END FUNCTION lower_bound

    ! Skip out-of-domain or row-wrapping offsets; each stream remains sorted.
    SUBROUTINE next_candidate(base,delta,direction,nx,ny,lo,hi,cursor,candidate)
        INTEGER, INTENT(IN) :: base(:),delta,direction,nx,ny,lo,hi
        INTEGER, INTENT(INOUT) :: cursor
        INTEGER, INTENT(OUT) :: candidate
        INTEGER :: source,p
        LOGICAL :: valid

        candidate=HUGE(candidate)
        DO WHILE (cursor<=SIZE(base))
            source=base(cursor)
            p=source+delta
            IF (p>hi) RETURN
            valid=p>=lo .AND. p>=1 .AND. p<=nx*ny
            IF (valid) THEN
                SELECT CASE (direction)
                CASE (2)
                    valid=MOD(source-1,nx)>0
                CASE (3)
                    valid=MOD(source-1,nx)<nx-1
                END SELECT
            END IF
            IF (valid) THEN
                candidate=p
                RETURN
            END IF
            cursor=cursor+1
        END DO
    END SUBROUTINE next_candidate

END MODULE domain_compaction_test_kernels
