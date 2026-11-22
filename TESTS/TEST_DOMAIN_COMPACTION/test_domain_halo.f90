!> @brief Compare OpenMP ping-pong domain halos with the serial production stencil.
PROGRAM test_domain_halo
    USE omp_lib, ONLY : omp_set_dynamic,omp_set_num_threads
    IMPLICIT NONE
    INTEGER, PARAMETER :: nx_cases(8)=[1,1,2,3,5,17,73,231]
    INTEGER, PARAMETER :: ny_cases(8)=[1,7,5,3,6,23,47,137]
    INTEGER, PARAMETER :: threads(4)=[1,2,4,8]
    LOGICAL, ALLOCATABLE :: seed(:,:),excluded(:,:),reference(:,:),mask(:,:),temp(:,:)
    INTEGER, ALLOCATABLE :: ref_count(:),row_count(:),offset(:),list_ref(:),list_got(:)
    INTEGER :: c,scenario,ti,nx,ny,j,k,expected_size,p,checks

    CALL omp_set_dynamic(.FALSE.)
    checks=0
    DO c=1,SIZE(nx_cases)
        nx=nx_cases(c)
        ny=ny_cases(c)
        ALLOCATE(seed(nx,ny),excluded(nx,ny),reference(nx,ny),mask(nx,ny),temp(nx,ny))
        ALLOCATE(ref_count(ny),row_count(ny),offset(ny+1),list_ref(nx*ny),list_got(nx*ny))
        DO scenario=1,10
            DO k=1,ny
                DO j=1,nx
                    SELECT CASE(scenario)
                    CASE(1)
                        seed(j,k)=.FALSE.
                        excluded(j,k)=.FALSE.
                    CASE(2)
                        seed(j,k)=.TRUE.
                        excluded(j,k)=.FALSE.
                    CASE(3)
                        seed(j,k)=j==1 .OR. j==nx .OR. k==1 .OR. k==ny
                        excluded(j,k)=.FALSE.
                    CASE(4)
                        seed(j,k)=j==(nx+1)/2 .AND. k==(ny+1)/2
                        excluded(j,k)=.FALSE.
                    CASE(5)
                        seed(j,k)=MOD(17*j+31*k,11)==0
                        excluded(j,k)=j==1 .OR. j==nx .OR. k==1 .OR. k==ny
                    CASE(6)
                        seed(j,k)=MOD(17*j+31*k,7)==0
                        excluded(j,k)=MOD(13*j+19*k,5)==0
                    CASE(7)
                        seed(j,k)=MOD(j+k,2)==0
                        excluded(j,k)=.TRUE.
                    CASE(8)
                        seed(j,k)=j==2 .OR. k==2
                        excluded(j,k)=MOD(23*j+7*k,13)==0
                    CASE(9)
                        seed(j,k)=MOD(83*j+29*k,37)==0
                        excluded(j,k)=MOD(j+k,7)==0
                    CASE DEFAULT
                        seed(j,k)=MOD(19*j+43*k,101)<13
                        excluded(j,k)=MOD(7*j+11*k,17)==0
                    END SELECT
                END DO
            END DO
            CALL reference_halo(seed,excluded,reference,ref_count)
            expected_size=COUNT(reference)
            p=0
            DO k=1,ny
                DO j=1,nx
                    IF (.NOT.reference(j,k)) CYCLE
                    p=p+1
                    list_ref(p)=j+(k-1)*nx
                END DO
            END DO
            DO ti=1,SIZE(threads)
                CALL omp_set_num_threads(threads(ti))
                mask=seed
                temp=.NOT.seed
                CALL fused_halo(mask,temp,excluded,row_count)
                IF (ANY(mask.NEQV.reference)) ERROR STOP 'OpenMP halo membership differs'
                IF (ANY(row_count/=ref_count)) ERROR STOP 'OpenMP row counts differ'
                offset(1)=0
                DO k=1,ny
                    offset(k+1)=offset(k)+row_count(k)
                END DO
                IF (offset(ny+1)/=expected_size) ERROR STOP 'OpenMP prefix sum differs'
                CALL emit_parallel(mask,offset,list_got)
                IF (ANY(list_got(1:expected_size)/=list_ref(1:expected_size))) &
                    ERROR STOP 'OpenMP ordered cell list differs'
                checks=checks+1
            END DO
        END DO
        DEALLOCATE(seed,excluded,reference,mask,temp,ref_count,row_count,offset,list_ref,list_got)
    END DO
    WRITE(*,'(A,I0)') 'PASS: exact halo, exclusion and ordered compaction comparisons: ',checks
CONTAINS
    !> @brief Evaluate the baseline two-copy, two-pass halo serially.
    !> @param[in] seeds Initial wetness and source membership.
    !> @param[in] exclusion Post-halo radial exclusion mask.
    !> @param[out] result Fully dilated, excluded solve mask.
    !> @param[out] counts Final active count per row.
    SUBROUTINE reference_halo(seeds,exclusion,result,counts)
        LOGICAL, INTENT(IN) :: seeds(:,:),exclusion(:,:)
        LOGICAL, INTENT(OUT) :: result(:,:)
        INTEGER, INTENT(OUT) :: counts(:)
        LOGICAL :: snapshot(SIZE(seeds,1),SIZE(seeds,2))
        INTEGER :: i,j,k,nx,ny
        nx=SIZE(seeds,1)
        ny=SIZE(seeds,2)
        result=seeds
        DO i=1,2
            snapshot=result
            DO k=2,ny-1
                DO j=2,nx-1
                    result(j,k)=result(j,k) .OR. snapshot(j-1,k) .OR. snapshot(j+1,k) .OR. &
                        snapshot(j,k-1) .OR. snapshot(j,k+1)
                END DO
            END DO
        END DO
        WHERE (exclusion) result=.FALSE.
        DO k=1,ny
            counts(k)=COUNT(result(:,k))
        END DO
    END SUBROUTINE reference_halo

    !> @brief Compute two OpenMP halo rings and the final row counts without full mask copies.
    !> @param[in,out] current Seed mask, replaced by final halo.
    !> @param[in,out] scratch Existing second mask used for the first halo.
    !> @param[in] exclusion Post-halo radial source exclusion.
    !> @param[out] counts Final row cardinalities.
    SUBROUTINE fused_halo(current,scratch,exclusion,counts)
        LOGICAL, INTENT(INOUT) :: current(:,:),scratch(:,:)
        LOGICAL, INTENT(IN) :: exclusion(:,:)
        INTEGER, INTENT(OUT) :: counts(:)
        LOGICAL :: active
        INTEGER :: j,k,n,nx,ny,km,kp
        nx=SIZE(current,1)
        ny=SIZE(current,2)
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,n,active,km,kp)
        !$OMP SINGLE
        scratch(1,:)=current(1,:)
        scratch(nx,:)=current(nx,:)
        scratch(:,1)=current(:,1)
        scratch(:,ny)=current(:,ny)
        !$OMP END SINGLE
        !$OMP DO SCHEDULE(STATIC)
        DO k=2,ny-1
            DO j=2,nx-1
                scratch(j,k)=current(j,k) .OR. current(j-1,k) .OR. current(j+1,k) .OR. &
                    current(j,k-1) .OR. current(j,k+1)
            END DO
        END DO
        !$OMP END DO
        !$OMP DO SCHEDULE(STATIC)
        DO k=1,ny
            n=0
            IF (k==1 .OR. k==ny .OR. nx<=2) THEN
                DO j=1,nx
                    active=scratch(j,k)
                    IF (exclusion(j,k)) active=.FALSE.
                    current(j,k)=active
                    IF (active) n=n+1
                END DO
            ELSE
                km=k-1
                kp=k+1
                active=scratch(1,k)
                IF (exclusion(1,k)) active=.FALSE.
                current(1,k)=active
                IF (active) n=n+1
                DO j=2,nx-1
                    active=scratch(j,k) .OR. scratch(j-1,k) .OR. scratch(j+1,k) .OR. &
                        scratch(j,km) .OR. scratch(j,kp)
                    IF (exclusion(j,k)) active=.FALSE.
                    current(j,k)=active
                    IF (active) n=n+1
                END DO
                active=scratch(nx,k)
                IF (exclusion(nx,k)) active=.FALSE.
                current(nx,k)=active
                IF (active) n=n+1
            END IF
            counts(k)=n
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE fused_halo

    !> @brief Emit ordered linear cell indices using disjoint OpenMP row intervals.
    !> @param[in] membership Final shared domain mask.
    !> @param[in] offset Exclusive prefix sum of final row counts.
    !> @param[out] indices Ordered cell indices.
    SUBROUTINE emit_parallel(membership,offset,indices)
        LOGICAL, INTENT(IN) :: membership(:,:)
        INTEGER, INTENT(IN) :: offset(:)
        INTEGER, INTENT(OUT) :: indices(:)
        INTEGER :: j,k,n,nx
        nx=SIZE(membership,1)
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,n) SCHEDULE(STATIC)
        DO k=1,SIZE(membership,2)
            n=offset(k)
            DO j=1,nx
                IF (.NOT.membership(j,k)) CYCLE
                n=n+1
                indices(n)=j+(k-1)*nx
            END DO
        END DO
        !$OMP END PARALLEL DO
    END SUBROUTINE emit_parallel
END PROGRAM test_domain_halo
