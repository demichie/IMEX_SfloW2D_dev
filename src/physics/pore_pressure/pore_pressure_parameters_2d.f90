!> @brief Canonical layer-local Darcy configuration and independently owned derived cache.
!> @details NAMELIST configuration contains scalar data only. Pascal coefficients
!>          are derived once during input parsing and passed read-only to cell kernels.
MODULE pore_pressure_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: pore_pressure_config_type, pore_pressure_l1, pore_pressure_l2
    PUBLIC :: pore_pressure_cache_type, pore_pressure_cache_l1, pore_pressure_cache_l2
    PUBLIC :: precompute_pascal_coefficient

    !> @brief Darcy material closure and source scheduling for one layer.
    TYPE :: pore_pressure_config_type
        REAL(wp) :: hydraulic_permeability=0.0_wp  !< Fixed hydraulic permeability [m2].
        REAL(wp) :: alpha_trans=0.0_wp             !< Solid volume fraction for inhibition onset.
        INTEGER :: N_inh=2                        !< Smoothstep polynomial half order.
        CHARACTER(LEN=8) :: f_inhibit_mode='OFF'  !< OFF, STEP, STATIC, DYNAMIC.
        LOGICAL :: dynamic_permeability_flag=.FALSE.
        CHARACTER(LEN=16) :: pore_pressure_method=''
        INTEGER :: pore_pressure_group=3
    END TYPE pore_pressure_config_type

    !> @brief Derived Pascal coefficients, never a second copy of physical parameters.
    TYPE :: pore_pressure_cache_type
        REAL(wp), ALLOCATABLE :: pascal_coeff(:)  !< Binomial products indexed 0:N_inh.
    END TYPE pore_pressure_cache_type

    ! Canonical layer objects: L2 remains inactive in the current solver.
    TYPE(pore_pressure_config_type), TARGET :: pore_pressure_l1
    TYPE(pore_pressure_config_type), TARGET :: pore_pressure_l2
    TYPE(pore_pressure_cache_type), TARGET :: pore_pressure_cache_l1
    TYPE(pore_pressure_cache_type), TARGET :: pore_pressure_cache_l2
CONTAINS
    !> @brief Precompute the historical Pascal products without changing their evaluation order.
    !> @param[in] config Canonical layer Darcy configuration.
    !> @param[in,out] cache Derived coefficient storage owned by the same layer.
    SUBROUTINE precompute_pascal_coefficient(config,cache)
        TYPE(pore_pressure_config_type), INTENT(IN) :: config
        TYPE(pore_pressure_cache_type), INTENT(INOUT) :: cache
        INTEGER :: N_order, n, k
        REAL(wp) :: pascal1, pascal2

        N_order=config%N_inh
        IF (N_order<0) THEN
            WRITE(*,*) 'ERROR: N must be >= 0. N =',N_order
            STOP
        END IF

        IF (ALLOCATED(cache%pascal_coeff)) THEN
            IF (LBOUND(cache%pascal_coeff,1)/=0 .OR. UBOUND(cache%pascal_coeff,1)/=N_order) THEN
                DEALLOCATE(cache%pascal_coeff)
                ALLOCATE(cache%pascal_coeff(0:N_order))
            END IF
        ELSE
            ALLOCATE(cache%pascal_coeff(0:N_order))
        END IF

        DO n=0,N_order
            pascal1=1.0_wp
            IF (n>0) THEN
                DO k=0,n-1
                    pascal1=pascal1*(-REAL(N_order,kind=wp)-1.0_wp-REAL(k,kind=wp))/REAL(k+1,kind=wp)
                END DO
            END IF

            pascal2=1.0_wp
            IF (N_order-n>0) THEN
                DO k=0,N_order-n-1
                    pascal2=pascal2*(2.0_wp*REAL(N_order,kind=wp)+1.0_wp-REAL(k,kind=wp))/REAL(k+1,kind=wp)
                END DO
            END IF

            cache%pascal_coeff(n)=pascal1*pascal2
        END DO

    END SUBROUTINE precompute_pascal_coefficient
END MODULE pore_pressure_parameters_2d
