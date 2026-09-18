!********************************************************************************
!> \brief Random-number support for the stochastic model
!
!> The intrinsic generator is seeded once during stochastic-workspace
!> initialization.  A non-negative stochastic_seed gives reproducible runs;
!> the default negative value requests processor-dependent initialization.
!********************************************************************************
MODULE stochastic_random_2d

  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  INTEGER, PUBLIC :: stochastic_seed = -1

  PUBLIC :: initialize_stochastic_rng
  PUBLIC :: gaussian_noise

CONTAINS

  SUBROUTINE initialize_stochastic_rng

    INTEGER :: i
    INTEGER :: seed_size
    INTEGER, ALLOCATABLE :: seed_state(:)
    INTEGER(int64) :: seed_modulus

    IF (stochastic_seed < 0) THEN
       CALL random_seed
       RETURN
    END IF

    CALL random_seed(SIZE=seed_size)
    ALLOCATE(seed_state(seed_size))

    seed_modulus = INT(HUGE(0), int64)
    DO i = 1, seed_size
       seed_state(i) = INT(MODULO(INT(stochastic_seed, int64) +             &
            104729_int64 * INT(i, int64), seed_modulus))
    END DO

    CALL random_seed(PUT=seed_state)
    DEALLOCATE(seed_state)

  END SUBROUTINE initialize_stochastic_rng

  FUNCTION gaussian_noise(noise_size) RESULT(noise)

    INTEGER, INTENT(IN) :: noise_size
    REAL(wp) :: noise(noise_size)

    REAL(wp) :: uniform_samples(noise_size, 2)
    REAL(wp), PARAMETER :: two_pi = 2.0_wp * ACOS(-1.0_wp)

    CALL random_number(uniform_samples)

    ! Box-Muller transform. RANDOM_NUMBER returns values in [0,1), so
    ! 1-uniform_samples(:,1) is strictly positive and LOG is well-defined.
    noise = SQRT(-2.0_wp * LOG(1.0_wp - uniform_samples(:, 1))) *           &
         COS(two_pi * uniform_samples(:, 2))

  END FUNCTION gaussian_noise

END MODULE stochastic_random_2d
