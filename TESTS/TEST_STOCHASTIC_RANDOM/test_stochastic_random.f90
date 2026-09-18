PROGRAM test_stochastic_random

  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE parameters_2d, ONLY : wp
  USE stochastic_random_2d, ONLY : stochastic_seed,                         &
       initialize_stochastic_rng, gaussian_noise

  IMPLICIT NONE

  INTEGER, PARAMETER :: sample_size = 64
  INTEGER :: seed_size
  INTEGER, ALLOCATABLE :: saved_state(:)
  REAL(wp) :: first_sample(sample_size)
  REAL(wp) :: second_sample(sample_size)
  REAL(wp) :: repeated_sample(sample_size)
  REAL(wp) :: restored_sample(sample_size)

  stochastic_seed = 24681357
  CALL initialize_stochastic_rng
  first_sample = gaussian_noise(sample_size)

  CALL random_seed(SIZE=seed_size)
  ALLOCATE(saved_state(seed_size))
  CALL random_seed(GET=saved_state)
  second_sample = gaussian_noise(sample_size)

  IF (bitwise_equal(first_sample, second_sample))                           &
       ERROR STOP 'successive samples unexpectedly repeat'

  stochastic_seed = 24681357
  CALL initialize_stochastic_rng
  repeated_sample = gaussian_noise(sample_size)

  IF (.NOT. bitwise_equal(first_sample, repeated_sample))                   &
       ERROR STOP 'fixed seed does not reproduce the first sample'

  CALL random_seed(PUT=saved_state)
  restored_sample = gaussian_noise(sample_size)

  IF (.NOT. bitwise_equal(second_sample, restored_sample))                  &
       ERROR STOP 'restored generator state does not reproduce continuation'

  DEALLOCATE(saved_state)

  WRITE(*,*) 'PASS: stochastic random stream is reproducible and continuous'

CONTAINS

  LOGICAL FUNCTION bitwise_equal(first, second)

    REAL(wp), INTENT(IN) :: first(:)
    REAL(wp), INTENT(IN) :: second(:)
    INTEGER(int64) :: first_bits(SIZE(first))
    INTEGER(int64) :: second_bits(SIZE(second))

    first_bits = TRANSFER(first, first_bits)
    second_bits = TRANSFER(second, second_bits)
    bitwise_equal = ALL(first_bits == second_bits)

  END FUNCTION bitwise_equal

END PROGRAM test_stochastic_random
