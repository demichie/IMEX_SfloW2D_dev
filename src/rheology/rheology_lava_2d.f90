!********************************************************************************
!> \brief Temperature-dependent lava rheology model terms.
!********************************************************************************

MODULE rheology_lava_2d

  USE parameters_2d, ONLY : wp
  USE constitutive_parameters_2d, ONLY : convective_term_coeff,                &
       radiative_term_coeff, T_env
  USE rheology_parameters_2d, ONLY : T_ref

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nu_ref, visc_par, tau0
  PUBLIC :: eval_lava_implicit_source, eval_lava_semi_impl_source

  REAL(wp) :: nu_ref
  REAL(wp) :: visc_par
  REAL(wp) :: tau0

CONTAINS

  SUBROUTINE eval_lava_implicit_source(flow_depth, temperature, mixture_density, &
       velocity_x, velocity_y, speed, moving, source_x, source_y, source_energy)

    COMPLEX(wp), INTENT(IN) :: flow_depth, temperature, mixture_density
    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y, source_energy

    COMPLEX(wp) :: gamma
    REAL(wp), PARAMETER :: h_threshold = 1.0E-10_wp

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (REAL(flow_depth, wp) .GT. h_threshold) THEN
       gamma = 3.0_wp * nu_ref / flow_depth * EXP(-visc_par * (temperature - T_ref))
    ELSE
       gamma = 3.0_wp * nu_ref / h_threshold * EXP(-visc_par * (temperature - T_ref))
    END IF

    IF (moving) THEN
       source_x = -mixture_density * gamma * velocity_x
       source_y = -mixture_density * gamma * velocity_y
    END IF

    source_energy = -radiative_term_coeff * (temperature**4 - T_env**4)       &
         - convective_term_coeff * (temperature - T_env)

  END SUBROUTINE eval_lava_implicit_source

  SUBROUTINE eval_lava_semi_impl_source(velocity_x, velocity_y, speed, source_x, source_y)

    REAL(wp), INTENT(IN) :: velocity_x, velocity_y, speed
    REAL(wp), INTENT(OUT) :: source_x, source_y

    source_x = 0.0_wp
    source_y = 0.0_wp
    IF (speed .LE. 0.0_wp) RETURN

    source_x = -tau0 * velocity_x / speed
    source_y = -tau0 * velocity_y / speed

  END SUBROUTINE eval_lava_semi_impl_source

END MODULE rheology_lava_2d
