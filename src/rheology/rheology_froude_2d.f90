!********************************************************************************
!> \brief Froude-dependent Coulomb model 9 source terms.
!********************************************************************************

MODULE rheology_froude_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : mu_0, mu_inf, Fr_0

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_froude_source

CONTAINS

  !******************************************************************************
  !> \brief Evaluate model 9's Froude-dependent Coulomb friction.
  !> \details Follows Zhu et al. 2020 (doi:10.1007/s10035-020-01053-7).
  !******************************************************************************

  SUBROUTINE eval_froude_source( flow_depth, velocity_x, velocity_y, speed,  &
       horizontal_speed, mixture_density, reduced_gravity, gravity_coefficient, &
       curvature_term, stochastic, stochastic_value, source_x, source_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: froude_number
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp

    IF ( speed .GT. 0.0_wp .AND. reduced_gravity .GT. 0.0_wp ) THEN
       froude_number = speed / SQRT(reduced_gravity * flow_depth)

       IF ( stochastic ) THEN
          IF ( flow_depth .GT. 0.01_wp ) froude_number = froude_number + stochastic_value
          IF ( froude_number .LT. 0.0_wp ) froude_number = 0.0_wp
       END IF

       friction_coefficient = mu_inf + (mu_0 - mu_inf) * EXP(-froude_number/Fr_0)
       friction_stress = mixture_density * ( friction_coefficient * flow_depth * &
            gravity_coefficient * ( reduced_gravity + curvature_term ) )
       source_x = source_x - friction_stress * velocity_x / horizontal_speed
       source_y = source_y - friction_stress * velocity_y / horizontal_speed
    END IF

  END SUBROUTINE eval_froude_source

END MODULE rheology_froude_2d
