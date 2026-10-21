!********************************************************************************
!> \brief Lucas velocity-weakening model 10 source terms.
!********************************************************************************

MODULE rheology_lucas_2d

  USE parameters_2d, ONLY : wp
  USE rheology_parameters_2d, ONLY : mu_0, mu_inf, U_w

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eval_lucas_source

CONTAINS

  !******************************************************************************
  !> \brief Evaluate model 10's velocity-dependent Coulomb friction.
  !> \details Follows Lucas et al. 2014 (doi:10.1038/ncomms4417).
  !******************************************************************************

  SUBROUTINE eval_lucas_source( flow_depth, velocity_x, velocity_y, speed,   &
       horizontal_speed, mixture_density, reduced_gravity, gravity_coefficient, &
       curvature_term, stochastic, stochastic_value, source_x, source_y )

    REAL(wp), INTENT(IN) :: flow_depth, velocity_x, velocity_y, speed
    REAL(wp), INTENT(IN) :: horizontal_speed, mixture_density, reduced_gravity
    REAL(wp), INTENT(IN) :: gravity_coefficient, curvature_term
    LOGICAL, INTENT(IN) :: stochastic
    REAL(wp), INTENT(IN) :: stochastic_value
    REAL(wp), INTENT(OUT) :: source_x, source_y

    REAL(wp) :: friction_speed
    REAL(wp) :: friction_coefficient
    REAL(wp) :: friction_stress

    source_x = 0.0_wp
    source_y = 0.0_wp

    IF ( speed .GT. 0.0_wp ) THEN
       IF ( stochastic ) THEN
          friction_speed = speed + stochastic_value
       ELSE
          friction_speed = speed
       END IF

       IF ( friction_speed .LT. U_w ) THEN
          friction_coefficient = mu_0
       ELSE
          IF ( flow_depth .LT. 0.1_wp ) THEN
             friction_coefficient = mu_0
          ELSE
             friction_coefficient = mu_inf + (mu_0-mu_inf)/(friction_speed/U_w)
          END IF
       END IF

       friction_stress = mixture_density * ( friction_coefficient * flow_depth * &
            gravity_coefficient * ( reduced_gravity + curvature_term ) )
       source_x = source_x - friction_stress * velocity_x / horizontal_speed
       source_y = source_y - friction_stress * velocity_y / horizontal_speed
    END IF

  END SUBROUTINE eval_lucas_source

END MODULE rheology_lucas_2d
