!> @file rheology_lahar_viscosity.inc
!> @brief Shared temperature- and solids-dependent lahar dynamic viscosity.
!> @details Uses the configured alpha1 and beta1 coefficients and the water viscosity correlation.
!> Shared REAL/COMPLEX lahar carrier viscosity correlation.
!> The caller declares T, alphas_tot, fluid_viscosity, temperature_celsius, exponent_a, exponent_b and alpha1
!> using the same REAL or COMPLEX kind, with the layer-specific config in scope.
temperature_celsius = T - 273.15_wp

IF (REAL(temperature_celsius,wp) .LT. 20.0_wp) THEN
    exponent_a = 1301.0_wp / (998.333_wp + 8.1855_wp *                                                                 &
        (temperature_celsius - 20.0_wp) + 0.00585_wp *                                                                 &
        (temperature_celsius - 20.0_wp)**2) - 1.30223_wp
    alpha1 = config%alpha1_coeff * 1.0E-3_wp * 10.0_wp**exponent_a
ELSE
    exponent_b = (1.3272_wp * (20.0_wp - temperature_celsius)                                                          &
        - 0.001053_wp * (temperature_celsius - 20.0_wp)**2)                                                            &
        / (temperature_celsius + 105.0_wp)
    alpha1 = config%alpha1_coeff * 1.002E-3_wp * 10.0_wp**exponent_b
END IF

fluid_viscosity = alpha1 * EXP(config%beta1 * alphas_tot)
