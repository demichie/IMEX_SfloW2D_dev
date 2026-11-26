!> @file thermal_cooling_rate.inc
!> @brief Shared convective and radiative cooling rate.
!> @details Evaluates the Stefan-Boltzmann-like fourth-power term and linear ambient-temperature exchange.
        rate=-radiative*(temperature**4-ambient**4)-convective*(temperature-ambient)
