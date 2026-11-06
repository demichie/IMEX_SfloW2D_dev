!> @file reduced_gravity.inc
!> @brief Compute reduced gravity relative to ambient air density.
    red_grav=(rho_m-rho_a_amb)*inv_rho_m*grav
