!> @file rheology_newmodel_linear.inc
!> @brief Linear momentum drag shared by the real and complex rheology kernels.
!> @details Uses the current mixture density and the caller-supplied drag coefficient.
            rate_x=-drag_coefficient*state%rho_m*state%u
            rate_y=-drag_coefficient*state%rho_m*state%v
