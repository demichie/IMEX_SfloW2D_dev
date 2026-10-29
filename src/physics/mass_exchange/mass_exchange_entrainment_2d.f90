!> @brief Evaluate entrainment process contributions to mass, energy and bed evolution.
!> @details Explicit treatment is supported; other treatment interfaces reject
!> unsupported configurations. Shared resource budgets belong to the coordinator.
MODULE mass_exchange_entrainment_2d
    USE parameters_2d, ONLY : wp, n_vars, slope_correction_flag, mass_exchange_rates_type, reset_mass_exchange_rates
    USE diagnostics_2d, ONLY : fatal_error
    USE state_conversion_2d, ONLY : mixt_var
    USE constitutive_parameters_2d, ONLY : entrainment_flag, rho_a_amb, &
              T_ambient, sp_heat_a
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_mass_exchange_entrainment_rates_explicit
    PUBLIC :: eval_mass_exchange_entrainment_explicit
    PUBLIC :: eval_mass_exchange_entrainment_semi_implicit
    PUBLIC :: eval_mass_exchange_entrainment_implicit
    PUBLIC :: solve_mass_exchange_entrainment_analytic
    PUBLIC :: eval_air_entrainment_rate
    ! Four-entry physical treatment contract; only EXPLICIT is supported.
CONTAINS
    !> @brief Evaluate explicit ambient-air entrainment.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] u Horizontal x velocity of the mobile material [m/s].
    !> @param[in] v Horizontal y velocity of the mobile material [m/s].
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_explicit(qp,u,v,Bx,By,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2),u,v,Bx,By
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        REAL(wp) :: vertical_velocity,speed_squared,Ri,rhom,rhoc,redgrav,cpc,cpm
        CALL reset_mass_exchange_rates(r)
        vertical_velocity=0.0_wp
        IF (slope_correction_flag) vertical_velocity=u*Bx+v*By
        speed_squared=u*u+v*v+vertical_velocity*vertical_velocity
        CALL mixt_var(qp,Ri,rhom,rhoc,redgrav,cpc,cpm)
        ! Convert entrained ambient-air volume to mixture mass at ambient density.
        r%mixture_mass_rate=rho_a_amb*eval_entrainment_volume_kernel(qp(1),Ri,speed_squared)
        r%energy_rate=T_ambient*sp_heat_a*r%mixture_mass_rate
    END SUBROUTINE
    !> @brief Reject unsupported semi-implicit treatment for entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] u Horizontal x velocity of the mobile material [m/s].
    !> @param[in] v Horizontal y velocity of the mobile material [m/s].
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_semi_implicit(qp,u,v,Bx,By,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2),u,v,Bx,By
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_entrainment_semi_implicit
    !> @brief Reject unsupported implicit treatment for ambient-air entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] u Horizontal x velocity of the mobile material [m/s].
    !> @param[in] v Horizontal y velocity of the mobile material [m/s].
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE eval_mass_exchange_entrainment_implicit(qp,u,v,Bx,By,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2),u,v,Bx,By
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE eval_mass_exchange_entrainment_implicit
    !> @brief Reject unsupported analytic treatment for ambient-air entrainment.
    !> @details This entry point raises fatal_error; the rates are reset first.
    !> @param[in] qp Local primitive flow state; component order follows model configuration.
    !> @param[in] u Horizontal x velocity of the mobile material [m/s].
    !> @param[in] v Horizontal y velocity of the mobile material [m/s].
    !> @param[in] Bx Bed elevation gradient in the x direction [dimensionless].
    !> @param[in] By Bed elevation gradient in the y direction [dimensionless].
    !> @param[in,out] r Process rate record; reset or overwritten by the selected treatment.
    SUBROUTINE solve_mass_exchange_entrainment_analytic(qp,u,v,Bx,By,r)
        REAL(wp), INTENT(IN) :: qp(n_vars+2),u,v,Bx,By
        TYPE(mass_exchange_rates_type), INTENT(INOUT) :: r
        CALL reset_mass_exchange_rates(r)
        CALL fatal_error('entrainment: unsupported physical treatment')
    END SUBROUTINE solve_mass_exchange_entrainment_analytic
    !> @brief Switch air entrainment off in dry or non-buoyant conditions.
    !> @param[in] r_h Mobile flow thickness [m].
    !> @param[in] r_Ri Mixture Richardson number [dimensionless].
    !> @param[in] mod_vel2 Squared slope-corrected speed [m2/s2].
    !> @return Physical rate from the indicated constitutive closure.
    REAL(wp) FUNCTION eval_entrainment_volume_kernel(r_h, r_Ri, mod_vel2)
        REAL(wp), INTENT(IN) :: r_h, r_Ri, mod_vel2
        IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
                  ( r_Ri .GT. 0.0_wp ) ) THEN

              eval_entrainment_volume_kernel = eval_air_entrainment_rate(r_Ri, mod_vel2)

        ELSE

              eval_entrainment_volume_kernel = 0.0_wp

        END IF
    END FUNCTION eval_entrainment_volume_kernel
    !> @brief Evaluate the Richardson-number-dependent ambient-air entrainment law.
    !> @param[in] richardson_number Richardson number controlling entrainment [dimensionless].
    !> @param[in] speed_squared Squared three-dimensional flow speed [m2/s2].
    !> @return Physical rate from the indicated constitutive closure.
    PURE FUNCTION eval_air_entrainment_rate(richardson_number, speed_squared)  &
              RESULT(entrainment_rate)

        REAL(wp), INTENT(IN) :: richardson_number, speed_squared
        REAL(wp) :: entrainment_coefficient, entrainment_rate

        entrainment_coefficient = 0.075_wp / SQRT(1.0_wp + 718.0_wp *            &
                  richardson_number**2.4_wp)
        entrainment_rate = entrainment_coefficient * SQRT(speed_squared)

    END FUNCTION eval_air_entrainment_rate
END MODULE mass_exchange_entrainment_2d
