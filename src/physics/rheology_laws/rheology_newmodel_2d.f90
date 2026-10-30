!> \brief D17 extension fixture: simple linear Newton-only drag.
!> \details Fixed coefficient k=0.1 s^-1. This demonstration law is intentionally
!> independent of model-specific namelist parameters and numerical assemblers.
MODULE rheology_newmodel_2d
    USE parameters_2d, ONLY : wp
    USE diagnostics_2d, ONLY : fatal_error
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: eval_rheology_newmodel_explicit
    PUBLIC :: eval_rheology_newmodel_semi_implicit
    PUBLIC :: eval_rheology_newmodel_implicit
    PUBLIC :: solve_rheology_newmodel_analytic
CONTAINS
    !> \brief Evaluate the explicit newmodel contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the configured newmodel rheology, with its currently supported source treatment, using the treatment-specific velocity and effective-stress inputs. Outputs are conservative x/y momentum source rates.
    !> \param[in] h Flow depth [m].
    !> \param[in] T Absolute mixture temperature [K].
    !> \param[in] u Horizontal x velocity [m s^-1].
    !> \param[in] v Horizontal y velocity [m s^-1].
    !> \param[in] mod_vel Three-dimensional speed [m s^-1].
    !> \param[in] mod_hor_vel Horizontal speed [m s^-1].
    !> \param[in] mod_vel2 Squared three-dimensional speed [m^2 s^-2].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[in] inv_rhom Reciprocal mixture density [m^3 kg^-1].
    !> \param[in] red_grav Reduced gravitational acceleration [m s^-2].
    !> \param[in] grav_coeff Local dimensionless gravitational correction.
    !> \param[in] alphas Solid volume fractions by class.
    !> \param[in] curvature_term Enable curvature correction.
    !> \param[in] B_second_xx Topographic second derivative in x [m^-1].
    !> \param[in] B_second_xy Mixed topographic derivative [m^-1].
    !> \param[in] B_second_yy Topographic second derivative in y [m^-1].
    !> \param[in] use_pore_pressure Enable pore-pressure correction.
    !> \param[in] exc_pore_pres Excess pore pressure [Pa].
    !> \param[in] moving True for a moving local state.
    !> \param[in] use_stochastic Enable model stochastic parameter.
    !> \param[in] stochastic_value Stochastic parameter value.
    !> \param[out] rate Two-component conservative momentum source.
    SUBROUTINE eval_rheology_newmodel_explicit( &
            h, T, u, v, mod_vel, mod_hor_vel, &
            mod_vel2, rhom, inv_rhom, red_grav, &
            grav_coeff, alphas, curvature_term, B_second_xx, B_second_xy, &
            B_second_yy, use_pore_pressure, exc_pore_pres, moving, use_stochastic, &
            stochastic_value, rate )
        REAL(wp), INTENT(IN) :: h, T, u, v
        REAL(wp), INTENT(IN) :: mod_vel, mod_hor_vel, mod_vel2
        REAL(wp), INTENT(IN) :: rhom, inv_rhom, red_grav
        REAL(wp), INTENT(IN) :: grav_coeff, alphas(:)
        REAL(wp), INTENT(IN) :: B_second_xx, B_second_xy, B_second_yy
        REAL(wp), INTENT(IN) :: exc_pore_pres, stochastic_value
        LOGICAL, INTENT(IN) :: curvature_term, use_pore_pressure, moving, use_stochastic
        REAL(wp), INTENT(OUT) :: rate(2)
        ! Keep the common four-treatment interface, but this fixture is Newton-only.

        rate=0.0_wp
        CALL fatal_error('rheology_newmodel: explicit treatment unsupported')
    END SUBROUTINE eval_rheology_newmodel_explicit
    !> \brief Evaluate the semi-implicit newmodel contribution to the horizontal momentum source.
    !> \details Return the contribution assigned to the semi-implicit treatment. This implementation supplies a zero source when the model has no term in this treatment.
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_newmodel_semi_implicit(rate_x,rate_y)
        REAL(wp), INTENT(OUT) :: rate_x,rate_y

        rate_x=0.0_wp
        rate_y=0.0_wp
    END SUBROUTINE eval_rheology_newmodel_semi_implicit
    !> \brief Evaluate the implicit newmodel contribution to the horizontal momentum source.
    !> \details Compute the momentum source under the configured newmodel rheology, with its currently supported source treatment, using the treatment-specific velocity and effective-stress inputs. Outputs are conservative x/y momentum source rates.
    !> \param[in] u Horizontal x velocity [m s^-1].
    !> \param[in] v Horizontal y velocity [m s^-1].
    !> \param[in] rhom Mixture density [kg m^-3].
    !> \param[out] rate_x x momentum source.
    !> \param[out] rate_y y momentum source.
    SUBROUTINE eval_rheology_newmodel_implicit(u,v, &
            rhom,rate_x,rate_y)
        COMPLEX(wp), INTENT(IN) :: u,v,rhom
        COMPLEX(wp), INTENT(OUT) :: rate_x,rate_y
        REAL(wp), PARAMETER :: drag_coefficient=0.1_wp

        rate_x=-drag_coefficient*rhom*u
        rate_y=-drag_coefficient*rhom*v
    END SUBROUTINE eval_rheology_newmodel_implicit
    !> \brief Integrate the newmodel velocity source in time.
    !> \details Advance the scalar speed under the configured newmodel rheology, with its currently supported source treatment with model coefficients fixed over dt; the horizontal direction is retained in the returned velocity vector.
    SUBROUTINE solve_rheology_newmodel_analytic()

        CALL fatal_error('rheology_newmodel: analytic treatment unsupported')
    END SUBROUTINE solve_rheology_newmodel_analytic

END MODULE rheology_newmodel_2d
