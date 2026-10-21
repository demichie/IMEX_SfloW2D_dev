!********************************************************************************
!> \brief Constant-strength plastic basal-friction model.
!> \details The model applies a speed-independent resisting strength
!> \f$\tau\f$ opposite to the horizontal flow direction:
!> \f[
!> \boldsymbol{S}=-\rho_m\tau\,\boldsymbol{u}_h/|\boldsymbol{u}|.
!> \f]
!> The supplied speed is the normalization used by the implementation. The
!> source is zero when the caller indicates that the flow is not moving.
!> \ref eval_plastic_implicit_source evaluates the complex-valued x/y sources
!> for the implicit solve; \ref tau stores the model strength.
!********************************************************************************

MODULE rheology_plastic_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tau, eval_plastic_implicit_source

  REAL(wp) :: tau

CONTAINS

  !******************************************************************************
  !> \brief Evaluate the implicit plastic-friction source.
  !> \details Returns \f$-\rho_m\tau(u_x,u_y)/|\boldsymbol{u}|\f$ for a
  !> moving state and zero otherwise. Complex arithmetic preserves derivatives
  !> with respect to the state in the implicit Jacobian.
  !> \param velocity_x x component of horizontal velocity.
  !> \param velocity_y y component of horizontal velocity.
  !> \param speed Speed used to normalize the direction.
  !> \param mixture_density Mixture density \f$\rho_m\f$.
  !> \param moving Whether the state is moving.
  !> \param source_x Returned x momentum source.
  !> \param source_y Returned y momentum source.
  !******************************************************************************
  SUBROUTINE eval_plastic_implicit_source(velocity_x, velocity_y, speed,       &
       mixture_density, moving, source_x, source_y)

    COMPLEX(wp), INTENT(IN) :: velocity_x, velocity_y, speed, mixture_density
    LOGICAL, INTENT(IN) :: moving
    COMPLEX(wp), INTENT(OUT) :: source_x, source_y

    source_x = CMPLX(0.0_wp, 0.0_wp, wp)
    source_y = CMPLX(0.0_wp, 0.0_wp, wp)
    IF (.NOT. moving) RETURN

    source_x = -mixture_density * tau * (velocity_x / speed)
    source_y = -mixture_density * tau * (velocity_y / speed)

  END SUBROUTINE eval_plastic_implicit_source

END MODULE rheology_plastic_2d
