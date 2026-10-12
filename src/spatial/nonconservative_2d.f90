!********************************************************************************
!> \brief Generic quadrature engine for nonconservative path products
!********************************************************************************

MODULE nonconservative_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: PATH_DIR_X = 1
  INTEGER, PARAMETER, PUBLIC :: PATH_DIR_Y = 2

  PUBLIC :: eval_path_contribution
  PUBLIC :: nonconservative_integrand_callback
  PUBLIC :: path_state_callback

  ABSTRACT INTERFACE

     !***************************************************************************
     !> \brief Declare a model callback that evaluates the complete path integrand.
     !>
     !> \param[in] direction Cartesian path direction: PATH_DIR_X or PATH_DIR_Y.
     !> \param[in] path_state Extended state evaluated at the current path quadrature coordinate.
     !> \param[in] dstate_ds Derivative of the extended path state with respect to its unit-interval
     !>                      coordinate.
     !> \param[out] integrand Model one-form contracted with the path tangent, one entry per
     !>                       equation.
     !***************************************************************************

     SUBROUTINE nonconservative_integrand_callback( direction, path_state,  &
          dstate_ds, integrand )
       IMPORT :: wp
       INTEGER, INTENT(IN) :: direction
       REAL(wp), INTENT(IN) :: path_state(:)
       REAL(wp), INTENT(IN) :: dstate_ds(:)
       REAL(wp), INTENT(OUT) :: integrand(:)
     END SUBROUTINE nonconservative_integrand_callback

     !***************************************************************************
     !> \brief Declare a callback that evaluates a path state and its tangent.
     !>
     !> \param[in] stateL Extended state at the beginning of the model path.
     !> \param[in] stateR Extended state at the end of the model path.
     !> \param[in] s Path coordinate in the unit interval [0,1].
     !> \param[out] path_state Extended state evaluated at the current path quadrature coordinate.
     !> \param[out] dstate_ds Derivative of the extended path state with respect to its
     !>                       unit-interval coordinate.
     !***************************************************************************

     SUBROUTINE path_state_callback( stateL, stateR, s, path_state,          &
          dstate_ds )
       IMPORT :: wp
       REAL(wp), INTENT(IN) :: stateL(:)
       REAL(wp), INTENT(IN) :: stateR(:)
       REAL(wp), INTENT(IN) :: s
       REAL(wp), INTENT(OUT) :: path_state(:)
       REAL(wp), INTENT(OUT) :: dstate_ds(:)
     END SUBROUTINE path_state_callback

  END INTERFACE

CONTAINS

  !******************************************************************************
  !> \brief Integrate a model-provided nonconservative one-form along an endpoint path.
  !>
  !> The extended states are opaque to this module. By default they are joined
  !> by a straight segment. A model may supply a different path callback.
  !>
  !> \param[in] direction Cartesian path direction: PATH_DIR_X or PATH_DIR_Y.
  !> \param[in] stateL Extended state at the beginning of the model path.
  !> \param[in] stateR Extended state at the end of the model path.
  !> \param[in] integrand_provider Callback evaluating all nonconservative products at a quadrature
  !>                               state.
  !> \param[out] path_contribution Integrated nonconservative path contribution, one entry per
  !>                               equation.
  !> \param[in] n_quad Optional quadrature size (1,2,3); the integration default is three points.
  !> \param[in] path_provider Optional path/tangent callback; omitted to use the straight endpoint
  !>                          segment.
  !******************************************************************************

  SUBROUTINE eval_path_contribution( direction, stateL, stateR, integrand_provider, &
       path_contribution, n_quad, path_provider )

    INTEGER, INTENT(IN) :: direction
    REAL(wp), INTENT(IN) :: stateL(:)
    REAL(wp), INTENT(IN) :: stateR(:)
    PROCEDURE(nonconservative_integrand_callback) :: integrand_provider
    REAL(wp), INTENT(OUT) :: path_contribution(:)
    INTEGER, INTENT(IN), OPTIONAL :: n_quad
    PROCEDURE(path_state_callback), OPTIONAL :: path_provider

    REAL(wp) :: dstate_ds(SIZE(stateL))
    REAL(wp) :: gauss_abscissa(3)
    REAL(wp) :: gauss_weight(3)
    REAL(wp) :: integrand(SIZE(path_contribution))
    REAL(wp) :: path_state(SIZE(stateL))

    INTEGER :: i_quad
    INTEGER :: n_quad_local

    IF ( SIZE(stateR) .NE. SIZE(stateL) ) THEN
       ERROR STOP 'eval_path_contribution: inconsistent state sizes'
    END IF

    IF ( direction .NE. PATH_DIR_X .AND. direction .NE. PATH_DIR_Y ) THEN
       ERROR STOP 'eval_path_contribution: invalid direction'
    END IF

    n_quad_local = 3
    IF ( PRESENT(n_quad) ) n_quad_local = n_quad

    CALL path_quadrature_rule( n_quad_local, gauss_abscissa, gauss_weight )

    path_contribution = 0.0_wp

    DO i_quad = 1, n_quad_local

       IF ( PRESENT(path_provider) ) THEN
          CALL path_provider( stateL, stateR, gauss_abscissa(i_quad),       &
               path_state, dstate_ds )
       ELSE
          path_state = stateL + gauss_abscissa(i_quad) * (stateR - stateL)
          dstate_ds = stateR - stateL
       END IF

       integrand = 0.0_wp
       CALL integrand_provider( direction, path_state, dstate_ds, integrand )
       path_contribution = path_contribution                               &
            + gauss_weight(i_quad) * integrand

    END DO

  END SUBROUTINE eval_path_contribution


  !******************************************************************************
  !> \brief Return a one-, two- or three-point Gauss rule on the unit interval.
  !>
  !> \param[in] n_quad Number of Gauss points: one, two or three.
  !> \param[out] abscissa Quadrature coordinates on [0,1]; unused entries are zero.
  !> \param[out] weight Quadrature weights on [0,1]; unused entries are zero.
  !******************************************************************************

  SUBROUTINE path_quadrature_rule( n_quad, abscissa, weight )

    INTEGER, INTENT(IN) :: n_quad
    REAL(wp), INTENT(OUT) :: abscissa(3)
    REAL(wp), INTENT(OUT) :: weight(3)

    REAL(wp) :: offset

    abscissa = 0.0_wp
    weight = 0.0_wp

    SELECT CASE ( n_quad )

    CASE ( 1 )
       abscissa(1) = 0.5_wp
       weight(1) = 1.0_wp

    CASE ( 2 )
       offset = 0.5_wp / SQRT(3.0_wp)
       abscissa(1:2) = [ 0.5_wp - offset, 0.5_wp + offset ]
       weight(1:2) = 0.5_wp

    CASE ( 3 )
       offset = 0.5_wp * SQRT(3.0_wp / 5.0_wp)
       abscissa = [ 0.5_wp - offset, 0.5_wp, 0.5_wp + offset ]
       weight = [ 5.0_wp / 18.0_wp, 4.0_wp / 9.0_wp, 5.0_wp / 18.0_wp ]

    CASE DEFAULT
       ERROR STOP 'path_quadrature_rule: supported orders are 1, 2 and 3'

    END SELECT

  END SUBROUTINE path_quadrature_rule

END MODULE nonconservative_2d
