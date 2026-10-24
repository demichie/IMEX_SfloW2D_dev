!> \brief Common physical source assembly entry points.
!> \details Owns the explicit STEP source-rate bridge independently of the
!> spatial flux and hydrostatic interfaces in equation_terms_2d.
!> The bounded mass-exchange transaction and all inventory updates remain
!> the responsibility of step_source_executor_2d.
MODULE source_terms_2d
  USE parameters_2d, ONLY : wp, n_vars, n_eqns, n_solid
  USE mass_exchange_assembly_2d, ONLY : eval_mass_exchange_terms
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_step_source_explicit
CONTAINS
  !> \brief Evaluate the atomic explicit STEP mass-exchange rate.
  !> \details Returns rates only; it never mutates conservative state, bed or inventory.
  !> \param[in] qpj Cell-local primitive state.
  !> \param[in] B_zone Connected-water-region indicator.
  !> \param[in] B_prime_x Filtered bed x derivative.
  !> \param[in] B_prime_y Filtered bed y derivative.
  !> \param[in] erodible Per-class available substrate inventory.
  !> \param[in] dt Time-step length.
  !> \param[out] erosion_term Per-class substrate erosion rate.
  !> \param[out] deposition_term Per-class deposition rate.
  !> \param[out] continuous_phase_erosion_term Carrier entrainment rate.
  !> \param[out] continuous_phase_loss_term Limited carrier loss rate.
  !> \param[out] eqns_term Conservative equation source-rate vector.
  !> \param[out] topo_term Proposed bed-elevation rate.
  SUBROUTINE eval_step_source_explicit(qpj, B_zone, B_prime_x, B_prime_y, &
       erodible, dt, erosion_term, deposition_term, &
       continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2), B_prime_x, B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid), dt
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(OUT) :: erosion_term(n_solid), deposition_term(n_solid)
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term, continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns), topo_term

    CALL eval_mass_exchange_terms(qpj, B_zone, B_prime_x, B_prime_y, &
         erodible, dt, erosion_term, deposition_term, &
         continuous_phase_erosion_term, continuous_phase_loss_term, eqns_term, topo_term)
  END SUBROUTINE eval_step_source_explicit
END MODULE source_terms_2d
