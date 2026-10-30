!********************************************************************************
!> \brief State conversion and directly derived physical quantities
!>
!> Conservative states store mixture mass, momenta and thermal energy per unit area. Reconstructed
!> states store h, hu, hv, T, component mass fractions and appended u/v. Volume fractions are
!> derived only where constitutive laws require them.
!********************************************************************************

MODULE state_conversion_2d
  USE parameters_2d, ONLY : idx_hu, idx_hv, idx_totMassEqn, &
       idx_uEqn, idx_vEqn, idx_engyEqn

    USE parameters_2d, ONLY : idx_h, idx_T
  USE constitutive_parameters_2d
  USE diagnostics_2d, ONLY : fatal_error

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : n_vars, n_solid, n_add_gas
  USE parameters_2d, ONLY : eps_sing, eps_sing4, dry_thickness_tolerance
  USE parameters_2d, ONLY : liquid_flag, gas_flag,                              &
       stoch_transport_flag, pore_pressure_flag, sutherland_flag

  USE parameters_2d, ONLY : idx_solid_first, idx_solid_last, idx_add_gas_first,  &
       idx_add_gas_last, idx_stoch, idx_pore, idx_u, idx_v

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: mass_fraction_roundoff_tolerance =                 &
       4096.0_wp * EPSILON(1.0_wp)

  !> Shared local physical state for non-Newton source evaluation.
  !> Allocatable phase arrays accommodate the runtime model layout.
  TYPE, PUBLIC :: physical_state_type
     REAL(wp) :: h=0.0_wp, u=0.0_wp, v=0.0_wp, T=0.0_wp
     REAL(wp) :: rho_m=0.0_wp, inv_rhom=0.0_wp, red_grav=0.0_wp
     REAL(wp) :: Zs=0.0_wp, exc_pore_pres=0.0_wp
     REAL(wp), ALLOCATABLE :: alphas(:), alphag(:)
  END TYPE physical_state_type

  !> Corresponding complex-step state; no loss of imaginary perturbations.
  TYPE, PUBLIC :: complex_physical_state_type
     COMPLEX(wp) :: h=(0.0_wp,0.0_wp), u=(0.0_wp,0.0_wp), v=(0.0_wp,0.0_wp)
     COMPLEX(wp) :: T=(0.0_wp,0.0_wp), rho_m=(0.0_wp,0.0_wp)
     COMPLEX(wp) :: inv_rhom=(0.0_wp,0.0_wp), red_grav=(0.0_wp,0.0_wp)
     COMPLEX(wp) :: Zs=(0.0_wp,0.0_wp), exc_pore_pres=(0.0_wp,0.0_wp)
     COMPLEX(wp), ALLOCATABLE :: alphas(:), alphag(:)
  END TYPE complex_physical_state_type

  !> \brief Evaluate carrier and whole-mixture specific heat capacities.
  !> \param[in]  xs            Solid mass fractions.
  !> \param[in]  xg            Additional-gas mass fractions.
  !> \param[out] xc            Total carrier mass fraction.
  !> \param[out] sp_heat_c_mix Carrier-mixture specific heat capacity.
  !> \param[out] sp_heat_mix   Whole-mixture specific heat capacity.
  INTERFACE eval_mixture_heat_capacity
     MODULE PROCEDURE eval_mixture_heat_capacity_real
     MODULE PROCEDURE eval_mixture_heat_capacity_complex
  END INTERFACE eval_mixture_heat_capacity

  !> \brief Evaluate mixture properties from component mass fractions.
  !> \param[in]  T         Mixture temperature.
  !> \param[in]  xs        Solid mass fractions.
  !> \param[in]  xg        Additional-gas mass fractions.
  !> \param[out] rho_m     Mixture density.
  !> \param[out] inv_rhom  Inverse mixture density.
  !> \param[out] rho_c     Carrier density.
  !> \param[out] inv_rho_c Inverse carrier density.
  !> \param[out] alphas    Solid volume fractions.
  !> \param[out] alphag    Additional-gas volume fractions.
  INTERFACE eval_mixture_properties_from_mass_fractions
     MODULE PROCEDURE eval_mixture_properties_from_mass_real
     MODULE PROCEDURE eval_mixture_properties_from_mass_complex
  END INTERFACE eval_mixture_properties_from_mass_fractions

  !> Recover velocity from conservative mass and momenta with the same
  !> near-dry desingularization in the REAL and COMPLEX code paths.
  INTERFACE velocity_from_conservative
     MODULE PROCEDURE velocity_from_conservative_real
     MODULE PROCEDURE velocity_from_conservative_complex
  END INTERFACE velocity_from_conservative

  PUBLIC :: r_phys_var, c_phys_var
  PUBLIC :: physical_state_from_qp, physical_state_from_q, real_part_physical_state
  PUBLIC :: qc_to_qp, qp_to_qc, qp_to_qp2
  PUBLIC :: mixt_var
  PUBLIC :: sauter_diameter, average_density_solids
  PUBLIC :: eval_mixture_heat_capacity
  PUBLIC :: eval_mixture_properties_from_mass_fractions
  PUBLIC :: eval_mixture_properties_from_volume_fractions
  PUBLIC :: enforce_mass_fraction_closure
  PUBLIC :: enforce_primitive_mass_fraction_closure
  PUBLIC :: clip_roundoff_negative_component_masses
  PUBLIC :: primitive_to_volume_fractions
  PUBLIC :: velocity_from_conservative
  PUBLIC :: project_dry_cell_momenta

CONTAINS

  !> Build the shared REAL source state directly from the available primitive vector.
  !> No conservative-to-primitive conversion is performed on this path.
  SUBROUTINE physical_state_from_qp(qp, state)
    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    TYPE(physical_state_type), INTENT(OUT) :: state
    REAL(wp) :: Ri, rho_c, sp_heat_c, sp_heat_mix

    ALLOCATE(state%alphas(n_solid), state%alphag(n_add_gas))
    state%h = qp(idx_h)
    state%u = qp(idx_u)
    state%v = qp(idx_v)
    state%T = qp(idx_T)
    state%Zs = 0.0_wp
    state%exc_pore_pres = 0.0_wp
    IF (stoch_transport_flag) state%Zs = qp(idx_stoch)
    IF (pore_pressure_flag) state%exc_pore_pres = qp(idx_pore)
    CALL mixt_var(qp, Ri, state%rho_m, rho_c, state%red_grav, &
         sp_heat_c, sp_heat_mix)
    state%inv_rhom = 1.0_wp / state%rho_m
    CALL primitive_to_volume_fractions(qp, state%alphas, state%alphag)
  END SUBROUTINE physical_state_from_qp

  !> Build a COMPLEX physical source state through the existing complex-step kernel.
  SUBROUTINE physical_state_from_q(q, state)
    COMPLEX(wp), INTENT(IN) :: q(n_vars)
    TYPE(complex_physical_state_type), INTENT(OUT) :: state

    ALLOCATE(state%alphas(n_solid), state%alphag(n_add_gas))
    CALL c_phys_var(q, state%h, state%u, state%v, state%T, state%rho_m, &
         state%alphas, state%alphag, state%inv_rhom, state%Zs, &
         state%exc_pore_pres)
    state%red_grav = (state%rho_m - rho_a_amb) * state%inv_rhom * grav
  END SUBROUTINE physical_state_from_q

  !> Extract a real source state without re-evaluating any constitutive formula.
  SUBROUTINE real_part_physical_state(complex_state, state)
    TYPE(complex_physical_state_type), INTENT(IN) :: complex_state
    TYPE(physical_state_type), INTENT(OUT) :: state

    state%h = REAL(complex_state%h,wp)
    state%u = REAL(complex_state%u,wp)
    state%v = REAL(complex_state%v,wp)
    state%T = REAL(complex_state%T,wp)
    state%rho_m = REAL(complex_state%rho_m,wp)
    state%inv_rhom = REAL(complex_state%inv_rhom,wp)
    state%red_grav = REAL(complex_state%red_grav,wp)
    state%Zs = REAL(complex_state%Zs,wp)
    state%exc_pore_pres = REAL(complex_state%exc_pore_pres,wp)
    state%alphas = REAL(complex_state%alphas,wp)
    state%alphag = REAL(complex_state%alphag,wp)
  END SUBROUTINE real_part_physical_state

  !> \brief Remove unresolved cell momenta without discarding conservative mass.
  !> \param[in,out] qc Conservative cell state; only entries 2 and 3 may change.
  !> \param[in] h Depth recovered with the current mixture thermodynamic closure.
  !> \details At h <= the existing shared dry-depth threshold, reset both
  !>          momenta to zero. Mass, thermal energy and transported components
  !>          are bitwise unchanged. Wet momenta are also unchanged.
  !> \note Apply after observing/checking the raw final IMEX assembly, not inside
  !>       the differentiable real/complex closure or the implicit Newton solve.
  PURE SUBROUTINE project_dry_cell_momenta(qc, h)
    REAL(wp), INTENT(INOUT) :: qc(:)
    REAL(wp), INTENT(IN) :: h

    IF (h <= dry_thickness_tolerance) qc(idx_uEqn:idx_vEqn) = 0.0_wp
  END SUBROUTINE project_dry_cell_momenta

  !> \brief Recover real velocities from mixture mass and momenta with dry-state regularization.
  !>
  !> \param[in] qc Local conservative vector: mixture mass, momenta, thermal energy and transported
  !>               component masses.
  !> \param[out] u Cartesian depth-averaged x velocity [m s^-1].
  !> \param[out] v Cartesian depth-averaged y velocity [m s^-1].

  SUBROUTINE velocity_from_conservative_real(qc, u, v)

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: u, v

    INCLUDE 'velocity_from_conservative.inc'

  END SUBROUTINE velocity_from_conservative_real

  !> \brief Recover complex-step velocities using the same real-path regularization.
  !>
  !> \param[in] qc Local conservative vector: mixture mass, momenta, thermal energy and transported
  !>               component masses.
  !> \param[out] u Cartesian depth-averaged x velocity [m s^-1].
  !> \param[out] v Cartesian depth-averaged y velocity [m s^-1].

  SUBROUTINE velocity_from_conservative_complex(qc, u, v)

    COMPLEX(wp), INTENT(IN) :: qc(n_vars)
    COMPLEX(wp), INTENT(OUT) :: u, v

    INCLUDE 'velocity_from_conservative.inc'

  END SUBROUTINE velocity_from_conservative_complex

  !> \brief Clip roundoff negatives and proportionally close explicit component fractions.
  !>
  !> Enforce the admissible closure of the independently transported mass
  !> fractions. Roundoff-sized negative values are clipped; material negative
  !> values indicate a reconstruction/configuration error. If the explicit
  !> component sum exceeds one, all components are normalized proportionally.
  !>
  !> \param[in,out] xs Solid mass fractions, one entry per solid class.
  !> \param[in,out] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[in] context Optional diagnostic context included in closure errors.

  SUBROUTINE enforce_mass_fraction_closure(xs, xg, context)

    REAL(wp), INTENT(INOUT) :: xs(n_solid), xg(n_add_gas)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: context
    ! A few thousand machine epsilons allow accumulated cancellation noise in
    ! the explicit FV update to be clipped, while materially negative
    ! fractions still trigger the fatal check below.
    REAL(wp) :: explicit_mass_fraction

    IF (n_solid .GT. 0) THEN
       IF (ANY(xs .LT. -mass_fraction_roundoff_tolerance)) THEN
          CALL mass_fraction_closure_error(                                  &
               'Material negative solid mass fraction', context)
       END IF
       xs = MAX(xs, 0.0_wp)
    END IF

    IF (n_add_gas .GT. 0) THEN
       IF (ANY(xg .LT. -mass_fraction_roundoff_tolerance))                    &
            CALL mass_fraction_closure_error(                                &
                 'Material negative additional-gas mass fraction', context)
       xg = MAX(xg, 0.0_wp)
    END IF

    ! Only one continuous carrier is transported.

    explicit_mass_fraction = SUM(xs) + SUM(xg)
    IF (explicit_mass_fraction .GT. 1.0_wp) THEN
       xs = xs / explicit_mass_fraction
       xg = xg / explicit_mass_fraction
    END IF

  END SUBROUTINE enforce_mass_fraction_closure


  !> \brief Apply admissible mass-fraction closure to a physical state vector.
  !>
  !> Apply the mass-fraction closure directly to a primitive state vector.
  !>
  !> \param[in,out] qp Canonical physical states: h, hu, hv, T, component mass fractions, then
  !>                   appended u and v.
  !> \param[in] context Optional diagnostic context included in closure errors.

  SUBROUTINE enforce_primitive_mass_fraction_closure(qp, context)

    REAL(wp), INTENT(INOUT) :: qp(n_vars+2)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: context
    REAL(wp) :: xs(n_solid), xg(n_add_gas)

    xs = qp(idx_solid_first:idx_solid_last)
    xg = qp(idx_add_gas_first:idx_add_gas_last)

    IF (PRESENT(context)) THEN
       CALL enforce_mass_fraction_closure(xs, xg, context)
    ELSE
       CALL enforce_mass_fraction_closure(xs, xg)
    END IF

    qp(idx_solid_first:idx_solid_last) = xs
    qp(idx_add_gas_first:idx_add_gas_last) = xg

  END SUBROUTINE enforce_primitive_mass_fraction_closure

  !> \brief Clip only roundoff-sized negative conservative component masses.
  !> Material negatives are deliberately left unchanged for the existing
  !> admissibility checks to report as errors.
  !> \param[in,out] qc Local conservative state; only roundoff-sized negative
  !>                   solid, gas or liquid masses may be changed.
  SUBROUTINE clip_roundoff_negative_component_masses(qc)

    REAL(wp), INTENT(INOUT) :: qc(n_vars)
    REAL(wp) :: tolerance

    IF (qc(idx_totMassEqn) .LE. 0.0_wp) RETURN
    tolerance = mass_fraction_roundoff_tolerance *                         &
         MAX(ABS(qc(idx_totMassEqn)),TINY(1.0_wp))

    IF (n_solid .GT. 0) THEN
       WHERE ((qc(idx_solid_first:idx_solid_last) .LT. 0.0_wp) .AND.       &
            (qc(idx_solid_first:idx_solid_last) .GE. -tolerance))
          qc(idx_solid_first:idx_solid_last) = 0.0_wp
       END WHERE
    END IF
    IF (n_add_gas .GT. 0) THEN
       WHERE ((qc(idx_add_gas_first:idx_add_gas_last) .LT. 0.0_wp) .AND.  &
            (qc(idx_add_gas_first:idx_add_gas_last) .GE. -tolerance))
          qc(idx_add_gas_first:idx_add_gas_last) = 0.0_wp
       END WHERE
    END IF

  END SUBROUTINE clip_roundoff_negative_component_masses


  !> \brief Derive constitutive volume fractions from canonical physical mass fractions.
  !>
  !> Derive volume fractions from the canonical primitive mass fractions.
  !>
  !> \param[in] qp Canonical physical states: h, hu, hv, T, component mass fractions, then appended
  !>               u and v.
  !> \param[out] alphas Solid volume fractions, one entry per solid class.
  !> \param[out] alphag Additional-gas volume fractions, excluding ambient air.

  SUBROUTINE primitive_to_volume_fractions(qp, alphas, alphag)

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)

    REAL(wp) :: xs(n_solid), xg(n_add_gas)
    REAL(wp) :: rho_m, inv_rhom, rho_c, inv_rho_c

    alphas = 0.0_wp
    alphag = 0.0_wp

    IF (qp(idx_h) .LE. EPSILON(1.0_wp)) RETURN

    xs = qp(idx_solid_first:idx_solid_last)
    xg = qp(idx_add_gas_first:idx_add_gas_last)

    CALL enforce_mass_fraction_closure(xs, xg)
    CALL eval_mixture_properties_from_mass_fractions(                         &
         qp(idx_T), xs, xg, rho_m, inv_rhom, rho_c, inv_rho_c,               &
         alphas, alphag)

  END SUBROUTINE primitive_to_volume_fractions

  !> \brief Compute the volume-weighted Sauter diameter of the solid mixture.
  !>
  !> Function that calculates the Sauter diameter
  !>
  !> \param[in] alpha_solids Depth-averaged solid volume fractions used as mixture weights.
  !> \return Sauter diameter [m]; an unweighted mean diameter when solid volume is negligible.

  FUNCTION sauter_diameter(alpha_solids)
      !> Sauter diameter
      REAL(wp) :: sauter_diameter

      !> depth-averaged volumetric fractions of solid phases
      REAL(wp), INTENT(IN) :: alpha_solids(n_solid)

      IF ( SUM( alpha_solids ) .LT. EPSILON(1.0_wp) ) THEN
         sauter_diameter = SUM( diam_s ) / n_solid
      ELSE
         sauter_diameter = SUM( alpha_solids ) / SUM( alpha_solids /  diam_s )
      END IF

  END FUNCTION sauter_diameter

  !> \brief Compute the volume-weighted mean density of the solid phases.
  !>
   !> Function that calculates the average density of solid phases
  !>
  !> \param[in] alpha_solids Depth-averaged solid volume fractions used as mixture weights.
  !> \return Solid mixture density [kg m^-3]; an unweighted mean when solid volume is negligible.

  FUNCTION average_density_solids(alpha_solids)
      !> Average density
      REAL(wp) :: average_density_solids

      !> depth-averaged volumetric fractions of solid phases
      REAL(wp), INTENT(IN) :: alpha_solids(n_solid)

      IF ( SUM( alpha_solids ) .LT. EPSILON(1.0_wp) ) THEN
         average_density_solids = SUM( rho_s ) / n_solid
      ELSE
         average_density_solids = SUM( alpha_solids * rho_s ) / SUM( alpha_solids )
      END IF

  END FUNCTION average_density_solids

  !> \brief Evaluate real carrier and whole-mixture specific heat capacities.
  !>
  !> \param[in] xs Solid mass fractions, one entry per solid class.
  !> \param[in] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[out] xc Total carrier mass fraction, including any additional gas components.
  !> \param[out] sp_heat_c_mix Carrier-mixture specific heat capacity [J kg^-1 K^-1].
  !> \param[out] sp_heat_mix Whole-mixture specific heat capacity [J kg^-1 K^-1].

  SUBROUTINE eval_mixture_heat_capacity_real(xs, xg, xc, sp_heat_c_mix,     &
       sp_heat_mix)

    REAL(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
    REAL(wp), INTENT(OUT) :: xc, sp_heat_c_mix, sp_heat_mix
    REAL(wp) :: carrier_sp_heat_mass

    INCLUDE 'mixture_heat_capacity.inc'

  END SUBROUTINE eval_mixture_heat_capacity_real

  !> \brief Evaluate complex-step carrier and whole-mixture specific heat capacities.
  !>
  !> \param[in] xs Solid mass fractions, one entry per solid class.
  !> \param[in] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[out] xc Total carrier mass fraction, including any additional gas components.
  !> \param[out] sp_heat_c_mix Carrier-mixture specific heat capacity [J kg^-1 K^-1].
  !> \param[out] sp_heat_mix Whole-mixture specific heat capacity [J kg^-1 K^-1].

  SUBROUTINE eval_mixture_heat_capacity_complex(xs, xg, xc,                &
       sp_heat_c_mix, sp_heat_mix)

    COMPLEX(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
    COMPLEX(wp), INTENT(OUT) :: xc, sp_heat_c_mix, sp_heat_mix
    COMPLEX(wp) :: carrier_sp_heat_mass

    INCLUDE 'mixture_heat_capacity.inc'

  END SUBROUTINE eval_mixture_heat_capacity_complex

  !> \brief Derive real densities and volume fractions from component mass fractions.
  !>
  !> \param[in] T Mixture temperature [K].
  !> \param[in] xs Solid mass fractions, one entry per solid class.
  !> \param[in] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[out] rho_m Whole-mixture density [kg m^-3].
  !> \param[out] inv_rhom Reciprocal of the whole-mixture density [m^3 kg^-1].
  !> \param[out] rho_c Density of the carrier mixture [kg m^-3].
  !> \param[out] inv_rho_c Reciprocal of carrier-mixture density [m^3 kg^-1].
  !> \param[out] alphas Solid volume fractions, one entry per solid class.
  !> \param[out] alphag Additional-gas volume fractions, excluding ambient air.

  SUBROUTINE eval_mixture_properties_from_mass_real(T, xs, xg, rho_m,     &
       inv_rhom, rho_c, inv_rho_c, alphas, alphag)

    REAL(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
    REAL(wp), INTENT(OUT) :: rho_m, inv_rhom, rho_c, inv_rho_c
    REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    REAL(wp) :: xc, sp_gas_const_c, inv_rho_g(n_add_gas)

    INCLUDE 'mixture_from_mass.inc'

  END SUBROUTINE eval_mixture_properties_from_mass_real

  !> \brief Derive complex-step densities and volume fractions from mass fractions.
  !>
  !> \param[in] T Mixture temperature [K].
  !> \param[in] xs Solid mass fractions, one entry per solid class.
  !> \param[in] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[out] rho_m Whole-mixture density [kg m^-3].
  !> \param[out] inv_rhom Reciprocal of the whole-mixture density [m^3 kg^-1].
  !> \param[out] rho_c Density of the carrier mixture [kg m^-3].
  !> \param[out] inv_rho_c Reciprocal of carrier-mixture density [m^3 kg^-1].
  !> \param[out] alphas Solid volume fractions, one entry per solid class.
  !> \param[out] alphag Additional-gas volume fractions, excluding ambient air.

  SUBROUTINE eval_mixture_properties_from_mass_complex(T, xs, xg, rho_m,  &
       inv_rhom, rho_c, inv_rho_c, alphas, alphag)

    COMPLEX(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
    COMPLEX(wp), INTENT(OUT) :: rho_m, inv_rhom, rho_c, inv_rho_c
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    COMPLEX(wp) :: xc, sp_gas_const_c, inv_rho_g(n_add_gas)

    INCLUDE 'mixture_from_mass.inc'

  END SUBROUTINE eval_mixture_properties_from_mass_complex


  !> \brief Convert volume fractions to mass fractions, densities and heat capacities.
  !>
  !> \param[in] T Mixture temperature [K].
  !> \param[in] alphag Additional-gas volume fractions, excluding ambient air.
  !> \param[in,out] alphas Solid volume fractions; rescaled if the dispersed sum exceeds one,
  !>                       or reset if that sum is negative.
  !>                       otherwise reset to zero.
  !> \param[out] rho_m Whole-mixture density [kg m^-3].
  !> \param[out] inv_rhom Reciprocal of the whole-mixture density [m^3 kg^-1].
  !> \param[out] rho_c Density of the carrier mixture [kg m^-3].
  !> \param[out] xs Solid mass fractions, one entry per solid class.
  !> \param[out] xg Additional-gas mass fractions, excluding ambient air.
  !> \param[out] xc Total carrier mass fraction, including any additional gas components.
  !> \param[out] sp_heat_c_mix Carrier-mixture specific heat capacity [J kg^-1 K^-1].
  !> \param[out] sp_heat_mix Whole-mixture specific heat capacity [J kg^-1 K^-1].

  SUBROUTINE eval_mixture_properties_from_volume_fractions(T, alphag, alphas, &
        rho_m, inv_rhom, rho_c, xs, xg, xc, sp_heat_c_mix,          &
       sp_heat_mix)

    REAL(wp), INTENT(IN) :: T
    REAL(wp), INTENT(IN) :: alphag(n_add_gas)
    REAL(wp), INTENT(INOUT) :: alphas(n_solid)
    REAL(wp), INTENT(OUT) :: rho_m, inv_rhom, rho_c
    REAL(wp), INTENT(OUT) :: xs(n_solid), xg(n_add_gas), xc
    REAL(wp), INTENT(OUT) :: sp_heat_c_mix, sp_heat_mix

    REAL(wp) :: alphas_tot, dispersed_tot, alphac, alpha_air
    REAL(wp) :: rho_a, rho_g(n_add_gas)

    alphas_tot = SUM(alphas)
    dispersed_tot = alphas_tot

    IF ( dispersed_tot .GT. 1.0_wp ) THEN

       alphas = alphas / dispersed_tot

    ELSEIF ( dispersed_tot .LT. 0.0_wp ) THEN

       alphas = 0.0_wp

    END IF

    alphas_tot = SUM(alphas)
    alphac = 1.0_wp - alphas_tot

    IF ( gas_flag ) THEN

       rho_a = pres / ( sp_gas_const_a * T )
       rho_g = pres / ( sp_gas_const_g * T )
       alpha_air = alphac - SUM(alphag)

       IF ( alphac .GT. EPSILON(1.0_wp) ) THEN
          rho_c = ( alpha_air * rho_a + SUM(alphag * rho_g) ) / alphac
       ELSE
          rho_c = rho_a
       END IF

    ELSE

       rho_c = rho_l
       rho_g = 0.0_wp

    END IF

    rho_m = SUM(alphas * rho_s) + alphac * rho_c
    inv_rhom = 1.0_wp / rho_m

    xs = alphas * rho_s * inv_rhom
    xg = alphag * rho_g * inv_rhom

    CALL eval_mixture_heat_capacity_real(                                      &
         ! IN
         xs, xg,                                                           &
         ! OUT
         xc, sp_heat_c_mix, sp_heat_mix)

  END SUBROUTINE eval_mixture_properties_from_volume_fractions

  !******************************************************************************
  !> \brief Decode conservative variables into real thermodynamic and kinematic quantities.
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !>
  !> \param[in] qj Local conservative state vector; mass and momenta are per unit horizontal area.
  !> \param[out] h Flow depth [m].
  !> \param[out] u Cartesian depth-averaged x velocity [m s^-1].
  !> \param[out] v Cartesian depth-averaged y velocity [m s^-1].
  !> \param[out] alphas Solid volume fractions, one entry per solid class.
  !> \param[out] rho_m Whole-mixture density [kg m^-3].
  !> \param[out] T Mixture temperature [K].
  !> \param[out] alphag Additional-gas volume fractions, excluding ambient air.
  !> \param[out] red_grav Buoyancy-reduced gravitational acceleration [m s^-2].
  !> \param[out] p_dyn Dynamic pressure rho_m*(u^2+v^2)/2 [Pa].
  !> \param[out] Zs Transported scalar stochastic fluctuation decoded from its conservative mass.
  !> \param[out] exc_pore_pres Transported excess pore pressure decoded from its conservative
  !>                           variable.
  !******************************************************************************

  SUBROUTINE r_phys_var(qj, h, u, v, alphas, rho_m, T, alphag,          &
       red_grav, p_dyn, Zs, exc_pore_pres)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qj(n_vars)
    REAL(wp), INTENT(OUT) :: h, u, v, T, rho_m
    REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    REAL(wp), INTENT(OUT) :: red_grav, p_dyn, Zs, exc_pore_pres

    REAL(wp) :: inv_rhom, rho_c, inv_rho_c
    REAL(wp) :: xs(n_solid), xg(n_add_gas)
    REAL(wp) :: xc, inv_qj1
    REAL(wp) :: mixture_sp_heat_c, sp_heat_mix
    LOGICAL :: is_wet

    ! The shared fragment returns immediately for a dry state, so initialize
    ! the REAL-only derived outputs before entering it.
    red_grav = 0.0_wp
    p_dyn = 0.0_wp

    INCLUDE 'state_conversion_phys_var.inc'

    red_grav = ( rho_m - rho_a_amb ) * inv_rhom * grav

    p_dyn = 0.5_wp * rho_m * ( u**2 + v**2 )

  END SUBROUTINE r_phys_var


  !******************************************************************************
  !> \brief Decode a perturbed conservative state for complex-step source differentiation.
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the complex thermodynamic and kinematic variables used for Jacobian evaluation.
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !>
  !> \param[in] qj Local conservative state vector; mass and momenta are per unit horizontal area.
  !> \param[out] h Flow depth [m].
  !> \param[out] u Cartesian depth-averaged x velocity [m s^-1].
  !> \param[out] v Cartesian depth-averaged y velocity [m s^-1].
  !> \param[out] T Mixture temperature [K].
  !> \param[out] rho_m Whole-mixture density [kg m^-3].
  !> \param[out] alphas Solid volume fractions, one entry per solid class.
  !> \param[out] alphag Additional-gas volume fractions, excluding ambient air.
  !> \param[out] inv_rhom Reciprocal of the whole-mixture density [m^3 kg^-1].
  !> \param[out] Zs Transported scalar stochastic fluctuation decoded from its conservative mass.
  !> \param[out] exc_pore_pres Transported excess pore pressure decoded from its conservative
  !>                           variable.
  !******************************************************************************

  SUBROUTINE c_phys_var(qj, h, u, v, T, rho_m, alphas, alphag, inv_rhom, Zs,    &
       exc_pore_pres)

    IMPLICIT NONE

    COMPLEX(wp), INTENT(IN) :: qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h, u, v, T, rho_m
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
    COMPLEX(wp), INTENT(OUT) :: inv_rhom, Zs, exc_pore_pres

    COMPLEX(wp) :: rho_c, inv_rho_c
    COMPLEX(wp) :: xs(n_solid), xg(n_add_gas)
    COMPLEX(wp) :: xc, inv_qj1
    COMPLEX(wp) :: mixture_sp_heat_c, sp_heat_mix
    LOGICAL :: is_wet

    INCLUDE 'state_conversion_phys_var.inc'

  END SUBROUTINE c_phys_var


  !******************************************************************************
  !> \brief Evaluate density, reduced gravity, heat capacities and Richardson number from a physical
  !>        state.
  !
  !> This subroutine evaluates from the physical real-value local variables qpj,
  !> some mixture variable.
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !> \date 10/10/2019
  !>
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[out] r_Ri Local dimensionless Richardson number.
  !> \param[out] r_rho_m Whole-mixture density [kg m^-3].
  !> \param[out] r_rho_c Carrier-mixture density [kg m^-3].
  !> \param[out] r_red_grav Buoyancy-reduced gravitational acceleration [m s^-2].
  !> \param[out] r_sp_heat_c Carrier-mixture specific heat capacity [J kg^-1 K^-1].
  !> \param[out] r_sp_heat_mix Whole-mixture specific heat capacity [J kg^-1 K^-1].
  !******************************************************************************

  SUBROUTINE mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,r_sp_heat_c,         &
       r_sp_heat_mix)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(wp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    REAL(wp), INTENT(OUT) :: r_sp_heat_c
    REAL(wp), INTENT(OUT) :: r_sp_heat_mix


    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness
    REAL(wp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)       !< real-value add.gas volume fractions
    REAL(wp) :: r_T                       !< real-value temperature [K]
    REAL(wp) :: r_inv_rhom, r_inv_rho_c
    REAL(wp) :: r_xs(n_solid), r_xg(n_add_gas)
    REAL(wp) :: r_xc
    REAL(wp) :: speed_squared, buoyancy_depth

    r_h = qpj(idx_h)

    IF ( qpj(idx_h) .LE. EPSILON(1.0_wp) ) THEN

       r_red_grav = 0.0_wp
       r_rho_m = rho_a_amb
       r_rho_c = rho_a_amb

       IF ( liquid_flag ) THEN

          r_sp_heat_c = sp_heat_l
          r_sp_heat_mix = sp_heat_l

       ELSE

          r_sp_heat_c = sp_heat_a
          r_sp_heat_mix = sp_heat_a


       END IF

       r_Ri = 0.0_wp

       RETURN

    END IF

    r_u = qpj(idx_u)
    r_v = qpj(idx_v)
    r_T = qpj(idx_T)

    r_xs = qpj(idx_solid_first:idx_solid_last)
    r_xg = qpj(idx_add_gas_first:idx_add_gas_last)
    CALL enforce_mass_fraction_closure(r_xs, r_xg)

    CALL eval_mixture_properties_from_mass_fractions(                          &
         ! IN
         r_T, r_xs, r_xg,                                                &
         ! OUT
         r_rho_m, r_inv_rhom, r_rho_c, r_inv_rho_c,                           &
         r_alphas, r_alphag)

    CALL eval_mixture_heat_capacity(r_xs, r_xg, r_xc,                   &
         r_sp_heat_c, r_sp_heat_mix)

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! Apply the existing upper cap BEFORE division: almost-resting wet cells
    ! can have a positive subnormal speed_squared whose reciprocal overflows.
    ! MIN(cap,numerator/denominator) evaluates that overflowing quotient first.
    speed_squared = r_u**2 + r_v**2
    buoyancy_depth = r_red_grav * r_h
    IF ( speed_squared .GT. 0.0_wp ) THEN

       IF ( buoyancy_depth / 1.E15_wp .GT. speed_squared ) THEN
          r_Ri = 1.E15_wp
       ELSE
          r_Ri = MIN(1.E15_wp,buoyancy_depth / speed_squared)
       END IF

    ELSE

       r_Ri = 0.0_wp

    END IF

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Convert a local conservative state into the canonical physical reconstruction state.
  !
  !> This subroutine evaluates from the conservative variables qc the
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_solid_first:idx_solid_last) = \f$ Ys(1:n_solid) \f$
  !> - qp(idx_add_gas_first:idx_add_gas_last) = \f$ Yg(1:n_add_gas) \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !
  !> \date 2019/11/11
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !>
  !> \param[in] qc Local conservative vector: mixture mass, momenta, thermal energy and transported
  !>               component masses.
  !> \param[out] qp Canonical physical states: h, hu, hv, T, component mass fractions, then appended
  !>                u and v.
  !> \param[out] p_dyn Dynamic pressure rho_m*(u^2+v^2)/2 [Pa].
  !> \param[in] context Optional diagnostic context included in conversion errors.
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp,p_dyn,context)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: p_dyn
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: context

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_Zs             !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres  !< real-value pore pressure
    REAL(wp) :: r_xs(n_solid), r_xg(n_add_gas)

    qp = 0.0_wp

    CALL r_phys_var( qc , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,          &
         r_alphag , r_red_grav , p_dyn , r_Zs , r_exc_pore_pres )

    qp(idx_h) = r_h

    qp(idx_hu) = r_h*r_u
    qp(idx_hv) = r_h*r_v

    qp(idx_T) = r_T

    r_xs = 0.0_wp
    r_xg = 0.0_wp
    IF ( qc(idx_totMassEqn) .GT. EPSILON(1.0_wp) ) THEN
       r_xs = qc(idx_solid_first:idx_solid_last) / qc(idx_totMassEqn)
       r_xg = qc(idx_add_gas_first:idx_add_gas_last) / qc(idx_totMassEqn)
       IF (PRESENT(context)) THEN
          CALL enforce_mass_fraction_closure(r_xs, r_xg, context)
       ELSE
          CALL enforce_mass_fraction_closure(r_xs, r_xg)
       END IF
    END IF

    qp(idx_solid_first:idx_solid_last) = r_xs
    qp(idx_add_gas_first:idx_add_gas_last) = r_xg

    IF ( stoch_transport_flag) qp(idx_stoch) = r_Zs

    IF ( pore_pressure_flag) qp(idx_pore) = r_exc_pore_pres

    qp(idx_u) = r_u
    qp(idx_v) = r_v

    RETURN

  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Convert a physical mass-fraction state into conservative mass, momenta and thermal
  !>        energy.
  !
  !> This subroutine evaluates the conservative real_value variables qc from the
  !> array of real_valued physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_solid_first:idx_solid_last) = \f$ Ys(1:n_s) \f$
  !> - qp(idx_add_gas_first:idx_add_gas_last) = \f$ Yg(1:n_g) \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !
  !> \date 2019/11/18
  !
  !> @author
  !> Mattia de' Michieli Vitturi
  !
  !>
  !> \param[in] qp Canonical physical states: h, hu, hv, T, component mass fractions, then appended
  !>               u and v.
  !> \param[out] qc Local conservative vector: mixture mass, momenta, thermal energy and transported
  !>                component masses.
  !> \param[in] context Optional diagnostic context included in conversion errors.
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,qc,context)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: qc(n_vars)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: context

    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_hu              !< real-value volumetric x-flow
    REAL(wp) :: r_hv              !< real-value volumetric y-flow
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)!< real-value add.gas volume fractions
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)   !< real-value add.gas mass fractions

    REAL(wp) :: r_Zs              !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value pore pressure

    REAL(wp) :: r_sp_heat_c

    REAL(wp) :: r_inv_rhom, r_inv_rho_c
    r_Zs = 0.0_wp
    r_exc_pore_pres = 0.0_wp

    r_h = qp(idx_h)

    IF ( r_h .GT. EPSILON(1.0_wp) ) THEN

       r_hu = qp(idx_hu)
       r_hv = qp(idx_hv)

       r_u = qp(idx_u)
       r_v = qp(idx_v)

    ELSE

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    r_T  = qp(idx_T)

    r_xs = qp(idx_solid_first:idx_solid_last)
    r_xg = qp(idx_add_gas_first:idx_add_gas_last)

    IF (PRESENT(context)) THEN
       CALL enforce_mass_fraction_closure(r_xs, r_xg, context)
    ELSE
       CALL enforce_mass_fraction_closure(r_xs, r_xg)
    END IF

    CALL eval_mixture_properties_from_mass_fractions(                          &
         ! IN
         r_T, r_xs, r_xg,                                                &
         ! OUT
         r_rho_m, r_inv_rhom, r_rho_c, r_inv_rho_c,                           &
         r_alphas, r_alphag)

    CALL eval_mixture_heat_capacity(r_xs, r_xg, r_xc, r_sp_heat_c,      &
         r_sp_heat_mix)

    IF ( stoch_transport_flag) r_Zs = qp(idx_stoch)

    IF ( pore_pressure_flag ) r_exc_pore_pres = qp(idx_pore)

    qc(idx_totMassEqn) = r_rho_m * r_h

    qc(idx_uEqn) = r_rho_m * r_hu
    qc(idx_vEqn) = r_rho_m * r_hv

    ! Thermal energy per unit horizontal area. Kinetic energy is not part of
    ! the fourth conservative variable.
    qc(idx_engyEqn) = r_h * r_rho_m * r_sp_heat_mix * r_T

    qc(idx_solid_first:idx_solid_last) = r_xs * qc(idx_totMassEqn)
    qc(idx_add_gas_first:idx_add_gas_last) = r_xg * qc(idx_totMassEqn)

    IF ( stoch_transport_flag ) qc(idx_stoch) = r_Zs * qc(idx_totMassEqn)

    IF ( pore_pressure_flag ) qc(idx_pore) = r_exc_pore_pres * qc(idx_totMassEqn)


    RETURN

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Convert a physical state to the auxiliary free-surface reconstruction representation.
  !
  !> This subroutine evaluates from the physical local variables qpj, the three
  !> additional local variables qp2j = (h+B,u,v).
  !> @author
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !>
  !> \param[in] qpj Local physical state of length n_vars+2 with component mass fractions.
  !> \param[in] Bj Local bed elevation [m].
  !> \param[out] qp2j Three-component auxiliary vector [eta,u,v], with eta=h+B.
  !******************************************************************************

  SUBROUTINE qp_to_qp2(qpj,Bj,qp2j)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Bj
    REAL(wp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(idx_h) + Bj

    IF ( qpj(idx_h) .LE. 0.0_wp ) THEN

       qp2j(2) = 0.0_wp
       qp2j(3) = 0.0_wp

    ELSE

       qp2j(2) = qpj(idx_hu)/qpj(idx_h)
       qp2j(3) = qpj(idx_hv)/qpj(idx_h)

    END IF

    RETURN

  END SUBROUTINE qp_to_qp2

  !******************************************************************************
  !> \brief Attach optional state context to a fatal mass-fraction closure error.
  !> \param[in] message Error description.
  !> \param[in] context Optional cell, face, or stage context.
  !******************************************************************************
  SUBROUTINE mass_fraction_closure_error(message, context)

    CHARACTER(LEN=*), INTENT(IN) :: message
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: context

    IF (PRESENT(context)) THEN
       CALL fatal_error(TRIM(context)//': '//TRIM(message))
    ELSE
       CALL fatal_error(message)
    END IF

  END SUBROUTINE mass_fraction_closure_error

END MODULE state_conversion_2d
