!> @brief Mixture thermodynamic properties shared by REAL and COMPLEX state paths.
!*******************************************************************************
!> Shared REAL/COMPLEX mixture kernels for the new state architecture.
!*******************************************************************************
MODULE state_mixture

    USE parameters_2d, ONLY : wp, gas_flag
    USE parameters_2d, ONLY : n_solid, n_add_gas, n_vars, idx_h, idx_T,                                                &
        idx_solid_first, idx_solid_last, idx_add_gas_first, idx_add_gas_last
    USE constitutive_parameters_2d, ONLY : sp_heat_a, sp_heat_g, sp_heat_l
    USE solid_transport_parameters_2d, ONLY : sp_heat_s
    USE constitutive_parameters_2d, ONLY : sp_gas_const_a, sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : inv_pres, pres, rho_l, inv_rho_l
    USE constitutive_parameters_2d, ONLY : rho_a_amb, grav
    USE solid_transport_parameters_2d, ONLY : inv_rho_s, diam_s, rho_s
    USE state_composition, ONLY : normalize_face_composition

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: evaluate_mixture
    PUBLIC :: evaluate_mixture_heat_capacity
    PUBLIC :: evaluate_mixture_density
    PUBLIC :: evaluate_reduced_gravity
    PUBLIC :: evaluate_gas_specific_volume
    PUBLIC :: evaluate_mixture_phase_properties
    PUBLIC :: sauter_diameter, average_density_solids
    PUBLIC :: evaluate_mixture_from_volume, primitive_to_volume_fractions

    INTERFACE evaluate_mixture
        MODULE PROCEDURE evaluate_mixture_real
        MODULE PROCEDURE evaluate_mixture_complex
    END INTERFACE evaluate_mixture

    !> Complete density and phase-volume properties when heat capacity is
    !> already available or is not consumed. Shares the exact density formula
    !> with evaluate_mixture, without repeating heat-capacity work.
    INTERFACE evaluate_mixture_phase_properties
        MODULE PROCEDURE evaluate_mixture_phase_properties_real
        MODULE PROCEDURE evaluate_mixture_phase_properties_complex
    END INTERFACE evaluate_mixture_phase_properties

    INTERFACE evaluate_mixture_heat_capacity
        MODULE PROCEDURE evaluate_mixture_heat_capacity_real
        MODULE PROCEDURE evaluate_mixture_heat_capacity_complex
    END INTERFACE evaluate_mixture_heat_capacity

    INTERFACE evaluate_mixture_density
        MODULE PROCEDURE evaluate_mixture_density_real
        MODULE PROCEDURE evaluate_mixture_density_complex
    END INTERFACE evaluate_mixture_density

    INTERFACE evaluate_gas_specific_volume
        MODULE PROCEDURE gas_specific_volume_real
        MODULE PROCEDURE gas_specific_volume_complex
        MODULE PROCEDURE gas_specific_volume_complex_real_const
    END INTERFACE evaluate_gas_specific_volume

    INTERFACE evaluate_reduced_gravity
        MODULE PROCEDURE evaluate_reduced_gravity_real
        MODULE PROCEDURE evaluate_reduced_gravity_complex
    END INTERFACE evaluate_reduced_gravity

CONTAINS

    !> @brief Evaluate all REAL mixture thermodynamic and phase-volume properties.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @param[out] rho_c Carrier-phase density (kg/m3).
    !> @param[out] inv_rho_c Reciprocal carrier density (m3/kg).
    !> @param[out] cp_m Mass-specific mixture heat capacity (J/kg/K).
    !> @param[out] cp_c Carrier heat capacity (J/kg/K).
    !> @param[out] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    SUBROUTINE evaluate_mixture_real(T, xs, xg, rho_m, inv_rho_m, rho_c, inv_rho_c,                                    &
            cp_m, cp_c, alphas, alphag)
        REAL(wp), INTENT(IN) :: T
        REAL(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m, inv_rho_m, rho_c, inv_rho_c
        REAL(wp), INTENT(OUT) :: cp_m, cp_c
        REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: xc, gas_const_c
        REAL(wp) :: inv_rho_g(n_add_gas)
        REAL(wp) :: carrier_sp_heat_mass

        INCLUDE 'inc/mixture_heat_capacity.inc'
        INCLUDE 'inc/mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_real

    !> @brief Evaluate all COMPLEX mixture thermodynamic and phase-volume properties.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @param[out] rho_c Carrier-phase density (kg/m3).
    !> @param[out] inv_rho_c Reciprocal carrier density (m3/kg).
    !> @param[out] cp_m Mass-specific mixture heat capacity (J/kg/K).
    !> @param[out] cp_c Carrier heat capacity (J/kg/K).
    !> @param[out] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    SUBROUTINE evaluate_mixture_complex(T, xs, xg, rho_m, inv_rho_m, rho_c, inv_rho_c,                                 &
            cp_m, cp_c, alphas, alphag)
        COMPLEX(wp), INTENT(IN) :: T
        COMPLEX(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m, inv_rho_m, rho_c, inv_rho_c
        COMPLEX(wp), INTENT(OUT) :: cp_m, cp_c
        COMPLEX(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
        COMPLEX(wp) :: xc, gas_const_c
        COMPLEX(wp) :: inv_rho_g(n_add_gas)
        COMPLEX(wp) :: carrier_sp_heat_mass

        INCLUDE 'inc/mixture_heat_capacity.inc'
        INCLUDE 'inc/mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_complex

    !> @brief Evaluate REAL bulk density from composition and temperature.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    SUBROUTINE evaluate_mixture_density_real(T,xs,xg,rho_m,inv_rho_m)
        REAL(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m, inv_rho_m
        REAL(wp) :: xc, gas_const_c, inv_rho_c
        INCLUDE 'inc/density_from_mass.inc'
    END SUBROUTINE evaluate_mixture_density_real

    !> @brief Evaluate COMPLEX bulk density from composition and temperature.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    SUBROUTINE evaluate_mixture_density_complex(T,xs,xg,rho_m,inv_rho_m)
        COMPLEX(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m, inv_rho_m
        COMPLEX(wp) :: xc, gas_const_c, inv_rho_c
        INCLUDE 'inc/density_from_mass.inc'
    END SUBROUTINE evaluate_mixture_density_complex

    !> @brief Evaluate REAL mixture and carrier heat capacities.
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] cp_m Mass-specific mixture heat capacity (J/kg/K).
    !> @param[out] cp_c Carrier heat capacity (J/kg/K).
    SUBROUTINE evaluate_mixture_heat_capacity_real(xs,xg,cp_m,cp_c)
        REAL(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: cp_m, cp_c
        REAL(wp) :: xc, carrier_sp_heat_mass
        INCLUDE 'inc/mixture_heat_capacity.inc'
    END SUBROUTINE evaluate_mixture_heat_capacity_real

    !> @brief Evaluate COMPLEX mixture and carrier heat capacities.
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] cp_m Mass-specific mixture heat capacity (J/kg/K).
    !> @param[out] cp_c Carrier heat capacity (J/kg/K).
    SUBROUTINE evaluate_mixture_heat_capacity_complex(xs,xg,cp_m,cp_c)
        COMPLEX(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: cp_m, cp_c
        COMPLEX(wp) :: xc, carrier_sp_heat_mass
        INCLUDE 'inc/mixture_heat_capacity.inc'
    END SUBROUTINE evaluate_mixture_heat_capacity_complex

    !> @brief Evaluate REAL mixture and phase-volume properties.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @param[out] rho_c Carrier-phase density (kg/m3).
    !> @param[out] inv_rho_c Reciprocal carrier density (m3/kg).
    !> @param[out] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    SUBROUTINE evaluate_mixture_phase_properties_real(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,alphas,alphag)
        REAL(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m,inv_rho_m,rho_c,inv_rho_c
        REAL(wp), INTENT(OUT) :: alphas(n_solid),alphag(n_add_gas)
        REAL(wp) :: xc,gas_const_c,inv_rho_g(n_add_gas)

        INCLUDE 'inc/mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_phase_properties_real

    !> @brief Evaluate COMPLEX mixture and phase-volume properties.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[in] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @param[out] rho_c Carrier-phase density (kg/m3).
    !> @param[out] inv_rho_c Reciprocal carrier density (m3/kg).
    !> @param[out] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    SUBROUTINE evaluate_mixture_phase_properties_complex(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,alphas,alphag)
        COMPLEX(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m,inv_rho_m,rho_c,inv_rho_c
        COMPLEX(wp), INTENT(OUT) :: alphas(n_solid),alphag(n_add_gas)
        COMPLEX(wp) :: xc,gas_const_c,inv_rho_g(n_add_gas)

        INCLUDE 'inc/mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_phase_properties_complex

    !> @brief Convert prescribed phase volumes to mass-based mixture properties.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    !> @param[inout] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] rho_m Bulk mixture density (kg/m3).
    !> @param[out] inv_rhom Reciprocal mixture density (m3/kg).
    !> @param[out] rho_c Carrier-phase density (kg/m3).
    !> @param[out] xs Solid mass fractions, with n_solid entries; may be clipped or normalized.
    !> @param[out] xg Additional-gas mass fractions, with n_add_gas entries; may be normalized.
    !> @param[out] xc Total carrier mass fraction, including additional gases.
    !> @param[out] cp_c Carrier heat capacity (J/kg/K).
    !> @param[out] cp_m Mass-specific mixture heat capacity (J/kg/K).
    SUBROUTINE evaluate_mixture_from_volume(T,alphag,alphas,rho_m,inv_rhom,                                            &
            rho_c,xs,xg,xc,cp_c,cp_m)
        REAL(wp), INTENT(IN) :: T,alphag(n_add_gas)
        REAL(wp), INTENT(INOUT) :: alphas(n_solid)
        REAL(wp), INTENT(OUT) :: rho_m,inv_rhom,rho_c
        REAL(wp), INTENT(OUT) :: xs(n_solid),xg(n_add_gas),xc,cp_c,cp_m
        REAL(wp) :: alphas_tot,alphac,alpha_air,rho_a,rho_g(n_add_gas)
        REAL(wp) :: rho_from_volume,inv_rho_c,phase_alphas(n_solid),phase_alphag(n_add_gas)

        ! Only the volume-to-mass conversion belongs here. All resulting
        ! thermodynamic properties come from the canonical mass-based kernel.
        alphas_tot=SUM(alphas)

        IF (alphas_tot>1.0_wp) THEN
            alphas=alphas/alphas_tot
        ELSEIF (alphas_tot<0.0_wp) THEN
            alphas=0.0_wp
        END IF

        alphac=1.0_wp-SUM(alphas)

        IF (gas_flag) THEN
            rho_a=1.0_wp/evaluate_gas_specific_volume(T,sp_gas_const_a)
            rho_g=1.0_wp/evaluate_gas_specific_volume(T,sp_gas_const_g)
            alpha_air=alphac-SUM(alphag)
            rho_from_volume=SUM(alphas*rho_s)+alpha_air*rho_a+SUM(alphag*rho_g)
        ELSE
            rho_g=0.0_wp
            rho_from_volume=SUM(alphas*rho_s)+alphac*rho_l
        END IF

        xs=alphas*rho_s/rho_from_volume
        xg=alphag*rho_g/rho_from_volume

        CALL evaluate_mixture(T,xs,xg,rho_m,inv_rhom,rho_c,inv_rho_c,                                                  &
            cp_m,cp_c,phase_alphas,phase_alphag)
        xc=1.0_wp-SUM(xs)
    END SUBROUTINE evaluate_mixture_from_volume

    !> @brief Convert primitive mass fractions to phase-volume fractions.
    !> @param[in] qp Primitive cell vector: thickness, velocity, temperature and fractions.
    !> @param[out] alphas Solid phase volume fractions, with n_solid entries.
    !> @param[out] alphag Additional-gas phase volume fractions, with n_add_gas entries.
    SUBROUTINE primitive_to_volume_fractions(qp,alphas,alphag)
        REAL(wp), INTENT(IN) :: qp(n_vars)
        REAL(wp), INTENT(OUT) :: alphas(n_solid),alphag(n_add_gas)
        REAL(wp) :: xs(n_solid),xg(n_add_gas)
        REAL(wp) :: rho_m,inv_rhom,rho_c,inv_rho_c

        alphas=0.0_wp
        alphag=0.0_wp
        IF (qp(idx_h)<=EPSILON(1.0_wp)) RETURN
        xs=qp(idx_solid_first:idx_solid_last)
        xg=qp(idx_add_gas_first:idx_add_gas_last)
        IF (.NOT.normalize_face_composition(xs,xg))                                                                    &
            ERROR STOP 'Invalid primitive composition in volume conversion'
        CALL evaluate_mixture_phase_properties(qp(idx_T),xs,xg,rho_m,inv_rhom,rho_c,inv_rho_c,                         &
            alphas,alphag)
    END SUBROUTINE primitive_to_volume_fractions

    !> @brief Evaluate ideal-gas specific volume for REAL temperature.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] gas_const Specific ideal-gas constant (J/kg/K).
    !> @return Gas specific volume (m3/kg).
    PURE ELEMENTAL REAL(wp) FUNCTION gas_specific_volume_real(T,gas_const) RESULT(inv_rho)
        REAL(wp), INTENT(IN) :: T,gas_const
        INCLUDE 'inc/gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_real

    !> @brief Evaluate ideal-gas specific volume for COMPLEX inputs.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] gas_const Specific ideal-gas constant (J/kg/K).
    !> @return Gas specific volume (m3/kg).
    PURE ELEMENTAL COMPLEX(wp) FUNCTION gas_specific_volume_complex(T,gas_const) RESULT(inv_rho)
        COMPLEX(wp), INTENT(IN) :: T,gas_const
        INCLUDE 'inc/gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_complex

    !> @brief Evaluate COMPLEX ideal-gas specific volume with a REAL gas constant.
    !> @param[in] T Mixture temperature (K).
    !> @param[in] gas_const Specific ideal-gas constant (J/kg/K).
    !> @return Gas specific volume (m3/kg).
    PURE ELEMENTAL COMPLEX(wp) FUNCTION gas_specific_volume_complex_real_const(T,gas_const) RESULT(inv_rho)
        COMPLEX(wp), INTENT(IN) :: T
        REAL(wp), INTENT(IN) :: gas_const
        INCLUDE 'inc/gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_complex_real_const

    !> @brief Evaluate REAL reduced gravity relative to ambient air.
    !> @param[in] rho_m Bulk mixture density (kg/m3).
    !> @param[in] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @return Reduced gravity (m/s2).
    PURE REAL(wp) FUNCTION evaluate_reduced_gravity_real(rho_m,inv_rho_m) RESULT(red_grav)
        REAL(wp), INTENT(IN) :: rho_m,inv_rho_m
        INCLUDE 'inc/reduced_gravity.inc'
    END FUNCTION evaluate_reduced_gravity_real

    !> @brief Evaluate COMPLEX reduced gravity relative to ambient air.
    !> @param[in] rho_m Bulk mixture density (kg/m3).
    !> @param[in] inv_rho_m Reciprocal bulk mixture density (m3/kg).
    !> @return Reduced gravity (m/s2).
    PURE COMPLEX(wp) FUNCTION evaluate_reduced_gravity_complex(rho_m,inv_rho_m) RESULT(red_grav)
        COMPLEX(wp), INTENT(IN) :: rho_m,inv_rho_m
        INCLUDE 'inc/reduced_gravity.inc'
    END FUNCTION evaluate_reduced_gravity_complex

    !> @brief Evaluate the Sauter mean diameter from solid volume fractions.
    !> @param[in] alpha_solids Solid volume fractions used as weights over size or density classes.
    !> @return Effective solid particle diameter (m).
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

    !> @brief Evaluate the volume-weighted density of solid phases.
    !> @param[in] alpha_solids Solid volume fractions used as weights over size or density classes.
    !> @return Effective solid density (kg/m3).
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

END MODULE state_mixture
