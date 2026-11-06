!*******************************************************************************
!> Shared REAL/COMPLEX mixture kernels for the new state architecture.
!*******************************************************************************
MODULE state_mixture

    USE parameters_2d, ONLY : wp, gas_flag
    USE parameters_2d, ONLY : n_solid, n_add_gas, n_vars, idx_h, idx_T, &
        idx_solid_first, idx_solid_last, idx_add_gas_first, idx_add_gas_last
    USE constitutive_parameters_2d, ONLY : sp_heat_a, sp_heat_g, sp_heat_s, sp_heat_l
    USE constitutive_parameters_2d, ONLY : sp_gas_const_a, sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : inv_pres, pres, rho_l, inv_rho_l
    USE constitutive_parameters_2d, ONLY : inv_rho_s, diam_s, rho_s, rho_a_amb, grav
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

    !> Canonical ideal-gas specific volume. The COMPLEX variants retain
    !> complex-step perturbations of temperature and mixture gas constant.
    PURE ELEMENTAL REAL(wp) FUNCTION gas_specific_volume_real(T,gas_const) RESULT(inv_rho)
        REAL(wp), INTENT(IN) :: T,gas_const
        INCLUDE 'gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_real

    PURE ELEMENTAL COMPLEX(wp) FUNCTION gas_specific_volume_complex(T,gas_const) RESULT(inv_rho)
        COMPLEX(wp), INTENT(IN) :: T,gas_const
        INCLUDE 'gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_complex

    PURE ELEMENTAL COMPLEX(wp) FUNCTION gas_specific_volume_complex_real_const(T,gas_const) RESULT(inv_rho)
        COMPLEX(wp), INTENT(IN) :: T
        REAL(wp), INTENT(IN) :: gas_const
        INCLUDE 'gas_specific_volume.inc'
    END FUNCTION gas_specific_volume_complex_real_const


    PURE REAL(wp) FUNCTION evaluate_reduced_gravity_real(rho_m,inv_rho_m) RESULT(red_grav)
        REAL(wp), INTENT(IN) :: rho_m,inv_rho_m
        INCLUDE 'reduced_gravity.inc'
    END FUNCTION evaluate_reduced_gravity_real

    PURE COMPLEX(wp) FUNCTION evaluate_reduced_gravity_complex(rho_m,inv_rho_m) RESULT(red_grav)
        COMPLEX(wp), INTENT(IN) :: rho_m,inv_rho_m
        INCLUDE 'reduced_gravity.inc'
    END FUNCTION evaluate_reduced_gravity_complex

    SUBROUTINE evaluate_mixture_density_real(T,xs,xg,rho_m,inv_rho_m)
        REAL(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m, inv_rho_m
        REAL(wp) :: xc, gas_const_c, inv_rho_c
        INCLUDE 'density_from_mass.inc'
    END SUBROUTINE evaluate_mixture_density_real

    SUBROUTINE evaluate_mixture_density_complex(T,xs,xg,rho_m,inv_rho_m)
        COMPLEX(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m, inv_rho_m
        COMPLEX(wp) :: xc, gas_const_c, inv_rho_c
        INCLUDE 'density_from_mass.inc'
    END SUBROUTINE evaluate_mixture_density_complex

    SUBROUTINE evaluate_mixture_phase_properties_real(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,alphas,alphag)
        REAL(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m,inv_rho_m,rho_c,inv_rho_c
        REAL(wp), INTENT(OUT) :: alphas(n_solid),alphag(n_add_gas)
        REAL(wp) :: xc,gas_const_c,inv_rho_g(n_add_gas)

        INCLUDE 'mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_phase_properties_real

    SUBROUTINE evaluate_mixture_phase_properties_complex(T,xs,xg,rho_m,inv_rho_m,rho_c,inv_rho_c,alphas,alphag)
        COMPLEX(wp), INTENT(IN) :: T, xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m,inv_rho_m,rho_c,inv_rho_c
        COMPLEX(wp), INTENT(OUT) :: alphas(n_solid),alphag(n_add_gas)
        COMPLEX(wp) :: xc,gas_const_c,inv_rho_g(n_add_gas)

        INCLUDE 'mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_phase_properties_complex

    SUBROUTINE evaluate_mixture_heat_capacity_real(xs,xg,cp_m,cp_c)
        REAL(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: cp_m, cp_c
        REAL(wp) :: xc, carrier_sp_heat_mass
        INCLUDE 'mixture_heat_capacity.inc'
    END SUBROUTINE evaluate_mixture_heat_capacity_real

    SUBROUTINE evaluate_mixture_heat_capacity_complex(xs,xg,cp_m,cp_c)
        COMPLEX(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: cp_m, cp_c
        COMPLEX(wp) :: xc, carrier_sp_heat_mass
        INCLUDE 'mixture_heat_capacity.inc'
    END SUBROUTINE evaluate_mixture_heat_capacity_complex

    SUBROUTINE evaluate_mixture_real(T, xs, xg, rho_m, inv_rho_m, rho_c, inv_rho_c, &
            cp_m, cp_c, alphas, alphag)
        REAL(wp), INTENT(IN) :: T
        REAL(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        REAL(wp), INTENT(OUT) :: rho_m, inv_rho_m, rho_c, inv_rho_c
        REAL(wp), INTENT(OUT) :: cp_m, cp_c
        REAL(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
        REAL(wp) :: xc, gas_const_c
        REAL(wp) :: inv_rho_g(n_add_gas)
        REAL(wp) :: carrier_sp_heat_mass

        INCLUDE 'mixture_heat_capacity.inc'
        INCLUDE 'mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_real

    SUBROUTINE evaluate_mixture_complex(T, xs, xg, rho_m, inv_rho_m, rho_c, inv_rho_c, &
            cp_m, cp_c, alphas, alphag)
        COMPLEX(wp), INTENT(IN) :: T
        COMPLEX(wp), INTENT(IN) :: xs(n_solid), xg(n_add_gas)
        COMPLEX(wp), INTENT(OUT) :: rho_m, inv_rho_m, rho_c, inv_rho_c
        COMPLEX(wp), INTENT(OUT) :: cp_m, cp_c
        COMPLEX(wp), INTENT(OUT) :: alphas(n_solid), alphag(n_add_gas)
        COMPLEX(wp) :: xc, gas_const_c
        COMPLEX(wp) :: inv_rho_g(n_add_gas)
        COMPLEX(wp) :: carrier_sp_heat_mass

        INCLUDE 'mixture_heat_capacity.inc'
        INCLUDE 'mixture_from_mass.inc'
    END SUBROUTINE evaluate_mixture_complex

    !> Convert prescribed volume fractions into the canonical source material.
    !> This inverse composition operation cannot be replaced by a mass-based
    !> density evaluation: the input fractions refer to constituent volumes.
    SUBROUTINE evaluate_mixture_from_volume(T,alphag,alphas,rho_m,inv_rhom, &
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

        CALL evaluate_mixture(T,xs,xg,rho_m,inv_rhom,rho_c,inv_rho_c, &
            cp_m,cp_c,phase_alphas,phase_alphag)
        xc=1.0_wp-SUM(xs)
    END SUBROUTINE evaluate_mixture_from_volume

    !> Primitive-state solid/gas volume fractions, zero for dry cells.
    !> The caller's state is never modified during composition closure.
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
        IF (.NOT.normalize_face_composition(xs,xg)) &
            ERROR STOP 'Invalid primitive composition in volume conversion'
        CALL evaluate_mixture_phase_properties(qp(idx_T),xs,xg,rho_m,inv_rhom,rho_c,inv_rho_c, &
            alphas,alphag)
    END SUBROUTINE primitive_to_volume_fractions

    !> Sauter diameter and volume-weighted particle density for granular laws.
    !> Keep the historical limiting behavior and arithmetic exactly unchanged.
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

END MODULE state_mixture
