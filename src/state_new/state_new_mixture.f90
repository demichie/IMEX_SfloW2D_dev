!*******************************************************************************
!> Shared REAL/COMPLEX mixture kernels for the new state architecture.
!*******************************************************************************
MODULE state_new_mixture

    USE parameters_2d, ONLY : wp, gas_flag
    USE parameters_2d, ONLY : n_solid, n_add_gas
    USE constitutive_parameters_2d, ONLY : sp_heat_a, sp_heat_g, sp_heat_s, sp_heat_l
    USE constitutive_parameters_2d, ONLY : sp_gas_const_a, sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : inv_pres, rho_l, inv_rho_l
    USE constitutive_parameters_2d, ONLY : inv_rho_s, diam_s, rho_s

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: evaluate_mixture
    PUBLIC :: evaluate_mixture_heat_capacity
    PUBLIC :: evaluate_mixture_density
    PUBLIC :: sauter_diameter, average_density_solids

    INTERFACE evaluate_mixture
        MODULE PROCEDURE evaluate_mixture_real
        MODULE PROCEDURE evaluate_mixture_complex
    END INTERFACE evaluate_mixture

    INTERFACE evaluate_mixture_heat_capacity
        MODULE PROCEDURE evaluate_mixture_heat_capacity_real
        MODULE PROCEDURE evaluate_mixture_heat_capacity_complex
    END INTERFACE evaluate_mixture_heat_capacity

    INTERFACE evaluate_mixture_density
        MODULE PROCEDURE evaluate_mixture_density_real
        MODULE PROCEDURE evaluate_mixture_density_complex
    END INTERFACE evaluate_mixture_density

CONTAINS

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

END MODULE state_new_mixture
