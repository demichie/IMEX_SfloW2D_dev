!*******************************************************************************
!> Shared REAL/COMPLEX mixture kernels for the new state architecture.
!*******************************************************************************
MODULE state_new_mixture

    USE parameters_2d, ONLY : wp, gas_flag
    USE parameters_2d, ONLY : n_solid, n_add_gas
    USE constitutive_parameters_2d, ONLY : sp_heat_a, sp_heat_g, sp_heat_s, sp_heat_l
    USE constitutive_parameters_2d, ONLY : sp_gas_const_a, sp_gas_const_g
    USE constitutive_parameters_2d, ONLY : inv_pres, rho_l, inv_rho_l
    USE constitutive_parameters_2d, ONLY : inv_rho_s

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: evaluate_mixture

    INTERFACE evaluate_mixture
        MODULE PROCEDURE evaluate_mixture_real
        MODULE PROCEDURE evaluate_mixture_complex
    END INTERFACE evaluate_mixture

CONTAINS

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

END MODULE state_new_mixture
