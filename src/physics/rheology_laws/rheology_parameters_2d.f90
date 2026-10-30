!********************************************************************************
!> \brief Shared input storage for rheology coefficients used by multiple laws.
!>
!> This module contains no rheology law terms. Each law has its own module;
!> these values retain the existing shared namelist names.
!********************************************************************************

MODULE rheology_parameters_2d

    USE parameters_2d, ONLY : wp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: I_0, mu_2, mu_s, muI_inf
    PUBLIC :: mu_0, mu_inf, Fr_0, U_w
    PUBLIC :: xi, T_ref

    REAL(wp) :: I_0
    REAL(wp) :: mu_2
    REAL(wp) :: mu_s
    REAL(wp) :: muI_inf

    REAL(wp) :: mu_0
    REAL(wp) :: mu_inf
    REAL(wp) :: Fr_0
    REAL(wp) :: U_w

    REAL(wp) :: xi
    REAL(wp) :: T_ref

END MODULE rheology_parameters_2d
