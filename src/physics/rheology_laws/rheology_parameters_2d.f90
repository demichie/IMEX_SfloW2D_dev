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

    ! Canonical input configuration for each material layer. Every rheology
    ! coefficient is declared once here, not once per NAMELIST.
    PUBLIC :: rheology_config_type
    PUBLIC :: rheology_l1, rheology_l2
    TYPE :: rheology_config_type
        CHARACTER(LEN=32) :: rheology_model='none'
        CHARACTER(LEN=16) :: rheology_method=''
        INTEGER :: rheology_group=1
        REAL(wp) :: mu=-1.0_wp, xi=-1.0_wp, tau=-1.0_wp
        REAL(wp) :: nu_ref=-1.0_wp, visc_par=-1.0_wp
        REAL(wp) :: T_ref=-1.0_wp, tau0=-1.0_wp
        REAL(wp) :: alpha2=-1.0_wp, beta2=-1.0_wp
        REAL(wp) :: alpha1_ref=-1.0_wp, alpha1_coeff=-1.0_wp, beta1=-1.0_wp
        REAL(wp) :: Kappa=-1.0_wp, n_td=-1.0_wp
        REAL(wp) :: friction_factor=-1.0_wp
        REAL(wp) :: mu_0=-1.0_wp, mu_inf=-1.0_wp
        REAL(wp) :: Fr_0=-1.0_wp, U_w=-1.0_wp
        REAL(wp) :: mu_s=-1.0_wp, mu_2=-1.0_wp
        REAL(wp) :: I_0=-1.0_wp, muI_inf=-1.0_wp
    END TYPE rheology_config_type

    ! Canonical layer configurations shared between input initialization and
    ! the simulation context. The parser owns their initialization and validation.
    TYPE(rheology_config_type) :: rheology_l1, rheology_l2





END MODULE rheology_parameters_2d
