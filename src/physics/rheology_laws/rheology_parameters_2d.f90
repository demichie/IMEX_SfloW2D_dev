!********************************************************************************
!> \brief Canonical per-layer rheology configuration shared by the physical laws.
!>
!> Every rheology law receives an explicit read-only configuration argument.
!> The type contains the sole coefficient storage; no per-law global copies exist.
!********************************************************************************

MODULE rheology_parameters_2d

    USE parameters_2d, ONLY : wp

    IMPLICIT NONE

    PRIVATE

    ! One coefficient representation per layer, independent of input naming.
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
    TYPE(rheology_config_type), TARGET :: rheology_l1
    TYPE(rheology_config_type), TARGET :: rheology_l2

END MODULE rheology_parameters_2d
