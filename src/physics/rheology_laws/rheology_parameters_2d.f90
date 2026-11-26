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
    !> @brief Canonical coefficients and scheduling for one material layer's rheology.
    TYPE :: rheology_config_type
        LOGICAL :: stochastic_rheology_flag=.FALSE. !< Stochastic rheology activation for this layer.
        !> Selected basal resistance law, or 'none' when disabled.
        CHARACTER(LEN=32) :: rheology_model='none'
        !> Configured explicit, Newton or analytic source treatment.
        CHARACTER(LEN=16) :: rheology_method=''
        !> Ordered Lie group containing this layer's rheological source.
        INTEGER :: rheology_group=1
        !> Friction coefficient, turbulent resistance coefficient and yield stress by law.
        REAL(wp) :: mu=-1.0_wp, xi=-1.0_wp, tau=-1.0_wp
        !> Reference viscosity and viscosity-shape coefficient of configured law.
        REAL(wp) :: nu_ref=-1.0_wp, visc_par=-1.0_wp
        !> Reference temperature [K] and reference yield stress of viscous laws.
        REAL(wp) :: T_ref=-1.0_wp, tau0=-1.0_wp
        !> Second-branch empirical coefficients for the configured rheological law.
        REAL(wp) :: alpha2=-1.0_wp, beta2=-1.0_wp
        !> First-branch reference/coefficient and exponent for empirical drag.
        REAL(wp) :: alpha1_ref=-1.0_wp, alpha1_coeff=-1.0_wp, beta1=-1.0_wp
        !> Empirical law parameter and power exponent used by the selected rheology.
        REAL(wp) :: Kappa=-1.0_wp, n_td=-1.0_wp
        !> Dimensionless coefficient of quadratic resistance.
        REAL(wp) :: friction_factor=-1.0_wp
        !> Low- and high-speed friction coefficients of the Froude-based law.
        REAL(wp) :: mu_0=-1.0_wp, mu_inf=-1.0_wp
        !> Froude transition parameter and reference velocity [m/s] of the law.
        REAL(wp) :: Fr_0=-1.0_wp, U_w=-1.0_wp
        !> Static and dynamic friction limits of the mu(I) relation.
        REAL(wp) :: mu_s=-1.0_wp, mu_2=-1.0_wp
        !> Inertial-number transition parameter and asymptotic mu(I) coefficient.
        REAL(wp) :: I_0=-1.0_wp, muI_inf=-1.0_wp
    END TYPE rheology_config_type

    ! Canonical layer configurations shared between input initialization and
    ! the simulation context. The parser owns their initialization and validation.
    TYPE(rheology_config_type), TARGET :: rheology_l1
    TYPE(rheology_config_type), TARGET :: rheology_l2

END MODULE rheology_parameters_2d
