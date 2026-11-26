!> @brief Canonical layer-local configurations for the combined mass-exchange source.
!> @details Each layer owns closure coefficients and a single Lie METHOD/GROUP.
!>          Shared material, terrain and atmospheric properties are not stored here.
MODULE mass_exchange_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: mass_exchange_config_type, mass_exchange_l1, mass_exchange_l2

    !> @brief Empirical bed erosion closure for one layer.
    TYPE :: mass_erosion_config_type
        REAL(wp) :: coefficient=-1.0_wp            !< Erosion coefficient [m-1].
    END TYPE mass_erosion_config_type

    !> @brief Sediment deposition and optional collective-drag closure.
    TYPE :: mass_deposition_config_type
        LOGICAL :: settling_flag=.FALSE.          !< Enable particulate deposition.
        LOGICAL :: collective_settling_flag=.FALSE. !< Enable collective drag correction.
        REAL(wp) :: A_drag=-1.0_wp                !< Collective drag coefficient.
        REAL(wp) :: B_drag=-1.0_wp                !< Collective drag logarithmic offset.
    END TYPE mass_deposition_config_type

    !> @brief Ambient-air entrainment activation for one layer.
    TYPE :: mass_entrainment_config_type
        LOGICAL :: flag=.FALSE.                   !< Enable atmospheric entrainment.
    END TYPE mass_entrainment_config_type

    !> @brief Prescribed carrier-liquid withdrawal for one layer.
    TYPE :: mass_carrier_loss_config_type
        REAL(wp) :: loss_rate=-1.0_wp             !< Carrier withdrawal speed [m s-1].
    END TYPE mass_carrier_loss_config_type

    !> @brief Steam production from deposited hot particles.
    TYPE :: mass_vaporization_config_type
        LOGICAL :: liquid_vaporization_flag=.FALSE. !< Enable steam generation.
        REAL(wp) :: gamma_steam=-1.0_wp           !< Heat fraction transferred to steam.
    END TYPE mass_vaporization_config_type

    !> @brief Temperature/viscosity-dependent immobilization of mobile lava.
    TYPE :: mass_lava_deposition_config_type
        LOGICAL :: flag=.FALSE.                   !< Enable lava immobilization.
        CHARACTER(LEN=16) :: model='CONSTANT'     !< CONSTANT, TEMPERATURE or VISCOSITY.
        REAL(wp) :: rate=0.0_wp                   !< Prescribed immobilization speed [m s-1].
        REAL(wp) :: visc_ref=1.0E-2_wp            !< Reference viscosity [m2 s-1].
        REAL(wp) :: exponent=1.0_wp               !< Viscosity model exponent.
        REAL(wp) :: temp_ref=900.0_wp             !< Reference lava temperature [K].
        REAL(wp) :: temp_width=100.0_wp           !< Transition temperature width [K].
    END TYPE mass_lava_deposition_config_type

    !> @brief All active external mass-exchange laws share one Lie scheduling choice.
    TYPE :: mass_exchange_config_type
        CHARACTER(LEN=16) :: method='EXPLICIT'    !< Combined source integration method.
        INTEGER :: group=4                        !< Combined source Lie group.
        TYPE(mass_erosion_config_type) :: erosion  !< Substrate erosion.
        TYPE(mass_deposition_config_type) :: deposition !< Solid settling and drag.
        TYPE(mass_entrainment_config_type) :: entrainment !< Atmospheric entrainment.
        TYPE(mass_carrier_loss_config_type) :: carrier_liquid_loss !< Liquid withdrawal.
        TYPE(mass_vaporization_config_type) :: vaporization !< Steam formation.
        TYPE(mass_lava_deposition_config_type) :: lava_deposition !< Lava immobilization.
    END TYPE mass_exchange_config_type

    !> Canonical per-layer owners; L2 is configured but cannot execute yet.
    TYPE(mass_exchange_config_type), TARGET :: mass_exchange_l1
    TYPE(mass_exchange_config_type), TARGET :: mass_exchange_l2
END MODULE mass_exchange_parameters_2d
