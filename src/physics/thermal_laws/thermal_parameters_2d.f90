!> \brief Canonical per-layer thermal configuration and derived cooling coefficients.
!> \details Ground temperature belongs to the single shared substrate; atmosphere stays shared.
MODULE thermal_parameters_2d
    USE parameters_2d, ONLY : wp
    USE constitutive_parameters_2d, ONLY : SBconst
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: thermal_config_type, thermal_l1, thermal_l2
    PUBLIC :: thermal_radiative_coefficient, thermal_convective_coefficient

    !> @brief Layer-local thermal cooling law and source scheduling parameters.
    TYPE :: thermal_config_type
        !> Enable the selected layer's thermal cooling source.
        LOGICAL :: thermal_cooling_flag=.FALSE.
        !> Configured source treatment for the thermal cooling term.
        CHARACTER(LEN=16) :: thermal_method=''
        !> Ordered Lie group identifier for thermal cooling.
        INTEGER :: thermal_group=2
        !> Surface emissivity for Stefan-Boltzmann radiation (dimensionless).
        REAL(wp) :: emissivity=-1.0_wp
        !> Atmosphere-to-surface convective heat-transfer coefficient.
        REAL(wp) :: atm_heat_transf_coeff=-1.0_wp
        !> Fraction of layer surface effectively exposed to cooling.
        REAL(wp) :: exp_area_fract=-1.0_wp
        !> Configured thermal conductivity used by the cooling law.
        REAL(wp) :: thermal_conductivity=-1.0_wp
        !> Configured heat capacity of the thermal model [J/(kg K)].
        REAL(wp) :: c_p=-1.0_wp
        !> Empirical cooling-law exponent n.
        REAL(wp) :: enne=-1.0_wp
        !> Empirical cooling-law exponent m.
        REAL(wp) :: emme=-1.0_wp
    END TYPE thermal_config_type

    TYPE(thermal_config_type), TARGET :: thermal_l1
    TYPE(thermal_config_type), TARGET :: thermal_l2
CONTAINS
    !> \brief Reproduce the radiative coefficient without a cached global copy.
    !> @param[in] config Canonical thermal settings containing emissivity and exposure.
    !> @return Stefan-Boltzmann radiation coefficient scaled by emissivity/exposure.
    PURE REAL(wp) FUNCTION thermal_radiative_coefficient(config) RESULT(coefficient)
        TYPE(thermal_config_type), INTENT(IN) :: config
        coefficient=SBconst*config%emissivity*config%exp_area_fract
    END FUNCTION thermal_radiative_coefficient

    !> \brief Reproduce the convective coefficient without a cached global copy.
    !> @param[in] config Canonical thermal settings containing emissivity and exposure.
    !> @return Stefan-Boltzmann radiation coefficient scaled by emissivity/exposure.
    PURE REAL(wp) FUNCTION thermal_convective_coefficient(config) RESULT(coefficient)
        TYPE(thermal_config_type), INTENT(IN) :: config
        coefficient=config%atm_heat_transf_coeff*config%exp_area_fract
    END FUNCTION thermal_convective_coefficient
END MODULE thermal_parameters_2d
