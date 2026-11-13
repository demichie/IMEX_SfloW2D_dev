!> \brief Canonical per-layer thermal configuration and derived cooling coefficients.
!> \details Ground temperature belongs to the single shared substrate; atmosphere stays shared.
MODULE thermal_parameters_2d
    USE parameters_2d, ONLY : wp
    USE constitutive_parameters_2d, ONLY : SBconst
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: thermal_config_type, thermal_l1, thermal_l2
    PUBLIC :: thermal_radiative_coefficient, thermal_convective_coefficient

    TYPE :: thermal_config_type
        LOGICAL :: thermal_cooling_flag=.FALSE.
        CHARACTER(LEN=16) :: thermal_method=''
        INTEGER :: thermal_group=2
        REAL(wp) :: emissivity=-1.0_wp
        REAL(wp) :: atm_heat_transf_coeff=-1.0_wp
        REAL(wp) :: exp_area_fract=-1.0_wp
        REAL(wp) :: thermal_conductivity=-1.0_wp
        REAL(wp) :: c_p=-1.0_wp
        REAL(wp) :: enne=-1.0_wp
        REAL(wp) :: emme=-1.0_wp
    END TYPE thermal_config_type

    TYPE(thermal_config_type), TARGET :: thermal_l1
    TYPE(thermal_config_type) :: thermal_l2
CONTAINS
    !> \brief Reproduce the radiative coefficient without a cached global copy.
    PURE REAL(wp) FUNCTION thermal_radiative_coefficient(config) RESULT(coefficient)
        TYPE(thermal_config_type), INTENT(IN) :: config
        coefficient=SBconst*config%emissivity*config%exp_area_fract
    END FUNCTION thermal_radiative_coefficient

    !> \brief Reproduce the convective coefficient without a cached global copy.
    PURE REAL(wp) FUNCTION thermal_convective_coefficient(config) RESULT(coefficient)
        TYPE(thermal_config_type), INTENT(IN) :: config
        coefficient=config%atm_heat_transf_coeff*config%exp_area_fract
    END FUNCTION thermal_convective_coefficient
END MODULE thermal_parameters_2d
