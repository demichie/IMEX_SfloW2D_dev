!> \brief Canonical per-layer thermal configuration and derived cooling coefficients.
!> \details Layer configurations also contain the common ground temperature.
!>          Both configurations must agree when two layers are enabled; atmosphere stays shared.
MODULE thermal_parameters_2d
    USE parameters_2d, ONLY : wp
    USE constitutive_parameters_2d, ONLY : SBconst
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: thermal_config_type, thermal_l1, thermal_l2
    PUBLIC :: thermal_radiative_coefficient, thermal_convective_coefficient
    PUBLIC :: thermal_ground_temperatures_match

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
        !> Ground temperature (K), equal between L1 and L2 for two-layer runs.
        REAL(wp) :: T_ground=-1.0_wp
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
    !> \brief Test the shared ground-temperature invariant of a two-layer configuration.
    !> \details Exact equality is required: the two namelists describe one physical terrain.
    PURE LOGICAL FUNCTION thermal_ground_temperatures_match(config_l1,config_l2) RESULT(matches)
        TYPE(thermal_config_type), INTENT(IN) :: config_l1,config_l2
        matches=config_l1%T_ground==config_l2%T_ground
    END FUNCTION thermal_ground_temperatures_match
END MODULE thermal_parameters_2d
