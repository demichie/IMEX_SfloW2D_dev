!> @brief Canonical physical options for paired L1/L2 mass transfers.
!> @details The interlayer exchange belongs to the layer pair, never to either layer.
MODULE interlayer_mass_exchange_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: interlayer_mass_exchange_config_type, interlayer_exchange

    TYPE :: interlayer_mass_exchange_config_type
        LOGICAL :: settling_flag=.FALSE.                 !< Transfer solids from L2 to L1.
        LOGICAL :: elutriation_flag=.FALSE.              !< Transfer solids from L1 to L2.
        REAL(wp) :: elutriation_solid_volume_fraction=-1.0_wp !< Solid fraction in lifted parcel.
        CHARACTER(LEN=24) :: elutriation_carrier_mode='LIMIT_SOLIDS'
        CHARACTER(LEN=24) :: elutriation_law='LINEAR_USTAR'
        REAL(wp) :: elutriation_efficiency=1.0E-2_wp
        REAL(wp) :: shear_interlayer=1.0E-3_wp
        REAL(wp) :: shear_surface=1.0E-4_wp
        REAL(wp) :: upper_depth_gate=1.0E-3_wp       !< Smooth dry L2 gate [m].
        REAL(wp) :: rouse_critical=2.5_wp
        REAL(wp) :: rouse_exponent=1.5_wp
        REAL(wp) :: bed_availability_exponent=2.0_wp
        REAL(wp) :: sediment_diffusivity_ratio=1.0_wp
        REAL(wp) :: elutriation_donor_cap=0.5_wp    !< Maximum withdrawn per step.
        REAL(wp) :: shields_critical=-1.0_wp        !< Required for SHEAR_POWER.
        REAL(wp) :: shields_exponent=-1.0_wp        !< Required for SHEAR_POWER.
        REAL(wp) :: lift_length_diameters=-1.0_wp  !< Required for SHEAR_POWER.
        REAL(wp) :: settling_velocity_weight=1.0_wp !< Interfacial slip interpolation.
        REAL(wp) :: settling_plug_thickness=0.02_wp !< Thin-layer plug blend [m].
        CHARACTER(LEN=16) :: method='EXPLICIT'           !< Paired mass-transfer Lie method.
        INTEGER :: group=4                              !< Paired mass-transfer Lie group.
    END TYPE interlayer_mass_exchange_config_type

    TYPE(interlayer_mass_exchange_config_type), TARGET :: interlayer_exchange
END MODULE interlayer_mass_exchange_parameters_2d
