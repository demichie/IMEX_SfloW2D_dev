!> @brief Canonical physical options for paired L1/L2 mass transfers.
!> @details The interlayer exchange belongs to the layer pair, never to either layer.
MODULE interlayer_exchange_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: interlayer_exchange_config_type, interlayer_exchange

    TYPE :: interlayer_exchange_config_type
        LOGICAL :: settling_flag=.FALSE.                 !< Transfer solids from L2 to L1.
        LOGICAL :: elutriation_flag=.FALSE.              !< Transfer solids from L1 to L2.
        REAL(wp) :: elutriation_solid_volume_fraction=-1.0_wp !< Solid fraction in lifted parcel.
        CHARACTER(LEN=24) :: elutriation_carrier_mode='LIMIT_SOLIDS'
        CHARACTER(LEN=16) :: method='EXPLICIT'           !< Paired mass-transfer Lie method.
        INTEGER :: group=4                              !< Paired mass-transfer Lie group.
    END TYPE interlayer_exchange_config_type

    TYPE(interlayer_exchange_config_type) :: interlayer_exchange
END MODULE interlayer_exchange_parameters_2d
