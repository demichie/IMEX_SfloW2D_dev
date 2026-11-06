!> @brief Canonical state API facade, exposing the physical and numerical representations.
!*******************************************************************************
!> Public entry point for the new state architecture.
!*******************************************************************************
MODULE state
    USE state_types
    USE state_fields
    USE state_mixture
    USE state_conversions
    USE state_physical
    USE state_composition
    IMPLICIT NONE
    PUBLIC
END MODULE state
